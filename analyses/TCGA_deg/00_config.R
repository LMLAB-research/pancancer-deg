#### Libraries ####

library(AnnotationDbi)
library(ashr)
library(BiocParallel)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(org.Hs.eg.db)
library(SummarizedExperiment)
library(TCGAbiolinks)

#### Analysis parameters ####

analysis_name <- "TCGA_deg"

# During development, limit scripts to the first N project IDs for faster runs.
# Set to Inf when running the full analysis.
debug_project_limit <- Inf

# Parallel execution for DESeq2 model fitting and log2FC shrinkage.
use_biocparallel <- TRUE
detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores)) {
  detected_cores <- 1
}
# Each worker holds its own copy of the data (about 2.4 GB for TCGA-BRCA), so
# the cap is set by memory, not cores: 4 workers stay within ~10 GB at peak.
max_biocparallel_workers <- 4
biocparallel_workers <- max(1, min(max_biocparallel_workers, detected_cores - 1))

# Step 04 skips designs already fitted since step 03 (resume after an
# interruption). Set to FALSE to refit everything.
reuse_existing_fits <- TRUE

make_biocparallel_param <- function() {
  if (!use_biocparallel || biocparallel_workers <= 1) {
    return(BiocParallel::SerialParam())
  }

  if (.Platform$OS.type == "unix") {
    return(BiocParallel::MulticoreParam(workers = biocparallel_workers))
  }

  BiocParallel::SnowParam(workers = biocparallel_workers, type = "SOCK")
}

# GDC query settings shared by all TCGA-Biolinks downloads.
tcga_data_category <- "Transcriptome Profiling"
tcga_data_type <- "Gene Expression Quantification"
tcga_workflow_type <- "STAR - Counts"
tcga_experimental_strategy <- "RNA-Seq"
tcga_count_assay <- "unstranded"
gdc_download_method <- "api"
gdc_files_per_chunk <- 20

# Clinical metadata settings. Indexed clinical data is lightweight; BCR Biotab
# supplements are richer and downloaded only when this flag is TRUE.
tcga_clinical_indexed_type <- "clinical"
download_clinical_supplement <- TRUE
tcga_clinical_data_category <- "Clinical"
tcga_clinical_data_type <- "Clinical Supplement"
tcga_clinical_data_format <- "BCR Biotab"

# Known TCGA sample types observed in GDC RNA-seq metadata.
# If GDC adds a new sample type, script 01 stops so this list can be reviewed.
normal_sample_type <- "Solid Tissue Normal"

tumor_sample_types <- "Primary Tumor"

# Non-primary tumour samples are known but excluded from tumour-vs-normal models.
uncharacterized_sample_types <- c(
  "Metastatic",
  "Recurrent Tumor",
  "Primary Blood Derived Cancer - Peripheral Blood",
  "Additional Metastatic",
  "Additional - New Primary"
)

known_sample_types <- c(
  normal_sample_type,
  tumor_sample_types,
  uncharacterized_sample_types
)

limit_projects <- function(project_ids) {
  if (is.null(debug_project_limit) || is.infinite(debug_project_limit)) {
    return(project_ids)
  }

  head(project_ids, debug_project_limit)
}

# Conservative defaults for deciding whether a project/design is usable.
min_normal_samples <- 5
min_tumor_samples <- 5
min_count <- 10
min_count_samples <- 3
max_covariate_missing_fraction <- 0.25
min_covariate_group_size <- 3
deseq_alpha <- 0.05
deseq_lfc_threshold <- log2(1.5)

# Differential expression test (step 05; changing it needs no refit, rerun from
# deg05):
# - "threshold": H0: |beta| <= deseq_lfc_threshold (DESeq2 lfcThreshold with
#   altHypothesis = "greaterAbs"). The FDR then controls genes whose change
#   exceeds the threshold.
# - "zero": H0: beta = 0, then |log2FC| > deseq_lfc_threshold as a filter.
# In both cases a DEG has padj < deseq_alpha and |raw log2FC| > deseq_lfc_threshold.
deseq_test <- "threshold"
use_lfc_shrinkage <- TRUE
lfc_shrinkage_type <- "ashr"

# Candidate covariates are tested per cancer type before entering the model.
candidate_covariates <- c(
  "sex",
  "race",
  "smoking_status"
)

# Sample subsets/models crossed with the naive and adjusted designs:
# "all" = all samples (tumor effect); "female"/"male" = one sex only;
# "sex_interaction" = both sexes with sex x condition, testing the difference
# between the male and female tumor effects.
sex_models <- c("all", "female", "male", "sex_interaction")

#### Plot settings ####

pca_blind <- FALSE
pca_max_groups <- 12

plot_dpi <- 300
pca_plot_width <- 7
pca_plot_height <- 5
ma_plot_width <- 7
ma_plot_height <- 5
volcano_plot_width <- 8
volcano_plot_height <- 8

ma_plot_ylim <- c(-5, 5)
ma_point_alpha <- 0.5
ma_point_size <- 0.8
volcano_point_alpha <- 0.55
volcano_point_size <- 0.9

deg_up_color <- "red"
deg_down_color <- "blue"
deg_neutral_color <- "grey70"

#### Metadata Audit Settings ####

metadata_column_patterns <- list(
  sex = "gender|sex",
  race = "race",
  smoking_status = "smok|tobacco|cigarette|pack",
  age = "age",
  stage_grade = "stage|grade",
  survival = "vital|death|survival|follow|recurrence|progression"
)

metadata_missing_value_labels <- c(
  "",
  "NA",
  "N/A",
  "Unknown",
  "unknown",
  "Not Reported",
  "not reported",
  "Not Available",
  "not available",
  "Not Applicable",
  "not applicable",
  "Not Evaluated",
  "not evaluated",
  "[Unknown]",
  "[Not Available]",
  "[Not Applicable]",
  "[Not Evaluated]",
  "[Not Reported]",
  "[Discrepancy]"
)

#### Paths ####

analysis_dir <- file.path("analyses", analysis_name)

data_tcga_dir <- file.path("data", analysis_name)
data_raw_tcga_dir <- file.path(data_tcga_dir, "raw")
data_processed_tcga_dir <- file.path(data_tcga_dir, "processed")
tcga_gdc_download_dir <- "GDCdata"

results_table_dir <- file.path("results", "tables", analysis_name)
results_figure_dir <- file.path("results", "figures", analysis_name)

#### Shared Output Files ####

tcga_project_status_file <- file.path(
  results_table_dir,
  "01_TCGA_project_status.csv"
)

tcga_sample_counts_file <- file.path(
  results_table_dir,
  "01_TCGA_sample_type_counts.csv"
)

tcga_metadata_availability_file <- file.path(
  results_table_dir,
  "01_TCGA_metadata_availability.csv"
)

tcga_clinical_supplement_table_dir <- file.path(
  data_processed_tcga_dir,
  "clinical_supplement_tables"
)

tcga_project_manifest_file <- file.path(
  data_processed_tcga_dir,
  "01_TCGA_project_manifest.rds"
)

tcga_design_plan_file <- file.path(
  results_table_dir,
  "03_TCGA_deseq2_design_plan.csv"
)

tcga_design_plan_rds <- file.path(
  data_processed_tcga_dir,
  "03_TCGA_deseq2_design_plan.rds"
)

tcga_covariate_report_file <- file.path(
  results_table_dir,
  "03_TCGA_covariate_report.csv"
)

tcga_deseq_fit_status_file <- file.path(
  results_table_dir,
  "04_TCGA_deseq2_fit_status.csv"
)

tcga_deseq_test_status_file <- file.path(
  results_table_dir,
  "05_TCGA_deseq2_test_status.csv"
)

tcga_design_comparison_file <- file.path(
  results_table_dir,
  "06_TCGA_design_comparison.csv"
)

tcga_sex_model_summary_file <- file.path(
  results_table_dir,
  "06_TCGA_sex_model_summary.csv"
)

#### Per-Project Output Helpers ####

tcga_project_rse_file <- function(project_id) {
  file.path(data_processed_tcga_dir, paste0(project_id, "_DGE_ready_RSE.rds"))
}

tcga_project_query_metadata_file <- function(project_id) {
  file.path(data_processed_tcga_dir, paste0(project_id, "_query_metadata.csv"))
}

tcga_project_clinical_indexed_file <- function(project_id) {
  file.path(data_processed_tcga_dir, paste0(project_id, "_clinical_indexed.csv"))
}

tcga_project_clinical_supplement_file <- function(project_id) {
  file.path(data_processed_tcga_dir, paste0(project_id, "_clinical_supplement.rds"))
}

tcga_project_metadata_file <- function(project_id) {
  file.path(data_processed_tcga_dir, paste0(project_id, "_metadata.csv"))
}

tcga_project_metadata_report_file <- function(project_id) {
  file.path(results_table_dir, paste0(project_id, "_metadata_report.xlsx"))
}

tcga_project_design_file <- function(project_id, design_type) {
  file.path(
    data_processed_tcga_dir,
    paste0(project_id, "_", design_type, "_design.rds")
  )
}

# Fitted DESeq2 object (step 04); step 05 tests from it without refitting.
tcga_project_fit_file <- function(project_id, design_type) {
  file.path(
    data_processed_tcga_dir,
    paste0(project_id, "_", design_type, "_fit_dds.rds")
  )
}

# Unfitted DESeq2 object (step 03).
tcga_project_dds_file <- function(project_id, design_type) {
  file.path(
    data_processed_tcga_dir,
    paste0(project_id, "_", design_type, "_dds.rds")
  )
}

# Result file label of a design type, e.g. "adjusted" -> "Tumor_vs_Normal_Adjusted",
# "adjusted_female" -> "Tumor_vs_Normal_Adjusted_Female",
# "adjusted_sex_interaction" -> "Sex_Interaction_Adjusted".
tcga_design_file_label <- function(design_type) {
  parts <- strsplit(design_type, "_", fixed = TRUE)[[1]]
  adjustment <- tools::toTitleCase(parts[[1]])
  sex_model <- paste(parts[-1], collapse = "_")

  if (sex_model == "sex_interaction") {
    return(paste0("Sex_Interaction_", adjustment))
  }

  paste(
    c("Tumor_vs_Normal", adjustment, if (nzchar(sex_model)) tools::toTitleCase(sex_model)),
    collapse = "_"
  )
}

tcga_deseq_results_file <- function(project_id, design_type) {
  file.path(
    results_table_dir,
    paste0(project_id, "_DGE_Results_", tcga_design_file_label(design_type), ".csv")
  )
}

tcga_deseq_raw_results_file <- function(project_id, design_type) {
  file.path(
    results_table_dir,
    paste0(project_id, "_DGE_Raw_Results_", tcga_design_file_label(design_type), ".csv")
  )
}

tcga_pca_figure_file <- function(project_id) {
  file.path(results_figure_dir, paste0(project_id, "_PCA_plot.png"))
}

tcga_pca_covariate_figure_file <- function(project_id, covariate) {
  safe_covariate <- gsub("[^A-Za-z0-9_.-]+", "_", covariate)
  file.path(results_figure_dir, paste0(project_id, "_PCA_", safe_covariate, ".png"))
}

tcga_ma_figure_file <- function(project_id) {
  file.path(results_figure_dir, paste0(project_id, "_MA_plot.png"))
}

tcga_volcano_figure_file <- function(project_id) {
  file.path(results_figure_dir, paste0(project_id, "_Volcano_plot.png"))
}

#### Create expected directories ####

dir.create(data_raw_tcga_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_processed_tcga_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tcga_clinical_supplement_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_figure_dir, recursive = TRUE, showWarnings = FALSE)

#### Libraries ####

library(AnnotationDbi)
library(ashr)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
library(ggplot2)
library(openxlsx)
library(org.Hs.eg.db)
library(pheatmap)
library(readr)
library(recount3)
library(SummarizedExperiment)
library(TCGAbiolinks)

#### Analysis Parameters ####

analysis_name <- "TCGA_deg"

# During development, limit scripts to the first N project IDs for faster runs.
# Set to Inf when running the full analysis.
debug_project_limit <- 2

# GDC query settings shared by all TCGA-Biolinks downloads.
tcga_data_category <- "Transcriptome Profiling"
tcga_data_type <- "Gene Expression Quantification"
tcga_workflow_type <- "STAR - Counts"
tcga_experimental_strategy <- "RNA-Seq"

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

tumor_sample_types <- c(
  "Primary Tumor",
  "Metastatic",
  "Recurrent Tumor",
  "Primary Blood Derived Cancer - Peripheral Blood",
  "Additional Metastatic"
)

uncharacterized_sample_types <- c(
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

# Candidate covariates are tested per cancer type before entering the model.
candidate_covariates <- c(
  "sex",
  "race",
  "ethnicity",
  "smoking_status"
)

#### Metadata Audit Settings ####

metadata_column_patterns <- list(
  sex = "gender|sex",
  race = "race",
  ethnicity = "ethnic",
  smoking_status = "smok|tobacco|cigarette|pack",
  age = "age",
  stage_grade = "stage|grade",
  survival = "vital|death|survival|follow|recurrence|progression"
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

tcga_project_design_file <- function(project_id) {
  file.path(data_processed_tcga_dir, paste0(project_id, "_design.rds"))
}

tcga_project_dds_file <- function(project_id) {
  file.path(data_processed_tcga_dir, paste0(project_id, "_dds.rds"))
}

tcga_deseq_results_file <- function(project_id) {
  file.path(results_table_dir, paste0(project_id, "_DGE_Results_Tumor_vs_Normal.csv"))
}

tcga_pca_figure_file <- function(project_id) {
  file.path(results_figure_dir, paste0(project_id, "_PCA_plot.png"))
}

tcga_ma_figure_file <- function(project_id) {
  file.path(results_figure_dir, paste0(project_id, "_MA_plot.png"))
}

tcga_volcano_figure_file <- function(project_id) {
  file.path(results_figure_dir, paste0(project_id, "_Volcano_plot.png"))
}

#### Create Expected Directories ####

dir.create(data_raw_tcga_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_processed_tcga_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tcga_clinical_supplement_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_figure_dir, recursive = TRUE, showWarnings = FALSE)

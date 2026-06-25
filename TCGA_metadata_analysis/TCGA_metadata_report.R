# ============================================================
# TCGA metadata report
# ============================================================

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(openxlsx)

# 1. Choose cancer type
cancer_type <- "TCGA-HNSC"

# 2. Query RNA-seq counts
query <- GDCquery(
  project = cancer_type,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# 3. Download data
GDCdownload(
  query,
  files.per.chunk = 5
)

# 4. Prepare SummarizedExperiment
se <- GDCprepare(query)

# 5. Extract counts and metadata
counts_mat <- assay(se, "unstranded")
metadata <- as.data.frame(colData(se))

# 6. Basic checks
dim(counts_mat)
dim(metadata)

table(metadata$sample_type, useNA = "ifany")
table(metadata$gender, useNA = "ifany")
summary(metadata$age_at_diagnosis)

# 7. Search metadata columns of interest

stage_cols <- grep(
  "stage|grade",
  colnames(metadata),
  value = TRUE,
  ignore.case = TRUE
)

smoking_cols <- grep(
  "smok|tobacco|cigarette|pack|cigar",
  colnames(metadata),
  value = TRUE,
  ignore.case = TRUE
)

alcohol_cols <- grep(
  "alcohol|drink|ethanol",
  colnames(metadata),
  value = TRUE,
  ignore.case = TRUE
)

bmi_cols <- grep(
  "bmi|body_mass_index|body.mass.index|height|weight_at_diagnosis|initial_weight",
  colnames(metadata),
  value = TRUE,
  ignore.case = TRUE
)

race_ethnicity_cols <- grep(
  "race|ethnic",
  colnames(metadata),
  value = TRUE,
  ignore.case = TRUE
)

survival_cols <- grep(
  "vital|death|survival|follow|recurrence|progression|new_tumor|first_event",
  colnames(metadata),
  value = TRUE,
  ignore.case = TRUE
)

immune_cols <- grep(
  "immune|leukocyte|lymph|inflamm|purity|stromal",
  colnames(metadata),
  value = TRUE,
  ignore.case = TRUE
)

# 8. Collect all columns

interesting_cols <- unique(c(
  stage_cols,
  smoking_cols,
  alcohol_cols,
  bmi_cols,
  race_ethnicity_cols,
  survival_cols,
  immune_cols
))

# 9. Function: create summary for each found metadata column

summarise_metadata_column <- function(metadata, col_name) {
  
  x <- metadata[[col_name]]
  n_total <- length(x)
  n_missing <- sum(is.na(x))
  n_available <- sum(!is.na(x))
  missing_pct <- round(100 * n_missing / n_total, 1)
  
  x_non_na <- x[!is.na(x)]
  
  if (length(x_non_na) == 0) {
    values_summary <- "all NA"
  } else if (is.numeric(x_non_na) || is.integer(x_non_na)) {
    values_summary <- paste(
      "min =", round(min(x_non_na), 2),
      "; median =", round(median(x_non_na), 2),
      "; mean =", round(mean(x_non_na), 2),
      "; max =", round(max(x_non_na), 2)
    )
  } else {
    tab <- sort(table(x, useNA = "ifany"), decreasing = TRUE)
    values_summary <- paste(
      paste(names(tab), as.integer(tab), sep = ": "),
      collapse = "; "
    )
  }
  
  data.frame(
    variable = col_name,
    available = n_available,
    missing = n_missing,
    missing_pct = missing_pct,
    values_summary = values_summary
  )
}

# 10. Summary for all interesting metadata columns

variable_summary <- bind_rows(
  lapply(
    interesting_cols,
    function(col) summarise_metadata_column(metadata, col)
  )
)

variable_summary

# 11. Sample type summary

sample_type_summary <- as.data.frame(
  table(metadata$sample_type, useNA = "ifany")
)

colnames(sample_type_summary) <- c("sample_type", "n")

sample_type_summary

# 12. Main one-row summary csv

metadata_summary <- data.frame(
  cancer_type = cancer_type,
  n_samples = nrow(metadata),
  n_genes = nrow(counts_mat),
  
  sample_types = paste(
    paste(sample_type_summary$sample_type, sample_type_summary$n, sep = ": "),
    collapse = "; "
  ),
  
  female = if ("gender" %in% colnames(metadata)) {
    sum(metadata$gender == "female", na.rm = TRUE)
  } else {
    NA
  },
  
  male = if ("gender" %in% colnames(metadata)) {
    sum(metadata$gender == "male", na.rm = TRUE)
  } else {
    NA
  },
  
  age_available = if ("age_at_diagnosis" %in% colnames(metadata)) {
    sum(!is.na(metadata$age_at_diagnosis))
  } else {
    NA
  },
  
  age_missing = if ("age_at_diagnosis" %in% colnames(metadata)) {
    sum(is.na(metadata$age_at_diagnosis))
  } else {
    NA
  },
  
  age_missing_pct = if ("age_at_diagnosis" %in% colnames(metadata)) {
    round(100 * mean(is.na(metadata$age_at_diagnosis)), 1)
  } else {
    NA
  },
  
  has_smoking = length(smoking_cols) > 0,
  has_alcohol = length(alcohol_cols) > 0,
  has_bmi = length(bmi_cols) > 0,
  has_race_ethnicity = length(race_ethnicity_cols) > 0,
  has_survival = length(survival_cols) > 0,
  has_stage_grade = length(stage_cols) > 0,
  has_immune_related = length(immune_cols) > 0,
  
  stage_cols = paste(stage_cols, collapse = "; "),
  smoking_cols = paste(smoking_cols, collapse = "; "),
  alcohol_cols = paste(alcohol_cols, collapse = "; "),
  bmi_cols = paste(bmi_cols, collapse = "; "),
  race_ethnicity_cols = paste(race_ethnicity_cols, collapse = "; "),
  survival_cols = paste(survival_cols, collapse = "; "),
  immune_cols = paste(immune_cols, collapse = "; ")
)

metadata_summary

# 13. version of full metadata csv

metadata_for_excel <- metadata

metadata_for_excel[] <- lapply(
  metadata_for_excel,
  function(x) {
    if (is.list(x)) {
      sapply(x, function(y) paste(y, collapse = "; "))
    } else {
      x
    }
  }
)

metadata_for_excel <- as.data.frame(metadata_for_excel)

# 14. Save csv report

wb <- createWorkbook()

addWorksheet(wb, "summary")
writeData(wb, "summary", metadata_summary)

addWorksheet(wb, "variable_summary")
writeData(wb, "variable_summary", variable_summary)

addWorksheet(wb, "sample_types")
writeData(wb, "sample_types", sample_type_summary)

addWorksheet(wb, "metadata_colnames")
writeData(wb, "metadata_colnames", data.frame(columns = colnames(metadata)))

addWorksheet(wb, "sample_metadata")
writeData(wb, "sample_metadata", metadata_for_excel, na.string = "NA")

addWorksheet(wb, "library_sizes")
writeData(wb, "library_sizes", library_size_summary)

addWorksheet(wb, "library_size_overview")
writeData(wb, "library_size_overview", library_size_overview)

saveWorkbook(
  wb,
  paste0(cancer_type, "_metadata_report.xlsx"),
  overwrite = TRUE
)

write.csv(
  variable_summary,
  paste0(cancer_type, "_metadata_summary.csv"),
  row.names = FALSE
)
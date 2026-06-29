#### Prepare Count And Metadata Inputs ####
#
# This script downloads eligible cohorts, builds SummarizedExperiment objects,
# standardizes common metadata fields, and writes inspectable metadata reports.

source(file.path("analyses", "TCGA_deg", "00_config.R"))

#### Helper Functions ####

# Reuse the same explicit sample-type classifier as in project discovery.
classify_sample_type <- function(sample_type) {
  unknown_sample_types <- setdiff(unique(as.character(sample_type)), known_sample_types)

  if (length(unknown_sample_types) > 0) {
    stop(
      "Unknown TCGA sample_type value(s): ",
      paste(sort(unknown_sample_types), collapse = "; "),
      "\nReview GDC metadata and update known_sample_types in 00_config.R."
    )
  }

  condition <- rep("uncharacterized", length(sample_type))
  condition[sample_type == normal_sample_type] <- "normal"
  condition[sample_type %in% tumor_sample_types] <- "tumor"
  factor(condition, levels = c("normal", "tumor", "uncharacterized"))
}

# Merge possible aliases for the same biological/clinical covariate.
coalesce_columns <- function(metadata, candidates) {
  existing <- intersect(candidates, colnames(metadata))

  if (length(existing) == 0) {
    return(rep(NA_character_, nrow(metadata)))
  }

  value <- as.character(metadata[[existing[[1]]]])

  if (length(existing) > 1) {
    for (column_name in existing[-1]) {
      missing <- is.na(value) | value == ""
      value[missing] <- as.character(metadata[[column_name]])[missing]
    }
  }

  value[value == ""] <- NA_character_
  value
}

# Add analysis-friendly metadata columns without deleting original metadata.
standardise_tcga_metadata <- function(metadata) {
  metadata$condition <- classify_sample_type(metadata$sample_type)
  metadata$sex <- coalesce_columns(metadata, c("gender", "sex", "sex_at_birth"))
  metadata$race <- coalesce_columns(metadata, c("race"))
  metadata$ethnicity <- coalesce_columns(metadata, c("ethnicity"))
  metadata$smoking_status <- coalesce_columns(metadata, c(
    "tobacco_smoking_history",
    "tobacco_smoking_status",
    "smoking_status",
    "cigarettes_per_day"
  ))

  if ("age_at_diagnosis" %in% colnames(metadata)) {
    metadata$age_at_diagnosis_numeric <- suppressWarnings(
      as.numeric(metadata$age_at_diagnosis)
    )
  } else {
    metadata$age_at_diagnosis_numeric <- NA_real_
  }

  metadata
}

# Flatten list columns before writing metadata to CSV/XLSX.
flatten_metadata_for_export <- function(metadata) {
  metadata[] <- lapply(metadata, function(x) {
    if (is.list(x)) {
      vapply(x, function(y) paste(y, collapse = "; "), character(1))
    } else {
      x
    }
  })

  as.data.frame(metadata)
}

# Treat blank strings as missing values in exported metadata summaries.
is_missing_metadata_value <- function(x) {
  x <- as.character(x)
  is.na(x) | x == ""
}

# Keep only row annotations that are useful downstream and easy to inspect.
slim_row_data <- function(rse) {
  keep_row_cols <- intersect(
    c("gene_id", "gene_name", "gene_type", "seqnames", "start", "end", "strand"),
    colnames(rowData(rse))
  )

  rowData(rse) <- rowData(rse)[, keep_row_cols, drop = FALSE]
  rse
}

#### Select Projects From Discovery Step ####

project_status <- read.csv(tcga_project_status_file, stringsAsFactors = FALSE)
eligible_projects <- project_status |>
  dplyr::filter(eligible_sample_counts) |>
  dplyr::pull(project_id)
eligible_projects <- limit_projects(eligible_projects)

#### Download And Prepare Each Cohort ####

for (project_id in eligible_projects) {
  message("Preparing input object for ", project_id)

  query <- GDCquery(
    project = project_id,
    data.category = tcga_data_category,
    data.type = tcga_data_type,
    workflow.type = tcga_workflow_type,
    experimental.strategy = tcga_experimental_strategy
  )

  GDCdownload(query, method = "api", files.per.chunk = 20)
  rse <- GDCprepare(query, summarizedExperiment = TRUE)

  # Use raw unstranded counts as the DESeq2 assay.
  assays(rse) <- list(counts = assay(rse, "unstranded"))
  rse <- slim_row_data(rse)

  metadata <- as.data.frame(colData(rse))
  metadata <- standardise_tcga_metadata(metadata)

  colData(rse) <- S4Vectors::DataFrame(metadata)

  keep_genes <- rowSums(assay(rse, "counts")) > 0
  rse <- rse[keep_genes, ]

  #### Save Prepared Inputs And Metadata Reports ####

  metadata_export <- flatten_metadata_for_export(as.data.frame(colData(rse)))

  saveRDS(rse, tcga_project_rse_file(project_id))
  write.csv(
    metadata_export,
    tcga_project_metadata_file(project_id),
    row.names = FALSE
  )

  sample_type_summary <- as.data.frame(
    table(metadata_export$sample_type, metadata_export$condition, useNA = "ifany"),
    stringsAsFactors = FALSE
  )
  colnames(sample_type_summary) <- c("sample_type", "condition", "n")

  metadata_availability <- data.frame(
    variable = colnames(metadata_export),
    n_available = vapply(
      metadata_export,
      function(x) sum(!is_missing_metadata_value(x)),
      integer(1)
    ),
    n_missing = vapply(
      metadata_export,
      function(x) sum(is_missing_metadata_value(x)),
      integer(1)
    ),
    stringsAsFactors = FALSE
  )
  metadata_availability$missing_fraction <-
    metadata_availability$n_missing / nrow(metadata_export)

  wb <- createWorkbook()
  addWorksheet(wb, "sample_types")
  writeData(wb, "sample_types", sample_type_summary)
  addWorksheet(wb, "metadata_availability")
  writeData(wb, "metadata_availability", metadata_availability)
  addWorksheet(wb, "sample_metadata")
  writeData(wb, "sample_metadata", metadata_export, na.string = "NA")
  saveWorkbook(wb, tcga_project_metadata_report_file(project_id), overwrite = TRUE)
}

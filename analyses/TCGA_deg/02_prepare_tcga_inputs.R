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

  values <- lapply(existing, function(column_name) {
    normalise_missing_strings(metadata[[column_name]])
  })

  Reduce(function(current, next_value) {
    missing <- is.na(current)
    current[missing] <- next_value[missing]
    current
  }, values)
}

# Harmonize common missing-value labels before merging or modelling.
normalise_missing_strings <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c(
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
    "[Not Available]",
    "[Not Applicable]"
  )] <- NA_character_
  x
}

# TCGA patient barcodes are the first 12 characters of sample barcodes.
patient_barcode_from_sample <- function(sample_barcode) {
  sample_barcode <- normalise_missing_strings(sample_barcode)
  ifelse(
    !is.na(sample_barcode) & grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", sample_barcode),
    substr(sample_barcode, 1, 12),
    NA_character_
  )
}

derive_sample_barcode <- function(metadata) {
  coalesce_columns(metadata, c(
    "sample.submitter_id",
    "sample_submitter_id",
    "barcode",
    "sample",
    "submitter_id"
  ))
}

derive_patient_barcode <- function(metadata) {
  patient_barcode <- coalesce_columns(metadata, c(
    "patient_barcode",
    "cases.submitter_id",
    "case_submitter_id",
    "bcr_patient_barcode",
    "patient",
    "submitter_id"
  ))

  missing <- is.na(patient_barcode)
  if (any(missing)) {
    patient_barcode[missing] <- patient_barcode_from_sample(
      derive_sample_barcode(metadata)[missing]
    )
  }

  patient_barcode
}

# Collapse repeated clinical rows so each patient joins at most once.
collapse_duplicate_keys <- function(metadata, key_column) {
  metadata <- metadata[!is.na(metadata[[key_column]]), , drop = FALSE]

  if (nrow(metadata) == 0) {
    return(metadata)
  }

  grouped <- split(metadata, metadata[[key_column]], drop = TRUE)
  collapsed <- lapply(grouped, function(group) {
    values <- lapply(group, function(column) {
      column <- normalise_missing_strings(column)
      available <- column[!is.na(column)]

      if (length(available) == 0) {
        return(NA_character_)
      }

      unique(available)[[1]]
    })

    as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
  })

  result <- do.call(rbind, collapsed)
  rownames(result) <- NULL
  result
}

prefix_non_key_columns <- function(metadata, prefix, key_column) {
  columns_to_prefix <- setdiff(colnames(metadata), key_column)
  colnames(metadata)[match(columns_to_prefix, colnames(metadata))] <-
    paste0(prefix, "_", columns_to_prefix)
  metadata
}

left_join_preserve_order <- function(metadata, lookup, by_x, by_y) {
  metadata$.row_order <- seq_len(nrow(metadata))
  merged <- merge(
    metadata,
    lookup,
    by.x = by_x,
    by.y = by_y,
    all.x = TRUE,
    sort = FALSE
  )
  merged <- merged[order(merged$.row_order), , drop = FALSE]
  merged$.row_order <- NULL
  rownames(merged) <- NULL
  merged
}

# Add query metadata from script 01 when GDCprepare omits useful sample fields.
merge_rnaseq_query_metadata <- function(metadata, project_id) {
  query_file <- tcga_project_query_metadata_file(project_id)

  if (!file.exists(query_file)) {
    return(metadata)
  }

  query_metadata <- read.csv(
    query_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  query_metadata$query_sample_barcode <- derive_sample_barcode(query_metadata)
  query_metadata$query_patient_barcode <- derive_patient_barcode(query_metadata)
  query_metadata <- collapse_duplicate_keys(query_metadata, "query_sample_barcode")

  if (nrow(query_metadata) == 0) {
    return(metadata)
  }

  query_metadata <- prefix_non_key_columns(
    query_metadata,
    "query",
    c("query_sample_barcode", "query_patient_barcode")
  )

  metadata$sample_barcode <- derive_sample_barcode(metadata)
  metadata <- left_join_preserve_order(
    metadata,
    query_metadata,
    by_x = "sample_barcode",
    by_y = "query_sample_barcode"
  )

  if (!"sample_type" %in% colnames(metadata)) {
    metadata$sample_type <- NA_character_
  }

  if ("query_sample_type" %in% colnames(metadata)) {
    missing_sample_type <- is.na(normalise_missing_strings(metadata$sample_type))
    metadata$sample_type[missing_sample_type] <-
      metadata$query_sample_type[missing_sample_type]
  }

  metadata$patient_barcode <- derive_patient_barcode(metadata)

  if ("query_patient_barcode" %in% colnames(metadata)) {
    missing_patient <- is.na(metadata$patient_barcode)
    metadata$patient_barcode[missing_patient] <-
      metadata$query_patient_barcode[missing_patient]
  }

  metadata
}

load_indexed_clinical_metadata <- function(project_id) {
  clinical_file <- tcga_project_clinical_indexed_file(project_id)

  if (!file.exists(clinical_file)) {
    return(NULL)
  }

  clinical <- read.csv(
    clinical_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  clinical$patient_barcode <- derive_patient_barcode(clinical)
  clinical <- collapse_duplicate_keys(clinical, "patient_barcode")

  if (nrow(clinical) == 0) {
    return(NULL)
  }

  prefix_non_key_columns(clinical, "indexed", "patient_barcode")
}

load_patient_supplement_metadata <- function(project_id) {
  supplement_file <- tcga_project_clinical_supplement_file(project_id)

  if (!file.exists(supplement_file)) {
    return(NULL)
  }

  supplement <- readRDS(supplement_file)
  patient_table_names <- grep(
    "^clinical_patient",
    names(supplement),
    ignore.case = TRUE,
    value = TRUE
  )

  if (length(patient_table_names) == 0) {
    return(NULL)
  }

  patient_table_name <- patient_table_names[[1]]

  patient_metadata <- as.data.frame(
    supplement[[patient_table_name]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  patient_metadata$patient_barcode <- derive_patient_barcode(patient_metadata)
  patient_metadata <- collapse_duplicate_keys(patient_metadata, "patient_barcode")

  if (nrow(patient_metadata) == 0) {
    return(NULL)
  }

  prefix_non_key_columns(patient_metadata, "supplement", "patient_barcode")
}

merge_clinical_metadata <- function(metadata, project_id) {
  metadata <- merge_rnaseq_query_metadata(metadata, project_id)
  metadata$sample_barcode <- derive_sample_barcode(metadata)
  metadata$patient_barcode <- derive_patient_barcode(metadata)

  clinical_indexed <- load_indexed_clinical_metadata(project_id)
  clinical_supplement <- load_patient_supplement_metadata(project_id)

  if (!is.null(clinical_indexed)) {
    metadata <- left_join_preserve_order(
      metadata,
      clinical_indexed,
      by_x = "patient_barcode",
      by_y = "patient_barcode"
    )
  }

  if (!is.null(clinical_supplement)) {
    metadata <- left_join_preserve_order(
      metadata,
      clinical_supplement,
      by_x = "patient_barcode",
      by_y = "patient_barcode"
    )
  }

  metadata
}

# Add analysis-friendly metadata columns without deleting original metadata.
standardise_tcga_metadata <- function(metadata) {
  metadata$condition <- classify_sample_type(metadata$sample_type)
  metadata$sex <- coalesce_columns(metadata, c(
    "indexed_sex_at_birth",
    "indexed_gender",
    "supplement_gender",
    "gender",
    "sex",
    "sex_at_birth"
  ))
  metadata$race <- coalesce_columns(metadata, c(
    "indexed_race",
    "supplement_race",
    "race"
  ))
  metadata$ethnicity <- coalesce_columns(metadata, c(
    "indexed_ethnicity",
    "supplement_ethnicity",
    "ethnicity"
  ))
  metadata$smoking_status <- coalesce_columns(metadata, c(
    "indexed_tobacco_smoking_status",
    "supplement_tobacco_smoking_history_indicator",
    "tobacco_smoking_history",
    "tobacco_smoking_status",
    "smoking_status",
    "cigarettes_per_day"
  ))

  metadata$age_at_diagnosis_numeric <- suppressWarnings(as.numeric(
    coalesce_columns(metadata, c(
      "indexed_age_at_diagnosis",
      "supplement_age_at_diagnosis",
      "age_at_diagnosis"
    ))
  ))

  metadata
}

summarise_clinical_merge <- function(metadata) {
  data.frame(
    metric = c(
      "n_samples",
      "n_samples_with_sample_barcode",
      "n_samples_with_patient_barcode",
      "n_samples_matched_indexed_clinical",
      "n_samples_matched_patient_supplement"
    ),
    value = c(
      nrow(metadata),
      sum(!is.na(metadata$sample_barcode)),
      sum(!is.na(metadata$patient_barcode)),
      if ("indexed_submitter_id" %in% colnames(metadata)) {
        sum(!is.na(metadata$indexed_submitter_id))
      } else {
        0
      },
      if ("supplement_bcr_patient_barcode" %in% colnames(metadata)) {
        sum(!is.na(metadata$supplement_bcr_patient_barcode))
      } else {
        0
      }
    ),
    stringsAsFactors = FALSE
  )
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

# Run GDC download/prepare calls from the raw data directory so TCGAbiolinks
# writes GDCdata/ and MANIFEST.txt under data/TCGA_deg/raw.
with_raw_data_dir <- function(expr) {
  old_wd <- setwd(data_raw_tcga_dir)
  on.exit(setwd(old_wd), add = TRUE)
  force(expr)
}

#### Select Projects From Discovery Step ####

project_status <- read.csv(tcga_project_status_file, stringsAsFactors = FALSE)
eligible_projects <- project_status |>
  dplyr::filter(eligible_sample_counts) |>
  dplyr::pull(project_id)
eligible_projects <- limit_projects(eligible_projects)

#### Download And Prepare Each Cohort ####

prepare_project_inputs <- function(project_id) {
  message("Preparing input object for ", project_id)

  query <- GDCquery(
    project = project_id,
    data.category = tcga_data_category,
    data.type = tcga_data_type,
    workflow.type = tcga_workflow_type,
    experimental.strategy = tcga_experimental_strategy
  )

  rse <- with_raw_data_dir({
    GDCdownload(
      query,
      method = "api",
      directory = tcga_gdc_download_dir,
      files.per.chunk = 20
    )
    GDCprepare(
      query,
      directory = tcga_gdc_download_dir,
      summarizedExperiment = TRUE
    )
  })

  # Use raw unstranded counts as the DESeq2 assay.
  assays(rse) <- list(counts = assay(rse, "unstranded"))
  rse <- slim_row_data(rse)

  metadata <- as.data.frame(colData(rse))
  metadata <- merge_clinical_metadata(metadata, project_id)
  metadata <- standardise_tcga_metadata(metadata)
  clinical_merge_summary <- summarise_clinical_merge(metadata)

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
  addWorksheet(wb, "clinical_merge")
  writeData(wb, "clinical_merge", clinical_merge_summary)
  addWorksheet(wb, "sample_metadata")
  writeData(wb, "sample_metadata", metadata_export, na.string = "NA")
  saveWorkbook(wb, tcga_project_metadata_report_file(project_id), overwrite = TRUE)

  invisible(project_id)
}

prepared_projects <- lapply(eligible_projects, prepare_project_inputs)

#### Discover TCGA RNA-Seq Projects ####
#
# This script queries GDC RNA-seq file metadata and clinical metadata. It does
# not download count matrices. Outputs are meant to be inspected before deciding
# which cohorts proceed.

source(file.path("analyses", "TCGA_deg", "00_config.R"))

#### Metadata Fields To Audit ####

metadata_column_patterns <- list(
  sex = "gender|sex",
  race = "race",
  ethnicity = "ethnic",
  smoking_status = "smok|tobacco|cigarette|pack",
  age = "age",
  stage_grade = "stage|grade",
  survival = "vital|death|survival|follow|recurrence|progression"
)

#### Helper Functions ####

# Assign tumor/normal/uncharacterized labels from explicitly reviewed TCGA sample types.
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
  condition
}

# Report which broad metadata categories a column name appears to represent.
matching_column_categories <- function(column_name) {
  matched <- names(metadata_column_patterns)[vapply(
    metadata_column_patterns,
    function(pattern) grepl(pattern, column_name, ignore.case = TRUE),
    logical(1)
  )]

  if (length(matched) == 0) {
    return(NA_character_)
  }

  paste(matched, collapse = ";")
}

# Summarize missingness and value diversity for one metadata column.
summarise_column_availability <- function(
  metadata,
  project_id,
  data_source,
  table_name,
  column_name
) {
  x <- as.character(metadata[[column_name]])
  n_total <- length(x)
  n_missing <- sum(is.na(x) | x == "")
  n_available <- n_total - n_missing
  n_distinct <- length(unique(x[!(is.na(x) | x == "")]))

  data.frame(
    project_id = project_id,
    data_source = data_source,
    table_name = table_name,
    column = column_name,
    matched_categories = matching_column_categories(column_name),
    n_total = n_total,
    n_available = n_available,
    n_missing = n_missing,
    missing_fraction = n_missing / n_total,
    n_distinct = n_distinct,
    stringsAsFactors = FALSE
  )
}

# Summarize all columns in one flat metadata table.
summarise_table_availability <- function(
  metadata,
  project_id,
  data_source,
  table_name
) {
  if (is.null(metadata) || ncol(metadata) == 0) {
    return(NULL)
  }

  metadata <- flatten_metadata_for_export(metadata)

  dplyr::bind_rows(lapply(
    colnames(metadata),
    function(column_name) {
      summarise_column_availability(
        metadata = metadata,
        project_id = project_id,
        data_source = data_source,
        table_name = table_name,
        column_name = column_name
      )
    }
  ))
}

# TCGAbiolinks metadata can contain list columns; CSV output needs flat values.
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

# Make table names safe for filenames.
sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9_-]+", "_", x)
}

# Clinical supplements can be either one data frame or a named list of tables.
as_named_supplement_tables <- function(supplement) {
  if (is.data.frame(supplement)) {
    return(list(clinical_supplement = supplement))
  }

  supplement_tables <- supplement
  table_names <- names(supplement_tables)

  if (is.null(table_names)) {
    names(supplement_tables) <- paste0("table_", seq_along(supplement_tables))
  } else {
    missing_names <- table_names == "" | is.na(table_names)
    names(supplement_tables)[missing_names] <- paste0(
      "table_",
      which(missing_names)
    )
  }

  supplement_tables
}

# Run GDC download/prepare calls from the raw data directory so TCGAbiolinks
# writes GDCdata/ and MANIFEST.txt under data/TCGA_deg/raw.
with_raw_data_dir <- function(expr) {
  old_wd <- setwd(data_raw_tcga_dir)
  on.exit(setwd(old_wd), add = TRUE)
  force(expr)
}

# Query RNA-seq file metadata for one project and keep errors as data.
query_rnaseq_metadata <- function(project_id) {
  tryCatch({
    query <- GDCquery(
      project = project_id,
      data.category = tcga_data_category,
      data.type = tcga_data_type,
      workflow.type = tcga_workflow_type,
      experimental.strategy = tcga_experimental_strategy
    )

    metadata <- getResults(query)
    metadata$project_id <- project_id
    metadata$condition_candidate <- classify_sample_type(metadata$sample_type)

    list(query = query, metadata = metadata, error = NA_character_)
  }, error = function(error) {
    list(query = NULL, metadata = NULL, error = conditionMessage(error))
  })
}

# Query indexed clinical metadata; this is lightweight and usually one row per case.
query_indexed_clinical <- function(project_id) {
  tryCatch({
    clinical <- GDCquery_clinic(
      project = project_id,
      type = tcga_clinical_indexed_type
    )

    list(clinical = clinical, error = NA_character_)
  }, error = function(error) {
    list(clinical = NULL, error = conditionMessage(error))
  })
}

# Query richer BCR Biotab clinical supplements when enabled in config.
query_clinical_supplement <- function(project_id) {
  if (!download_clinical_supplement) {
    return(list(supplement = NULL, error = NA_character_))
  }

  tryCatch({
    query <- GDCquery(
      project = project_id,
      data.category = tcga_clinical_data_category,
      data.type = tcga_clinical_data_type,
      data.format = tcga_clinical_data_format
    )

    supplement <- with_raw_data_dir({
      GDCdownload(query, directory = tcga_gdc_download_dir)
      GDCprepare(query, directory = tcga_gdc_download_dir)
    })

    list(supplement = supplement, error = NA_character_)
  }, error = function(error) {
    list(supplement = NULL, error = conditionMessage(error))
  })
}

# Query one project and keep each data source separate.
query_project_metadata <- function(project_id) {
  message("Querying ", project_id)

  list(
    rnaseq = query_rnaseq_metadata(project_id),
    clinical_indexed = query_indexed_clinical(project_id),
    clinical_supplement = query_clinical_supplement(project_id)
  )
}

# Write clinical supplement tables separately for inspection.
write_clinical_supplement_tables <- function(project_id, supplement) {
  if (is.null(supplement)) {
    return(invisible(NULL))
  }

  supplement_tables <- as_named_supplement_tables(supplement)

  invisible(lapply(names(supplement_tables), function(table_name) {
    table_data <- supplement_tables[[table_name]]

    if (!is.data.frame(table_data)) {
      return(NULL)
    }

    write.csv(
      flatten_metadata_for_export(table_data),
      file.path(
        tcga_clinical_supplement_table_dir,
        paste0(project_id, "_", sanitize_filename(table_name), ".csv")
      ),
      row.names = FALSE
    )
  }))
}

# Summarize clinical supplement tables for the availability report.
summarise_clinical_supplement_availability <- function(project_id, supplement) {
  if (is.null(supplement)) {
    return(NULL)
  }

  supplement_tables <- as_named_supplement_tables(supplement)

  dplyr::bind_rows(lapply(names(supplement_tables), function(table_name) {
    table_data <- supplement_tables[[table_name]]

    if (!is.data.frame(table_data)) {
      return(NULL)
    }

    summarise_table_availability(
      metadata = table_data,
      project_id = project_id,
      data_source = "clinical_supplement",
      table_name = table_name
    )
  }))
}

# Convert one query result into all outputs needed by this discovery step.
summarise_project_metadata <- function(project_id, query_result) {
  rnaseq_result <- query_result$rnaseq
  indexed_clinical_result <- query_result$clinical_indexed
  clinical_supplement_result <- query_result$clinical_supplement

  if (!is.null(rnaseq_result$metadata)) {
    metadata <- rnaseq_result$metadata

    # Save the raw query metadata so each cohort can be checked manually.
    metadata_export <- flatten_metadata_for_export(metadata)

    write.csv(
      metadata_export,
      tcga_project_query_metadata_file(project_id),
      row.names = FALSE
    )

    if (!is.null(indexed_clinical_result$clinical)) {
      write.csv(
        flatten_metadata_for_export(indexed_clinical_result$clinical),
        tcga_project_clinical_indexed_file(project_id),
        row.names = FALSE
      )
    }

    if (!is.null(clinical_supplement_result$supplement)) {
      saveRDS(
        clinical_supplement_result$supplement,
        tcga_project_clinical_supplement_file(project_id)
      )
      write_clinical_supplement_tables(
        project_id,
        clinical_supplement_result$supplement
      )
    }

    # Count all sample types before deciding what is tumor or normal.
    sample_type_count <- as.data.frame(
      table(metadata$sample_type, useNA = "ifany"),
      stringsAsFactors = FALSE
    ) |>
      dplyr::rename(sample_type = Var1, n = Freq) |>
      dplyr::mutate(project_id = project_id, .before = 1)

    # Audit RNA-seq, indexed clinical, and clinical supplement metadata columns.
    metadata_availability <- dplyr::bind_rows(
      summarise_table_availability(
        metadata = metadata,
        project_id = project_id,
        data_source = "rnaseq_query",
        table_name = "getResults"
      ),
      summarise_table_availability(
        metadata = indexed_clinical_result$clinical,
        project_id = project_id,
        data_source = "clinical_indexed",
        table_name = tcga_clinical_indexed_type
      ),
      summarise_clinical_supplement_availability(
        project_id,
        clinical_supplement_result$supplement
      )
    )

    n_normal <- sum(metadata$condition_candidate == "normal", na.rm = TRUE)
    n_tumor <- sum(metadata$condition_candidate == "tumor", na.rm = TRUE)
    n_uncharacterized <- sum(
      metadata$condition_candidate == "uncharacterized",
      na.rm = TRUE
    )

    # Eligibility here is only based on sample counts; design quality is checked later.
    project_status <- data.frame(
      project_id = project_id,
      query_ok = TRUE,
      n_samples = nrow(metadata),
      n_normal = n_normal,
      n_tumor = n_tumor,
      n_uncharacterized = n_uncharacterized,
      eligible_sample_counts =
        n_normal >= min_normal_samples && n_tumor >= min_tumor_samples,
      rnaseq_error = NA_character_,
      clinical_indexed_error = indexed_clinical_result$error,
      clinical_supplement_error = clinical_supplement_result$error,
      stringsAsFactors = FALSE
    )

    return(list(
      project_status = project_status,
      sample_type_counts = sample_type_count,
      metadata_availability = metadata_availability
    ))
  } else {
    project_status <- data.frame(
      project_id = project_id,
      query_ok = FALSE,
      n_samples = NA_integer_,
      n_normal = NA_integer_,
      n_tumor = NA_integer_,
      n_uncharacterized = NA_integer_,
      eligible_sample_counts = FALSE,
      rnaseq_error = rnaseq_result$error,
      clinical_indexed_error = indexed_clinical_result$error,
      clinical_supplement_error = clinical_supplement_result$error,
      stringsAsFactors = FALSE
    )

    return(list(
      project_status = project_status,
      sample_type_counts = NULL,
      metadata_availability = NULL
    ))
  }
}

#### Query Project-Level Metadata ####

# Recover TCGA project IDs from GDC projects.
tcga_projects <- sort(grep(
  "^TCGA",
  TCGAbiolinks:::getGDCprojects()$project_id,
  value = TRUE
))
tcga_projects <- limit_projects(tcga_projects)

project_manifest <- setNames(
  lapply(tcga_projects, query_project_metadata),
  tcga_projects
)

project_summaries <- Map(
  summarise_project_metadata,
  names(project_manifest),
  project_manifest
)

#### Save Discovery Outputs ####

saveRDS(project_manifest, tcga_project_manifest_file)

write.csv(
  dplyr::bind_rows(lapply(project_summaries, `[[`, "project_status")),
  tcga_project_status_file,
  row.names = FALSE
)

write.csv(
  dplyr::bind_rows(lapply(project_summaries, `[[`, "sample_type_counts")),
  tcga_sample_counts_file,
  row.names = FALSE
)

write.csv(
  dplyr::bind_rows(lapply(project_summaries, `[[`, "metadata_availability")),
  tcga_metadata_availability_file,
  row.names = FALSE
)

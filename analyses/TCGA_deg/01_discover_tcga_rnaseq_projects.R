#### Discover TCGA RNA-Seq Projects ####
#
# This script queries GDC metadata only. It does not download count matrices.
# Outputs are meant to be inspected before deciding which cohorts proceed.

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

# Find columns whose names suggest a metadata category of interest.
find_matching_columns <- function(metadata, pattern) {
  grep(pattern, colnames(metadata), value = TRUE, ignore.case = TRUE)
}

# Summarize missingness and value diversity for one metadata column.
summarise_column_availability <- function(metadata, project_id, column_name) {
  x <- as.character(metadata[[column_name]])
  n_total <- length(x)
  n_missing <- sum(is.na(x) | x == "")
  n_available <- n_total - n_missing
  n_distinct <- length(unique(x[!(is.na(x) | x == "")]))

  data.frame(
    project_id = project_id,
    column = column_name,
    n_total = n_total,
    n_available = n_available,
    n_missing = n_missing,
    missing_fraction = n_missing / n_total,
    n_distinct = n_distinct,
    stringsAsFactors = FALSE
  )
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

# Query one project and keep errors as data.
query_project_metadata <- function(project_id) {
  message("Querying ", project_id)

  # Keep errors in the project status table instead of stopping the scan.
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

# Convert one query result into all outputs needed by this discovery step.
summarise_project_metadata <- function(project_id, query_result) {
  if (!is.null(query_result$metadata)) {
    metadata <- query_result$metadata

    # Save the raw query metadata so each cohort can be checked manually.
    metadata_export <- flatten_metadata_for_export(metadata)

    write.csv(
      metadata_export,
      tcga_project_query_metadata_file(project_id),
      row.names = FALSE
    )

    # Count all sample types before deciding what is tumor or normal.
    sample_type_count <- as.data.frame(
      table(metadata$sample_type, useNA = "ifany"),
      stringsAsFactors = FALSE
    ) |>
      dplyr::rename(sample_type = Var1, n = Freq) |>
      dplyr::mutate(project_id = project_id, .before = 1)

    # Audit metadata columns that may later become covariates.
    matching_metadata_columns <- unique(unlist(lapply(
      metadata_column_patterns,
      function(pattern) find_matching_columns(metadata, pattern)
    )))

    metadata_availability <- dplyr::bind_rows(lapply(
      matching_metadata_columns,
      function(column_name) {
        summarise_column_availability(metadata, project_id, column_name)
      }
    ))

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
      error = NA_character_,
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
      error = query_result$error,
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

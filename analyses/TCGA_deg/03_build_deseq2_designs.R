#### Build Per-Cancer DESeq2 Designs ####
#
# This script chooses usable covariates separately for each cohort, creates
# DESeq2 objects, and writes a design plan for manual review before analysis.

source(file.path("analyses", "TCGA_deg", "00_config.R"))

#### Helper Functions ####

# Standardize common "unknown" labels before testing categorical covariates.
normalise_covariate <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% metadata_missing_value_labels] <- NA
  x <- factor(x)
  levels(x) <- make.names(levels(x), unique = TRUE)
  x
}

# A categorical covariate must be present, sufficiently complete, and variable
# in groups large enough to estimate without fragile model terms.
evaluate_categorical_covariate <- function(metadata, covariate) {
  if (!covariate %in% colnames(metadata)) {
    return(data.frame(
      covariate = covariate,
      type = "categorical",
      selected = FALSE,
      reason = "missing column",
      missing_fraction = NA_real_,
      n_groups = NA_integer_,
      smallest_group_size = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }

  x <- normalise_covariate(metadata[[covariate]])
  missing_fraction <- mean(is.na(x))
  tab <- table(x, useNA = "no")
  n_groups <- length(tab)
  smallest_group_size <- if (n_groups == 0) NA_integer_ else min(tab)

  selected <- missing_fraction <= max_covariate_missing_fraction &&
    n_groups >= 2 &&
    smallest_group_size >= min_covariate_group_size

  reason <- if (selected) {
    NA_character_
  } else if (missing_fraction > max_covariate_missing_fraction) {
    "too many missing values"
  } else if (n_groups < 2) {
    "fewer than two observed groups"
  } else {
    "smallest group below threshold"
  }

  data.frame(
    covariate = covariate,
    type = "categorical",
    selected = selected,
    reason = reason,
    missing_fraction = missing_fraction,
    n_groups = n_groups,
    smallest_group_size = smallest_group_size,
    stringsAsFactors = FALSE
  )
}

# Numeric covariates need enough observed values and variation.
evaluate_numeric_covariate <- function(metadata, covariate) {
  if (!covariate %in% colnames(metadata)) {
    return(data.frame(
      covariate = covariate,
      type = "numeric",
      selected = FALSE,
      reason = "missing column",
      missing_fraction = NA_real_,
      n_groups = NA_integer_,
      smallest_group_size = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }

  x <- suppressWarnings(as.numeric(metadata[[covariate]]))
  missing_fraction <- mean(is.na(x))
  n_distinct <- length(unique(x[!is.na(x)]))

  selected <- missing_fraction <= max_covariate_missing_fraction &&
    n_distinct >= min_covariate_group_size

  reason <- if (selected) {
    NA_character_
  } else if (missing_fraction > max_covariate_missing_fraction) {
    "too many missing values"
  } else {
    "insufficient numeric variation"
  }

  data.frame(
    covariate = covariate,
    type = "numeric",
    selected = selected,
    reason = reason,
    missing_fraction = missing_fraction,
    n_groups = n_distinct,
    smallest_group_size = NA_integer_,
    stringsAsFactors = FALSE
  )
}

# DESeq2 requires a full-rank model matrix.
formula_is_full_rank <- function(metadata, design_formula) {
  tryCatch(
    {
      design_matrix <- model.matrix(design_formula, data = metadata)
      qr(design_matrix)$rank == ncol(design_matrix)
    },
    error = function(e) FALSE
  )
}

# Use scaled age in years; convert defensively if older metadata is day-scale.
add_scaled_age_covariate <- function(metadata) {
  if (!"age_at_diagnosis_numeric" %in% colnames(metadata)) {
    return(metadata)
  }

  age <- suppressWarnings(as.numeric(metadata$age_at_diagnosis_numeric))
  if (all(is.na(age))) {
    metadata$age_at_diagnosis_years_scaled <- NA_real_
    return(metadata)
  }

  age_years <- if (stats::median(age, na.rm = TRUE) > 365) age / 365.25 else age
  metadata$age_at_diagnosis_years_scaled <- as.numeric(scale(
    age_years,
    center = TRUE,
    scale = TRUE
  ))
  metadata
}

# Select candidate covariates that pass basic availability checks.
select_covariates <- function(metadata) {
  metadata <- add_scaled_age_covariate(metadata)

  age_report <- evaluate_numeric_covariate(
    metadata,
    "age_at_diagnosis_years_scaled"
  )
  covariate_report <- dplyr::bind_rows(
    age_report,
    lapply(candidate_covariates, function(covariate) {
      evaluate_categorical_covariate(metadata, covariate)
    })
  )

  selected <- covariate_report$covariate[covariate_report$selected]

  list(
    metadata = metadata,
    selected_covariates = selected,
    covariate_report = covariate_report
  )
}

prepare_selected_covariates <- function(metadata, selected_covariates) {
  if (length(selected_covariates) == 0) {
    return(metadata)
  }

  prepared <- lapply(selected_covariates, function(covariate) {
    if (covariate == "age_at_diagnosis_years_scaled") {
      metadata[[covariate]]
    } else {
      normalise_covariate(metadata[[covariate]])
    }
  })
  names(prepared) <- selected_covariates

  metadata[selected_covariates] <- prepared
  metadata
}

# Put condition last so the tumor-vs-normal contrast is explicit.
make_design_formula <- function(covariates) {
  if (length(covariates) == 0) {
    return(as.formula("~ condition"))
  }

  as.formula(paste("~", paste(c(covariates, "condition"), collapse = " + ")))
}

format_design_formula <- function(design_formula) {
  paste(deparse(design_formula), collapse = "")
}

filter_complete_design_samples <- function(rse, metadata, selected_covariates) {
  complete_columns <- c("condition", selected_covariates)
  complete_samples <- stats::complete.cases(metadata[, complete_columns, drop = FALSE])

  filtered_metadata <- metadata[complete_samples, , drop = FALSE]
  filtered_metadata$condition <- factor(
    filtered_metadata$condition,
    levels = c("normal", "tumor")
  )

  list(
    rse = rse[, complete_samples],
    metadata = filtered_metadata
  )
}

#### Select Projects From Discovery Step ####

project_status <- read.csv(tcga_project_status_file, stringsAsFactors = FALSE)
eligible_projects <- project_status |>
  dplyr::filter(eligible_sample_counts) |>
  dplyr::pull(project_id)
eligible_projects <- limit_projects(eligible_projects)

#### Build Designs And DESeq2 Objects ####

build_project_design <- function(project_id) {
  message("Building DESeq2 design for ", project_id)

  if (!file.exists(tcga_project_rse_file(project_id))) {
    return(list(
      design_plan = data.frame(
        project_id = project_id,
        eligible = FALSE,
        reason = "prepared RSE file not found",
        n_normal = NA_integer_,
        n_tumor = NA_integer_,
        n_genes = NA_integer_,
        covariates = NA_character_,
        design_formula = NA_character_,
        stringsAsFactors = FALSE
      ),
      covariate_report = data.frame()
    ))
  }

  rse <- readRDS(tcga_project_rse_file(project_id))
  metadata <- as.data.frame(colData(rse))

  keep_samples <- metadata$condition %in% c("normal", "tumor")
  rse <- rse[, keep_samples]
  metadata <- as.data.frame(colData(rse))

  metadata$condition <- factor(metadata$condition, levels = c("normal", "tumor"))

  n_normal <- sum(metadata$condition == "normal", na.rm = TRUE)
  n_tumor <- sum(metadata$condition == "tumor", na.rm = TRUE)

  if (n_normal < min_normal_samples || n_tumor < min_tumor_samples) {
    return(list(
      design_plan = data.frame(
        project_id = project_id,
        eligible = FALSE,
        reason = "insufficient tumor/normal samples after preparation",
        n_normal = n_normal,
        n_tumor = n_tumor,
        n_genes = NA_integer_,
        covariates = NA_character_,
        design_formula = NA_character_,
        stringsAsFactors = FALSE
      ),
      covariate_report = data.frame()
    ))
  }

  selection <- select_covariates(metadata)
  metadata <- prepare_selected_covariates(
    selection$metadata,
    selection$selected_covariates
  )
  selected_covariates <- selection$selected_covariates
  covariate_report <- selection$covariate_report
  covariate_report$project_id <- project_id
  base_rse <- rse
  base_metadata <- metadata

  # Complete-case filtering can remove samples, so re-check cohort size after it.
  filtered <- filter_complete_design_samples(base_rse, base_metadata, selected_covariates)
  rse <- filtered$rse
  metadata <- filtered$metadata

  n_normal_complete <- sum(metadata$condition == "normal", na.rm = TRUE)
  n_tumor_complete <- sum(metadata$condition == "tumor", na.rm = TRUE)

  if (
    n_normal_complete < min_normal_samples ||
      n_tumor_complete < min_tumor_samples
  ) {
    return(list(
      design_plan = data.frame(
        project_id = project_id,
        eligible = FALSE,
        reason = "insufficient tumor/normal samples after covariate complete-case filtering",
        n_normal = n_normal_complete,
        n_tumor = n_tumor_complete,
        n_genes = NA_integer_,
        covariates = paste(selected_covariates, collapse = ";"),
        design_formula = NA_character_,
        stringsAsFactors = FALSE
      ),
      covariate_report = covariate_report
    ))
  }

  design_formula <- make_design_formula(selected_covariates)

  # Remove covariates from the end until the model matrix is estimable.
  while (!formula_is_full_rank(metadata, design_formula) && length(selected_covariates) > 0) {
    dropped_covariate <- selected_covariates[[length(selected_covariates)]]
    covariate_report$selected[covariate_report$covariate == dropped_covariate] <- FALSE
    covariate_report$reason[covariate_report$covariate == dropped_covariate] <-
      "dropped to restore full-rank design"
    selected_covariates <- selected_covariates[-length(selected_covariates)]
    filtered <- filter_complete_design_samples(
      base_rse,
      base_metadata,
      selected_covariates
    )
    rse <- filtered$rse
    metadata <- filtered$metadata
    design_formula <- make_design_formula(selected_covariates)
  }

  if (!formula_is_full_rank(metadata, design_formula)) {
    return(list(
      design_plan = data.frame(
        project_id = project_id,
        eligible = FALSE,
        reason = "design matrix is not full rank",
        n_normal = sum(metadata$condition == "normal"),
        n_tumor = sum(metadata$condition == "tumor"),
        n_genes = NA_integer_,
        covariates = paste(selected_covariates, collapse = ";"),
        design_formula = format_design_formula(design_formula),
        stringsAsFactors = FALSE
      ),
      covariate_report = covariate_report
    ))
  }

  colData(rse) <- S4Vectors::DataFrame(metadata)

  # Apply a simple low-count filter before saving the DESeq2 object.
  dds <- DESeqDataSet(rse, design = design_formula)
  smallest_group_size <- min(table(colData(dds)$condition))
  required_samples <- max(min_count_samples, smallest_group_size)
  keep_genes <- rowSums(counts(dds) >= min_count) >= required_samples
  dds <- dds[keep_genes, ]

  saveRDS(
    list(
      project_id = project_id,
      design_formula = design_formula,
      covariates = selected_covariates,
      dds_file = tcga_project_dds_file(project_id)
    ),
    tcga_project_design_file(project_id)
  )
  saveRDS(dds, tcga_project_dds_file(project_id))

  list(
    design_plan = data.frame(
      project_id = project_id,
      eligible = TRUE,
      reason = NA_character_,
      n_normal = sum(colData(dds)$condition == "normal"),
      n_tumor = sum(colData(dds)$condition == "tumor"),
      n_genes = nrow(dds),
      covariates = paste(selected_covariates, collapse = ";"),
      design_formula = format_design_formula(design_formula),
      stringsAsFactors = FALSE
    ),
    covariate_report = covariate_report
  )
}

design_outputs <- lapply(eligible_projects, build_project_design)

#### Save Design Plan For Review ####

design_plan_df <- dplyr::bind_rows(lapply(design_outputs, `[[`, "design_plan"))
covariate_report_df <- dplyr::bind_rows(lapply(
  design_outputs,
  `[[`,
  "covariate_report"
))

write.csv(design_plan_df, tcga_design_plan_file, row.names = FALSE)
saveRDS(design_plan_df, tcga_design_plan_rds)

write.csv(
  covariate_report_df,
  tcga_covariate_report_file,
  row.names = FALSE
)

#### Build Per-Cancer DESeq2 Designs ####
#
# This script chooses usable covariates separately for each cohort, creates
# DESeq2 objects, and writes a design plan for manual review before analysis.

source(file.path("analyses", "TCGA_deg", "00_config.R"))

#### Helper Functions ####

# Standardize common "unknown" labels before testing categorical covariates.
normalise_covariate <- function(x) {
  x <- as.character(x)
  x[x == "" | x == "not reported" | x == "Not Reported" | x == "unknown"] <- NA
  factor(x)
}

# A categorical covariate must be present, sufficiently complete, and variable.
covariate_is_usable <- function(metadata, covariate) {
  if (!covariate %in% colnames(metadata)) {
    return(FALSE)
  }

  x <- normalise_covariate(metadata[[covariate]])
  missing_fraction <- mean(is.na(x))
  tab <- table(x, useNA = "no")

  missing_fraction <= max_covariate_missing_fraction &&
    length(tab) >= 2 &&
    min(tab) >= min_covariate_group_size
}

# Numeric covariates need enough observed values and variation.
numeric_covariate_is_usable <- function(metadata, covariate) {
  if (!covariate %in% colnames(metadata)) {
    return(FALSE)
  }

  x <- suppressWarnings(as.numeric(metadata[[covariate]]))
  missing_fraction <- mean(is.na(x))

  missing_fraction <= max_covariate_missing_fraction &&
    length(unique(x[!is.na(x)])) >= min_covariate_group_size
}

# DESeq2 requires a full-rank model matrix.
formula_is_full_rank <- function(metadata, design_formula) {
  design_matrix <- model.matrix(design_formula, data = metadata)
  qr(design_matrix)$rank == ncol(design_matrix)
}

# Select candidate covariates that pass basic availability checks.
select_covariates <- function(metadata) {
  selected <- character()

  if (numeric_covariate_is_usable(metadata, "age_at_diagnosis_numeric")) {
    metadata$age_at_diagnosis_centered <- as.numeric(scale(
      as.numeric(metadata$age_at_diagnosis_numeric),
      scale = FALSE
    ))
    selected <- c(selected, "age_at_diagnosis_centered")
  }

  for (covariate in candidate_covariates) {
    if (covariate_is_usable(metadata, covariate)) {
      selected <- c(selected, covariate)
    }
  }

  selected
}

# Put condition last so the tumor-vs-normal contrast is explicit.
make_design_formula <- function(covariates) {
  if (length(covariates) == 0) {
    return(as.formula("~ condition"))
  }

  as.formula(paste("~", paste(c(covariates, "condition"), collapse = " + ")))
}

#### Select Projects From Discovery Step ####

project_status <- read.csv(tcga_project_status_file, stringsAsFactors = FALSE)
eligible_projects <- project_status |>
  dplyr::filter(eligible_sample_counts) |>
  dplyr::pull(project_id)

design_plan <- list()

#### Build Designs And DESeq2 Objects ####

for (project_id in eligible_projects) {
  message("Building DESeq2 design for ", project_id)

  rse <- readRDS(tcga_project_rse_file(project_id))
  metadata <- as.data.frame(colData(rse))

  keep_samples <- metadata$condition %in% c("normal", "tumor")
  rse <- rse[, keep_samples]
  metadata <- as.data.frame(colData(rse))

  metadata$condition <- factor(metadata$condition, levels = c("normal", "tumor"))

  n_normal <- sum(metadata$condition == "normal", na.rm = TRUE)
  n_tumor <- sum(metadata$condition == "tumor", na.rm = TRUE)

  if (n_normal < min_normal_samples || n_tumor < min_tumor_samples) {
    design_plan[[project_id]] <- data.frame(
      project_id = project_id,
      eligible = FALSE,
      reason = "insufficient tumor/normal samples after preparation",
      n_normal = n_normal,
      n_tumor = n_tumor,
      covariates = NA_character_,
      design_formula = NA_character_,
      stringsAsFactors = FALSE
    )
    next
  }

  # Start broad, then drop covariates only if complete-case/full-rank checks fail.
  selected_covariates <- select_covariates(metadata)

  for (covariate in selected_covariates) {
    if (covariate == "age_at_diagnosis_centered") {
      metadata[[covariate]] <- as.numeric(scale(
        as.numeric(metadata$age_at_diagnosis_numeric),
        scale = FALSE
      ))
    } else {
      metadata[[covariate]] <- normalise_covariate(metadata[[covariate]])
    }
  }

  complete_columns <- c("condition", selected_covariates)
  complete_samples <- stats::complete.cases(metadata[, complete_columns, drop = FALSE])
  rse <- rse[, complete_samples]
  metadata <- metadata[complete_samples, , drop = FALSE]
  metadata$condition <- factor(metadata$condition, levels = c("normal", "tumor"))

  design_formula <- make_design_formula(selected_covariates)

  # Remove covariates from the end until the model matrix is estimable.
  while (!formula_is_full_rank(metadata, design_formula) && length(selected_covariates) > 0) {
    selected_covariates <- selected_covariates[-length(selected_covariates)]
    design_formula <- make_design_formula(selected_covariates)
  }

  if (!formula_is_full_rank(metadata, design_formula)) {
    design_plan[[project_id]] <- data.frame(
      project_id = project_id,
      eligible = FALSE,
      reason = "design matrix is not full rank",
      n_normal = sum(metadata$condition == "normal"),
      n_tumor = sum(metadata$condition == "tumor"),
      covariates = paste(selected_covariates, collapse = ";"),
      design_formula = deparse(design_formula),
      stringsAsFactors = FALSE
    )
    next
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

  design_plan[[project_id]] <- data.frame(
    project_id = project_id,
    eligible = TRUE,
    reason = NA_character_,
    n_normal = sum(colData(dds)$condition == "normal"),
    n_tumor = sum(colData(dds)$condition == "tumor"),
    n_genes = nrow(dds),
    covariates = paste(selected_covariates, collapse = ";"),
    design_formula = deparse(design_formula),
    stringsAsFactors = FALSE
  )
}

#### Save Design Plan For Review ####

design_plan_df <- dplyr::bind_rows(design_plan)

write.csv(design_plan_df, tcga_design_plan_file, row.names = FALSE)
saveRDS(design_plan_df, tcga_design_plan_rds)

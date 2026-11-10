#### Build Per-Cancer DESeq2 Designs ####
#
# This script chooses usable covariates separately for each cohort, creates
# DESeq2 objects, and writes a design plan for manual review before analysis.
#
# Each cohort gets a naive and an adjusted design for every sex model in
# `sex_models` (00_config.R): all samples, females only, males only, and the
# sex x condition interaction model.

source(file.path("analyses", "TCGA_deg", "00_config.R"))

source(file.path(analysis_dir, "utils.R"))

#### Select Projects From Discovery Step ####

project_status <- read.csv(tcga_project_status_file, stringsAsFactors = FALSE)
eligible_projects <- project_status |>
  dplyr::filter(eligible_sample_counts) |>
  dplyr::pull(project_id)
eligible_projects <- limit_projects(eligible_projects)

#### Design Helpers ####

# "naive"/"adjusted" for all samples; otherwise e.g. "adjusted_female".
make_design_type <- function(adjustment, sex_model) {
  if (sex_model == "all") adjustment else paste(adjustment, sex_model, sep = "_")
}

# DESeq2 coefficient tested in script 04. With female and normal as reference
# levels, the interaction is (male tumor effect) - (female tumor effect).
design_coefficient <- function(sex_model) {
  if (sex_model == "sex_interaction") {
    "sexmale.conditiontumor"
  } else {
    "condition_tumor_vs_normal"
  }
}

is_interaction <- function(sex_model) {
  sex_model == "sex_interaction"
}

make_design_plan_row <- function(
  project_id,
  design_type,
  sex_model,
  eligible,
  reason,
  n_normal = NA_integer_,
  n_tumor = NA_integer_,
  n_genes = NA_integer_,
  covariates = character(),
  design_formula = NULL
) {
  data.frame(
    project_id = project_id,
    design_type = design_type,
    sex_model = sex_model,
    eligible = eligible,
    reason = reason,
    n_normal = n_normal,
    n_tumor = n_tumor,
    n_genes = n_genes,
    covariates = paste(covariates, collapse = ";"),
    design_formula = if (is.null(design_formula)) {
      NA_character_
    } else {
      format_design_formula(design_formula)
    },
    stringsAsFactors = FALSE
  )
}

# One ineligible row for both the naive and adjusted designs of a sex model.
make_ineligible_rows <- function(project_id, sex_model, reason, n_normal, n_tumor) {
  dplyr::bind_rows(lapply(c("naive", "adjusted"), function(adjustment) {
    make_design_plan_row(
      project_id, make_design_type(adjustment, sex_model), sex_model,
      FALSE, reason, n_normal, n_tumor
    )
  }))
}

# Sample groups used for size checks and the low-count gene filter:
# condition, or sex x condition for the interaction model.
design_group_counts <- function(metadata, sex_model) {
  sex_group <- if (is_interaction(sex_model)) metadata$sex else rep("all", nrow(metadata))
  table(
    factor(metadata$condition, levels = c("normal", "tumor")),
    sex_group
  )
}

has_enough_samples <- function(metadata, sex_model) {
  counts <- design_group_counts(metadata, sex_model)
  ncol(counts) > 0 &&
    all(counts["normal", ] >= min_normal_samples) &&
    all(counts["tumor", ] >= min_tumor_samples)
}

# Columns that must be complete in the adjusted model, besides condition.
required_design_columns <- function(sex_model, covariates) {
  c(if (is_interaction(sex_model)) "sex", covariates)
}

# Standardise sex labels to a female/male factor; other labels are an error.
standardise_sex <- function(sex) {
  x <- tolower(trimws(as.character(sex)))
  x[x %in% tolower(metadata_missing_value_labels)] <- NA
  unexpected <- setdiff(unique(x[!is.na(x)]), c("female", "male"))
  if (length(unexpected) > 0) {
    stop("Unexpected sex labels: ", paste(unexpected, collapse = ", "))
  }
  factor(x, levels = c("female", "male"))
}

# Keep the samples of one sex model. "all" is returned unchanged.
select_sex_model_samples <- function(rse, metadata, sex_model) {
  if (sex_model == "all") {
    return(list(rse = rse, metadata = metadata))
  }

  metadata$sex <- standardise_sex(metadata$sex)
  keep <- if (is_interaction(sex_model)) {
    !is.na(metadata$sex)
  } else {
    !is.na(metadata$sex) & metadata$sex == sex_model
  }

  rse <- rse[, keep]
  rse$sex <- metadata$sex[keep]
  list(rse = rse, metadata = metadata[keep, , drop = FALSE])
}

build_and_save_dds <- function(
  rse,
  project_id,
  design_type,
  sex_model,
  design_formula,
  covariates
) {
  dds <- DESeqDataSet(rse, design = design_formula)
  smallest_group_size <- min(design_group_counts(colData(dds), sex_model))
  required_samples <- max(min_count_samples, smallest_group_size)
  keep_genes <- rowSums(counts(dds) >= min_count) >= required_samples
  dds <- dds[keep_genes, ]

  dds_file <- tcga_project_dds_file(project_id, design_type)
  design_file <- tcga_project_design_file(project_id, design_type)
  saveRDS(dds, dds_file)
  saveRDS(
    list(
      project_id = project_id,
      design_type = design_type,
      sex_model = sex_model,
      design_formula = design_formula,
      coefficient = design_coefficient(sex_model),
      covariates = covariates,
      n_normal = sum(colData(dds)$condition == "normal"),
      n_tumor = sum(colData(dds)$condition == "tumor"),
      dds_file = dds_file
    ),
    design_file
  )

  dds
}

#### Build Designs And DESeq2 Objects ####

# Build the naive and adjusted designs of one sex model.
build_sex_model_designs <- function(project_id, rse, metadata, sex_model) {
  naive_type <- make_design_type("naive", sex_model)
  adjusted_type <- make_design_type("adjusted", sex_model)
  interaction <- is_interaction(sex_model)

  selected <- select_sex_model_samples(rse, metadata, sex_model)
  rse <- selected$rse
  metadata <- selected$metadata

  n_normal <- sum(metadata$condition == "normal", na.rm = TRUE)
  n_tumor <- sum(metadata$condition == "tumor", na.rm = TRUE)

  if (!has_enough_samples(metadata, sex_model)) {
    return(list(
      design_plan = make_ineligible_rows(
        project_id, sex_model,
        "insufficient tumor/normal samples after preparation",
        n_normal, n_tumor
      ),
      covariate_report = data.frame()
    ))
  }

  # The naive model uses every available sample of this sex model and no
  # covariates.
  rse$condition <- metadata$condition
  naive_formula <- make_design_formula(character(), interaction)
  naive_dds <- build_and_save_dds(
    rse, project_id, naive_type, sex_model, naive_formula, character()
  )
  naive_plan <- make_design_plan_row(
    project_id = project_id,
    design_type = naive_type,
    sex_model = sex_model,
    eligible = TRUE,
    reason = NA_character_,
    n_normal = sum(colData(naive_dds)$condition == "normal"),
    n_tumor = sum(colData(naive_dds)$condition == "tumor"),
    n_genes = nrow(naive_dds),
    covariates = character(),
    design_formula = naive_formula
  )

  # Sex is a covariate only in the all-samples model; sex-specific models
  # either fix sex or model it through the interaction term.
  candidates <- if (sex_model == "all") {
    candidate_covariates
  } else {
    setdiff(candidate_covariates, "sex")
  }
  selection <- select_covariates(metadata, candidates)
  metadata <- prepare_selected_covariates(
    selection$metadata,
    selection$selected_covariates
  )
  selected_covariates <- selection$selected_covariates
  covariate_report <- selection$covariate_report
  covariate_report$project_id <- project_id
  covariate_report$design_type <- adjusted_type

  adjusted_ineligible <- function(reason, n_normal, n_tumor, ...) {
    list(
      design_plan = dplyr::bind_rows(
        naive_plan,
        make_design_plan_row(
          project_id, adjusted_type, sex_model, FALSE, reason,
          n_normal, n_tumor, ...
        )
      ),
      covariate_report = covariate_report
    )
  }

  # If no covariates are usable, the adjusted model would duplicate the naive
  # model exactly, so only the canonical naive result is produced.
  if (length(selected_covariates) == 0) {
    return(adjusted_ineligible(
      "no usable covariates; adjusted design equals naive design",
      n_normal, n_tumor,
      design_formula = naive_formula
    ))
  }

  base_rse <- rse
  base_metadata <- metadata

  # Complete-case filtering can remove samples, so re-check cohort size after it.
  filtered <- filter_complete_design_samples(
    base_rse, base_metadata,
    required_design_columns(sex_model, selected_covariates)
  )
  rse <- filtered$rse
  metadata <- filtered$metadata

  if (!has_enough_samples(metadata, sex_model)) {
    return(adjusted_ineligible(
      "insufficient tumor/normal samples after covariate complete-case filtering",
      sum(metadata$condition == "normal", na.rm = TRUE),
      sum(metadata$condition == "tumor", na.rm = TRUE),
      covariates = selected_covariates
    ))
  }

  design_formula <- make_design_formula(selected_covariates, interaction)

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
      required_design_columns(sex_model, selected_covariates)
    )
    rse <- filtered$rse
    metadata <- filtered$metadata
    design_formula <- make_design_formula(selected_covariates, interaction)
  }

  if (!formula_is_full_rank(metadata, design_formula)) {
    return(adjusted_ineligible(
      "design matrix is not full rank",
      sum(metadata$condition == "normal"),
      sum(metadata$condition == "tumor"),
      covariates = selected_covariates,
      design_formula = design_formula
    ))
  }

  # Full-rank restoration can remove every selected covariate. In that case,
  # the adjusted model is again identical to the already-saved naive model.
  if (length(selected_covariates) == 0) {
    return(adjusted_ineligible(
      "all covariates dropped; adjusted design equals naive design",
      n_normal, n_tumor,
      design_formula = naive_formula
    ))
  }

  colData(rse) <- S4Vectors::DataFrame(metadata)
  dds <- build_and_save_dds(
    rse, project_id, adjusted_type, sex_model, design_formula, selected_covariates
  )

  list(
    design_plan = dplyr::bind_rows(
      naive_plan,
      make_design_plan_row(
        project_id = project_id,
        design_type = adjusted_type,
        sex_model = sex_model,
        eligible = TRUE,
        reason = NA_character_,
        n_normal = sum(colData(dds)$condition == "normal"),
        n_tumor = sum(colData(dds)$condition == "tumor"),
        n_genes = nrow(dds),
        covariates = selected_covariates,
        design_formula = design_formula
      )
    ),
    covariate_report = covariate_report
  )
}

build_project_design <- function(project_id) {
  message("Building DESeq2 designs for ", project_id)

  if (!file.exists(tcga_project_rse_file(project_id))) {
    return(list(
      design_plan = make_design_plan_row(
        project_id, "naive", "all", FALSE, "prepared RSE file not found"
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

  outputs <- lapply(sex_models, function(sex_model) {
    if (sex_model != "all" && !"sex" %in% colnames(metadata)) {
      return(list(
        design_plan = make_ineligible_rows(
          project_id, sex_model, "missing sex column", NA_integer_, NA_integer_
        ),
        covariate_report = data.frame()
      ))
    }
    build_sex_model_designs(project_id, rse, metadata, sex_model)
  })

  list(
    design_plan = dplyr::bind_rows(lapply(outputs, `[[`, "design_plan")),
    covariate_report = dplyr::bind_rows(lapply(outputs, `[[`, "covariate_report"))
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

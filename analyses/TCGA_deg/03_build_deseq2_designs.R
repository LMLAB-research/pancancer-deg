#### Build Per-Cancer DESeq2 Designs ####
#
# This script chooses usable covariates separately for each cohort, creates
# DESeq2 objects, and writes a design plan for manual review before analysis.

source(file.path("analyses", "TCGA_deg", "00_config.R"))
source(file.path(analysis_dir, "utils.R"))

#### Select Projects From Discovery Step ####

project_status <- read.csv(tcga_project_status_file, stringsAsFactors = FALSE)
eligible_projects <- project_status |>
  dplyr::filter(eligible_sample_counts) |>
  dplyr::pull(project_id)
eligible_projects <- limit_projects(eligible_projects)

#### Build Designs And DESeq2 Objects ####

make_design_plan_row <- function(
  project_id,
  design_type,
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

build_and_save_dds <- function(
  rse,
  project_id,
  design_type,
  design_formula,
  covariates
) {
  dds <- DESeqDataSet(rse, design = design_formula)
  smallest_group_size <- min(table(colData(dds)$condition))
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
      design_formula = design_formula,
      covariates = covariates,
      n_normal = sum(colData(dds)$condition == "normal"),
      n_tumor = sum(colData(dds)$condition == "tumor"),
      dds_file = dds_file
    ),
    design_file
  )

  dds
}

build_project_design <- function(project_id) {
  message("Building DESeq2 design for ", project_id)

  if (!file.exists(tcga_project_rse_file(project_id))) {
    return(list(
      design_plan = make_design_plan_row(
        project_id, "naive", FALSE, "prepared RSE file not found"
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
      design_plan = make_design_plan_row(
        project_id, "naive", FALSE,
        "insufficient tumor/normal samples after preparation",
        n_normal, n_tumor
      ),
      covariate_report = data.frame()
    ))
  }

  # The naive model uses every available tumor/normal sample and no covariates.
  naive_formula <- make_design_formula(character())
  naive_dds <- build_and_save_dds(
    rse, project_id, "naive", naive_formula, character()
  )
  naive_plan <- make_design_plan_row(
    project_id = project_id,
    design_type = "naive",
    eligible = TRUE,
    reason = NA_character_,
    n_normal = sum(colData(naive_dds)$condition == "normal"),
    n_tumor = sum(colData(naive_dds)$condition == "tumor"),
    n_genes = nrow(naive_dds),
    covariates = character(),
    design_formula = naive_formula
  )

  selection <- select_covariates(metadata)
  metadata <- prepare_selected_covariates(
    selection$metadata,
    selection$selected_covariates
  )
  selected_covariates <- selection$selected_covariates
  covariate_report <- selection$covariate_report
  covariate_report$project_id <- project_id

  # If no covariates are usable, the adjusted model would duplicate the naive
  # model exactly, so only the canonical naive result is produced.
  if (length(selected_covariates) == 0) {
    return(list(
      design_plan = dplyr::bind_rows(
        naive_plan,
        make_design_plan_row(
          project_id, "adjusted", FALSE,
          "no usable covariates; adjusted design equals naive design",
          n_normal, n_tumor,
          design_formula = naive_formula
        )
      ),
      covariate_report = covariate_report
    ))
  }

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
      design_plan = dplyr::bind_rows(
        naive_plan,
        make_design_plan_row(
          project_id, "adjusted", FALSE,
          "insufficient tumor/normal samples after covariate complete-case filtering",
          n_normal_complete, n_tumor_complete,
          covariates = selected_covariates
        )
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
      design_plan = dplyr::bind_rows(
        naive_plan,
        make_design_plan_row(
          project_id, "adjusted", FALSE, "design matrix is not full rank",
          sum(metadata$condition == "normal"),
          sum(metadata$condition == "tumor"),
          covariates = selected_covariates,
          design_formula = design_formula
        )
      ),
      covariate_report = covariate_report
    ))
  }

  # Full-rank restoration can remove every selected covariate. In that case,
  # the adjusted model is again identical to the already-saved naive model.
  if (length(selected_covariates) == 0) {
    return(list(
      design_plan = dplyr::bind_rows(
        naive_plan,
        make_design_plan_row(
          project_id, "adjusted", FALSE,
          "all covariates dropped; adjusted design equals naive design",
          n_normal, n_tumor,
          design_formula = naive_formula
        )
      ),
      covariate_report = covariate_report
    ))
  }

  colData(rse) <- S4Vectors::DataFrame(metadata)
  dds <- build_and_save_dds(
    rse, project_id, "adjusted", design_formula, selected_covariates
  )

  list(
    design_plan = dplyr::bind_rows(
      naive_plan,
      make_design_plan_row(
        project_id = project_id,
        design_type = "adjusted",
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

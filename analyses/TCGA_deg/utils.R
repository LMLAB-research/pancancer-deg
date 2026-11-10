# Covariate and design helpers -----------------------------------------------

# Convert known missing-value labels to NA and standardize factor levels.
normalise_covariate <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% metadata_missing_value_labels] <- NA
  x <- factor(x)
  levels(x) <- make.names(levels(x), unique = TRUE)
  x
}

# Select categorical covariates with acceptable missingness, at least two
# observed levels, and sufficient samples in the smallest level.
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

# Select numeric covariates using missingness and distinct-value thresholds.
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

# Check whether the model matrix can be constructed and has full rank.
formula_is_full_rank <- function(metadata, design_formula) {
  tryCatch(
    {
      design_matrix <- model.matrix(design_formula, data = metadata)
      qr(design_matrix)$rank == ncol(design_matrix)
    },
    error = function(e) FALSE
  )
}

# Treat age as days when its median exceeds 365, then standardize in years.
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

# Evaluate age and candidate categorical covariates for model inclusion.
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

# Append condition or the sex-by-condition terms after adjustment covariates.
make_design_formula <- function(covariates, interaction = FALSE) {
  condition_term <- if (interaction) "sex * condition" else "condition"
  as.formula(paste("~", paste(c(covariates, condition_term), collapse = " + ")))
}

format_design_formula <- function(design_formula) {
  gsub("\\s+", " ", paste(deparse(design_formula), collapse = ""))
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


# Plotting and result helpers ------------------------------------------------

is_plot_group_usable <- function(metadata, variable) {
  if (!variable %in% colnames(metadata)) {
    return(FALSE)
  }
  
  x <- metadata[[variable]]
  if (is.numeric(x)) {
    return(FALSE)
  }
  
  n_groups <- length(unique(x[!is.na(x)]))
  n_groups >= 2 && n_groups <= pca_max_groups
}

classify_deg_direction <- function(results_df, lfc = results_df$log2FoldChange) {
  direction <- rep("not_significant", nrow(results_df))
  significant <- !is.na(results_df$padj) &
    results_df$padj < deseq_alpha &
    !is.na(lfc) &
    abs(lfc) > deseq_lfc_threshold
  
  direction[significant & lfc > 0] <- "up"
  direction[significant & lfc < 0] <- "down"
  
  factor(direction, levels = c("down", "not_significant", "up"))
}

deg_direction_colors <- c(
  down = deg_down_color,
  not_significant = deg_neutral_color,
  up = deg_up_color
)

# Join rowData annotations and recover missing symbols from org.Hs.eg.db.
annotate_results <- function(results_df, row_data) {
  results_df$gene_id <- rownames(results_df)
  
  if (!"gene_id" %in% colnames(row_data)) {
    row_data$gene_id <- rownames(row_data)
  }
  
  annotation_cols <- intersect(
    c("gene_id", "gene_name", "gene_type", "seqnames", "start", "end", "strand"),
    colnames(row_data)
  )
  row_annotation <- row_data[
    match(results_df$gene_id, row_data$gene_id),
    annotation_cols,
    drop = FALSE
  ]
  
  row_annotation <- row_annotation[
    ,
    setdiff(colnames(row_annotation), "gene_id"),
    drop = FALSE
  ]
  results_df <- cbind(results_df, row_annotation)
  
  if (!"gene_name" %in% colnames(results_df)) {
    results_df$gene_name <- NA_character_
  }
  
  missing_gene_name <- is.na(results_df$gene_name) | results_df$gene_name == ""
  if (any(missing_gene_name)) {
    ensembl_id <- gsub("\\..*", "", results_df$gene_id[missing_gene_name])
    results_df$gene_name[missing_gene_name] <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = ensembl_id,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  }
  
  results_df
}

plot_project_pca <- function(vst_data, project_id, variable, output_file) {
  pca_plot <- plotPCA(vst_data, intgroup = variable) +
    theme_minimal() +
    labs(
      title = paste(project_id, "VST PCA"),
      color = variable
    )
  
  ggsave(
    output_file,
    plot = pca_plot,
    width = pca_plot_width,
    height = pca_plot_height,
    dpi = plot_dpi
  )
}

plot_project_pcas <- function(dds, project_id, design_info) {
  metadata <- as.data.frame(colData(dds))
  pca_variables <- unique(c("condition", design_info$covariates))
  pca_variables <- pca_variables[vapply(
    pca_variables,
    function(variable) is_plot_group_usable(metadata, variable),
    logical(1)
  )]
  
  vst_data <- vst(dds, blind = pca_blind)
  
  lapply(pca_variables, function(variable) {
    output_file <- if (variable == "condition") {
      tcga_pca_figure_file(project_id)
    } else {
      tcga_pca_covariate_figure_file(project_id, variable)
    }
    
    plot_project_pca(vst_data, project_id, variable, output_file)
  })
  
  invisible(pca_variables)
}

plot_project_ma <- function(results_df, project_id) {
  plot_df <- results_df[
    !is.na(results_df$baseMean) &
      results_df$baseMean > 0 &
      !is.na(results_df$log2FoldChange),
    ,
    drop = FALSE
  ]
  
  if (nrow(plot_df) == 0) {
    return(FALSE)
  }
  
  plot_df$direction <- classify_deg_direction(plot_df)
  
  ma_plot <- ggplot(
    plot_df,
    aes(x = baseMean, y = log2FoldChange, color = direction)
  ) +
    geom_point(alpha = ma_point_alpha, size = ma_point_size) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    scale_x_log10() +
    scale_color_manual(values = deg_direction_colors, drop = FALSE) +
    coord_cartesian(ylim = ma_plot_ylim) +
    theme_minimal() +
    labs(
      title = paste(project_id, "MA plot"),
      x = "Mean normalized count",
      y = "Shrunken log2 fold change",
      color = "Direction"
    )
  
  ggsave(
    tcga_ma_figure_file(project_id),
    plot = ma_plot,
    width = ma_plot_width,
    height = ma_plot_height,
    dpi = plot_dpi
  )
  
  TRUE
}

plot_project_volcano <- function(results_df, project_id) {
  plot_df <- results_df[
    !is.na(results_df$padj) &
      !is.na(results_df$log2FoldChange) &
      !is.na(results_df$gene_name) &
      results_df$gene_name != "",
    ,
    drop = FALSE
  ]
  
  if (nrow(plot_df) == 0) {
    return(FALSE)
  }
  
  plot_df$direction <- classify_deg_direction(plot_df)
  plot_df$neg_log10_padj <- -log10(pmax(plot_df$padj, .Machine$double.xmin))
  
  volcano_plot <- ggplot(
    plot_df,
    aes(x = log2FoldChange, y = neg_log10_padj, color = direction)
  ) +
    geom_point(alpha = volcano_point_alpha, size = volcano_point_size) +
    geom_vline(
      xintercept = c(-deseq_lfc_threshold, deseq_lfc_threshold),
      linetype = "dashed",
      color = "grey40"
    ) +
    geom_hline(
      yintercept = -log10(deseq_alpha),
      linetype = "dashed",
      color = "grey40"
    ) +
    scale_color_manual(values = deg_direction_colors, drop = FALSE) +
    theme_minimal() +
    labs(
      title = paste(project_id, "tumor vs normal"),
      x = "Shrunken log2 fold change",
      y = "-log10 adjusted p-value",
      color = "Direction"
    )
  
  ggsave(
    tcga_volcano_figure_file(project_id),
    plot = volcano_plot,
    width = volcano_plot_width,
    height = volcano_plot_height,
    dpi = plot_dpi
  )
  
  TRUE
}

summarise_results <- function(project_id, design_info, results_df, raw_results_df) {
  data.frame(
    project_id = project_id,
    design_type = design_info$design_type,
    status = "completed",
    reason = NA_character_,
    n_normal = design_info$n_normal,
    n_tumor = design_info$n_tumor,
    n_results = nrow(results_df),
    n_raw_results = nrow(raw_results_df),
    deseq_alpha = deseq_alpha,
    deseq_lfc_threshold = deseq_lfc_threshold,
    use_lfc_shrinkage = use_lfc_shrinkage,
    lfc_shrinkage_type = if (use_lfc_shrinkage) lfc_shrinkage_type else NA_character_,
    n_padj_below_alpha = sum(results_df$padj < deseq_alpha, na.rm = TRUE),
    n_abs_lfc_above_threshold = sum(
      abs(results_df$log2FoldChange) > deseq_lfc_threshold,
      na.rm = TRUE
    ),
    n_padj_below_alpha_and_abs_lfc_above_threshold = sum(
      results_df$padj < deseq_alpha &
        abs(results_df$log2FoldChange) > deseq_lfc_threshold,
      na.rm = TRUE
    ),
    biocparallel_backend = class(design_info$BPPARAM)[[1]],
    biocparallel_workers = BiocParallel::bpnworkers(design_info$BPPARAM),
    covariates = paste(design_info$covariates, collapse = ";"),
    design_formula = format_design_formula(design_info$design_formula),
    stringsAsFactors = FALSE
  )
}

skip_status <- function(project_id, design_type, reason) {
  data.frame(
    project_id = project_id,
    design_type = design_type,
    status = "skipped",
    reason = reason,
    n_normal = NA_integer_,
    n_tumor = NA_integer_,
    n_results = NA_integer_,
    n_raw_results = NA_integer_,
    deseq_alpha = deseq_alpha,
    deseq_lfc_threshold = deseq_lfc_threshold,
    use_lfc_shrinkage = use_lfc_shrinkage,
    lfc_shrinkage_type = if (use_lfc_shrinkage) lfc_shrinkage_type else NA_character_,
    n_padj_below_alpha = NA_integer_,
    n_abs_lfc_above_threshold = NA_integer_,
    n_padj_below_alpha_and_abs_lfc_above_threshold = NA_integer_,
    biocparallel_backend = NA_character_,
    biocparallel_workers = NA_integer_,
    covariates = NA_character_,
    design_formula = NA_character_,
    stringsAsFactors = FALSE
  )
}

failure_status <- function(project_id, design_type, error) {
  data.frame(
    project_id = project_id,
    design_type = design_type,
    status = "failed",
    reason = conditionMessage(error),
    n_normal = NA_integer_,
    n_tumor = NA_integer_,
    n_results = NA_integer_,
    n_raw_results = NA_integer_,
    deseq_alpha = deseq_alpha,
    deseq_lfc_threshold = deseq_lfc_threshold,
    use_lfc_shrinkage = use_lfc_shrinkage,
    lfc_shrinkage_type = if (use_lfc_shrinkage) lfc_shrinkage_type else NA_character_,
    n_padj_below_alpha = NA_integer_,
    n_abs_lfc_above_threshold = NA_integer_,
    n_padj_below_alpha_and_abs_lfc_above_threshold = NA_integer_,
    biocparallel_backend = NA_character_,
    biocparallel_workers = NA_integer_,
    covariates = NA_character_,
    design_formula = NA_character_,
    stringsAsFactors = FALSE
  )
}

run_project_deseq2 <- function(project_id, design_type, create_plots) {
  message("Running ", design_type, " DESeq2 for ", project_id)
  
  dds_file <- tcga_project_dds_file(project_id, design_type)
  design_file <- tcga_project_design_file(project_id, design_type)
  
  if (!file.exists(dds_file)) {
    return(skip_status(project_id, design_type, "DDS file not found"))
  }
  
  if (!file.exists(design_file)) {
    return(skip_status(project_id, design_type, "design file not found"))
  }
  
  dds <- readRDS(dds_file)
  design_info <- readRDS(design_file)
  design_info$BPPARAM <- make_biocparallel_param()
  parallel_enabled <- use_biocparallel &&
    BiocParallel::bpnworkers(design_info$BPPARAM) > 1
  
  pca_variables <- if (create_plots) {
    plot_project_pcas(dds, project_id, design_info)
  } else {
    character()
  }
  
  dds <- DESeq(
    dds,
    parallel = parallel_enabled,
    BPPARAM = design_info$BPPARAM
  )
  
  # Use normal as the reference level for the tumor-versus-normal contrast.
  contrast <- c("condition", "tumor", "normal")
  raw_results <- results(dds, contrast = contrast, alpha = deseq_alpha)
  shrunken_results <- if (use_lfc_shrinkage) {
    lfcShrink(
      dds,
      contrast = contrast,
      res = raw_results,
      type = lfc_shrinkage_type,
      parallel = parallel_enabled,
      BPPARAM = design_info$BPPARAM
    )
  } else {
    raw_results
  }
  
  row_data <- as.data.frame(rowData(dds))
  raw_results_df <- annotate_results(as.data.frame(raw_results), row_data)
  shrunken_results_df <- annotate_results(as.data.frame(shrunken_results), row_data)
  
  raw_results_df <- raw_results_df |>
    dplyr::filter(!is.na(padj)) |>
    dplyr::arrange(padj)
  shrunken_results_df <- shrunken_results_df |>
    dplyr::filter(!is.na(padj)) |>
    dplyr::arrange(padj)
  
  # Save results and quality-control outputs -----------------------------------
  
  write.csv(
    raw_results_df,
    tcga_deseq_raw_results_file(project_id, design_type),
    row.names = FALSE
  )
  write.csv(
    shrunken_results_df,
    tcga_deseq_results_file(project_id, design_type),
    row.names = FALSE
  )
  
  ma_created <- create_plots && plot_project_ma(shrunken_results_df, project_id)
  volcano_created <- create_plots && plot_project_volcano(shrunken_results_df, project_id)
  
  metadata(dds)$tcga_deg <- list(
    pca_variables = pca_variables,
    ma_created = ma_created,
    volcano_created = volcano_created,
    use_lfc_shrinkage = use_lfc_shrinkage,
    lfc_shrinkage_type = if (use_lfc_shrinkage) lfc_shrinkage_type else NA_character_,
    biocparallel_backend = class(design_info$BPPARAM)[[1]],
    biocparallel_workers = BiocParallel::bpnworkers(design_info$BPPARAM),
    design_type = design_type,
    raw_results_file = tcga_deseq_raw_results_file(project_id, design_type),
    shrunken_results_file = tcga_deseq_results_file(project_id, design_type)
  )
  saveRDS(dds, dds_file)
  
  summarise_results(project_id, design_info, shrunken_results_df, raw_results_df)
}

run_project_deseq2_safe <- function(project_id, design_type, create_plots) {
  tryCatch(
    run_project_deseq2(project_id, design_type, create_plots),
    error = function(e) failure_status(project_id, design_type, e)
  )
}

# Sex-specific analysis helpers (step 05) ------------------------------------
# Reuse the shared design, annotation, and DEG utilities defined above.

# Run control ----------------------------------------------------------------

# Distinguish an ineligible model from a model that fails unexpectedly.
skip_model <- function(reason) {
  stop(structure(list(message = reason, call = NULL),
                 class = c("model_skip", "error", "condition")))
}

cell_counts <- function(meta) {
  counts <- table(
    factor(meta$sex, c("female", "male")),
    factor(meta$condition, c("normal", "tumor"))
  )
  c(
    female_normal = counts["female", "normal"],
    male_normal = counts["male", "normal"],
    female_tumor = counts["female", "tumor"],
    male_tumor = counts["male", "tumor"]
  )
}

validate_labels <- function(x, allowed, column) {
  # Normalize case and whitespace; reject unrecognized labels explicitly.
  x <- tolower(trimws(as.character(x)))
  unexpected <- setdiff(unique(x[!is.na(x)]), allowed)
  if (length(unexpected)) {
    stop(column, " has unexpected categories: ", paste(unexpected, collapse = ", "))
  }
  factor(x, levels = allowed)
}

# Project data ---------------------------------------------------------------

load_sex_project <- function(project_id, adjusted_projects) {
  rse <- readRDS(tcga_project_rse_file(project_id))
  meta <- as.data.frame(colData(rse))
  
  # Restore sample identifiers from the existing TCGA metadata.
  if (!"barcode" %in% names(meta)) {
    stop("Missing barcode column")
  }
  
  sample_ids <- trimws(as.character(meta$barcode))
  
  if (anyNA(sample_ids) ||
      any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids)) {
    stop("Sample barcodes must be present and unique")
  }
  
  colnames(rse) <- sample_ids
  rownames(meta) <- sample_ids
  if (!all(c("sex", "condition") %in% names(meta))) {
    stop("Missing sex/condition column")
  }
  meta$sex <- validate_labels(meta$sex, c("female", "male"), "sex")
  # Exclude samples without a normal or tumor classification.
  condition_values <- tolower(trimws(as.character(meta$condition)))
  
  unexpected <- setdiff(
    unique(condition_values[!is.na(condition_values)]),
    c("normal", "tumor", "uncharacterized")
  )
  
  if (length(unexpected)) {
    stop("Unexpected condition labels: ", paste(unexpected, collapse = ", "))
  }
  
  keep <- !is.na(condition_values) &
    condition_values %in% c("normal", "tumor")
  
  message(project_id, ": excluding ", sum(!keep), " unclassified samples")
  
  rse <- rse[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  meta$condition <- validate_labels(meta$condition, c("normal", "tumor"), "condition")
  if (is.null(rownames(rse)) || anyNA(rownames(rse)) || anyDuplicated(rownames(rse))) {
    stop("Gene IDs must be present and unique")
  }
  if (is.null(colnames(rse)) || anyNA(colnames(rse)) || anyDuplicated(colnames(rse))) {
    stop("Sample IDs must be present and unique")
  }
  covariates <- character()
  design_source <- "No eligible adjusted design; condition-only covariates"
  if (project_id %in% adjusted_projects) {
    design_path <- tcga_project_design_file(project_id, "adjusted")
    info <- readRDS(design_path)
    if (is.null(info$covariates)) {
      stop("Adjusted design has no covariates field")
    }
    covariates <- setdiff(info$covariates, c("sex", "condition"))
    design_source <- "Covariates inherited from eligible adjusted design"
  }
  meta <- add_scaled_age_covariate(meta)
  # Replace a constant age vector's undefined scaled values with zero so
  # the covariate is dropped as constant without excluding valid samples.
  age_name <- "age_at_diagnosis_years_scaled"
  if (age_name %in% covariates && "age_at_diagnosis_numeric" %in% names(meta)) {
    age <- as.numeric(as.character(meta$age_at_diagnosis_numeric))
    if (length(unique(age[is.finite(age)])) == 1L) {
      meta[[age_name]] <- ifelse(is.finite(age), 0, NA_real_)
    }
  }
  list(rse = rse, meta = meta, covariates = covariates,
       design_source = design_source)
}

# Sample preparation ---------------------------------------------------------

sample_audit <- function(meta, stage, required, total, previous,
                         n_before, n_after = NA_integer_) {
  n <- cell_counts(meta)
  data.frame(
    stage = stage,
    n_samples = nrow(meta),
    n_normal = sum(meta$condition == "normal", na.rm = TRUE),
    n_tumor = sum(meta$condition == "tumor", na.rm = TRUE),
    n_female_normal = n[["female_normal"]], n_male_normal = n[["male_normal"]],
    n_female_tumor = n[["female_tumor"]], n_male_tumor = n[["male_tumor"]],
    n_missing_sex = sum(is.na(meta$sex)),
    n_missing_condition = sum(is.na(meta$condition)),
    n_missing_required = sum(!complete.cases(meta[, required, drop = FALSE])),
    n_excluded_stage = previous - nrow(meta),
    n_excluded_total = total - nrow(meta),
    n_genes_before = n_before,
    n_genes_after = n_after
  )
}

prepare_model_samples <- function(model, project) {
  meta <- project$meta
  missing_columns <- setdiff(model$used, names(meta))
  if (length(missing_columns)) {
    stop("Missing covariates: ", paste(missing_columns, collapse = ", "))
  }
  
  meta <- prepare_selected_covariates(meta, model$used)
  for (variable in model$used) {
    if (is.numeric(meta[[variable]])) {
      meta[[variable]][!is.finite(meta[[variable]])] <- NA_real_
    }
  }
  
  required <- c("sex", "condition", model$used)
  total <- nrow(meta)
  n_before <- nrow(project$rse)
  input_audit <- sample_audit(meta, "input", required, total, total, n_before)
  
  selected <- if (model$analysis == "interaction") {
    !is.na(meta$sex)
  } else {
    meta$sex %in% model$analysis
  }
  meta <- meta[selected, , drop = FALSE]
  sex_audit <- sample_audit(meta, "sex_selection", required, total, total, n_before)
  
  # Apply complete-case filtering independently to each analysis, so missing
  # covariates in one sex do not exclude samples from the other sex's fit.
  previous <- nrow(meta)
  complete <- complete.cases(meta[, required, drop = FALSE])
  meta <- droplevels(meta[complete, , drop = FALSE])
  complete_audit <- sample_audit(
    meta, "complete_cases", required, total, previous, n_before
  )
  
  model$meta <- meta
  model$total <- total
  model$n_before <- n_before
  model$audit <- bind_rows(input_audit, sex_audit, complete_audit)
  model
}

# Model design and fitting ---------------------------------------------------

prepare_model_design <- function(model) {
  meta <- model$meta
  n <- cell_counts(meta)
  eligible <- c(
    female = n[["female_normal"]] >= min_normal_samples &
      n[["female_tumor"]] >= min_tumor_samples,
    male = n[["male_normal"]] >= min_normal_samples &
      n[["male_tumor"]] >= min_tumor_samples
  )
  enough <- if (model$analysis == "interaction") {
    all(eligible)
  } else {
    eligible[[model$analysis]]
  }
  if (!enough) {
    model$skip_reason <- "Insufficient normal/tumor samples after complete-case filtering"
    return(model)
  }
  
  # Drop covariates that become constant in this subset and record the reason.
  # Retain variable covariates for the subsequent rank check.
  constant <- vapply(meta[model$used], function(x) length(unique(x)) < 2L, logical(1))
  if (any(constant)) {
    model$dropped <- paste0(model$used[constant], ": constant after sample filtering")
  }
  model$used <- model$used[!constant]
  meta$condition <- factor(meta$condition, c("normal", "tumor"))
  if (model$analysis == "interaction") {
    meta$sex <- factor(meta$sex, c("female", "male"))
  }
  model$meta <- meta
  model$formula <- make_design_formula(
    model$used, interaction = model$analysis == "interaction"
  )
  
  # Skip rank-deficient designs rather than silently removing potentially
  # confounded covariates from the adjustment set.
  matrix <- model.matrix(model$formula, meta)
  decomposition <- qr(matrix)
  if (decomposition$rank < ncol(matrix)) {
    dependent <- seq.int(decomposition$rank + 1L, ncol(matrix))
    aliased <- colnames(matrix)[decomposition$pivot[dependent]]
    model$skip_reason <- paste(
      "Rank-deficient design; confounding/collinearity unresolved:",
      paste(aliased, collapse = ", ")
    )
  } else if (nrow(matrix) <= ncol(matrix)) {
    model$skip_reason <- "No residual degrees of freedom"
  }
  model
}

filter_model_genes <- function(model, rse) {
  meta <- model$meta
  rse <- rse[, match(rownames(meta), colnames(rse))]
  colData(rse) <- S4Vectors::DataFrame(meta)
  group_sizes <- if (model$analysis == "interaction") {
    cell_counts(meta)
  } else {
    as.integer(table(meta$condition))
  }
  
  # Set the minimum expressing-sample count to the larger of the configured
  # minimum and the smallest group size, capped by total samples.
  model$required_samples <- min(
    ncol(rse), max(min_count_samples, min(group_sizes))
  )
  model$keep_genes <- rowSums(assay(rse, "counts") >= min_count) >= model$required_samples
  model$n_after <- sum(model$keep_genes)
  model$rse <- rse[model$keep_genes, ]
  model$audit <- bind_rows(model$audit, sample_audit(
    meta, "gene_filter", c("sex", "condition", model$used),
    model$total, nrow(meta), model$n_before, model$n_after
  ))
  model
}

fit_sex_model <- function(model) {
  # Keep the requested unpaired design; warn if repeated donors are present.
  # Accounting for within-donor dependence requires a separate design review.
  if ("patient_barcode" %in% names(model$meta) &&
      anyDuplicated(na.omit(model$meta$patient_barcode))) {
    warning("Repeated patient IDs: requested model does not account for within-patient dependence")
  }
  
  dds <- DESeqDataSet(model$rse, design = model$formula)
  bp <- make_biocparallel_param()
  parallel_enabled <- use_biocparallel && BiocParallel::bpnworkers(bp) > 1
  dds <- DESeq(dds, parallel = parallel_enabled, BPPARAM = bp)
  list(dds = dds, bp = bp, parallel_enabled = parallel_enabled)
}

# Result tables --------------------------------------------------------------

result_effect <- function(tab) {
  if (use_lfc_shrinkage) {
    tab$log2FoldChange_shrunken
  } else {
    tab$log2FoldChange
  }
}

extract_model_results <- function(fit, analysis) {
  dds <- fit$dds
  coefficient <- if (analysis == "interaction") {
    "sexmale.conditiontumor"
  } else {
    "condition_tumor_vs_normal"
  }
  if (!coefficient %in% resultsNames(dds)) {
    stop("Expected coefficient absent: ", coefficient,
         "; available: ", paste(resultsNames(dds), collapse = ", "))
  }
  
  # With female and normal as references, the interaction estimates the male
  # minus female tumor effect; its Wald test is separate from both stratified fits.
  raw <- if (analysis == "interaction") {
    results(dds, name = coefficient, alpha = deseq_alpha)
  } else {
    results(dds, contrast = c("condition", "tumor", "normal"), alpha = deseq_alpha)
  }
  tab <- annotate_results(as.data.frame(raw), as.data.frame(rowData(dds)))
  tab$log2FoldChange_shrunken <- NA_real_
  tab$posterior_sd <- NA_real_
  if (use_lfc_shrinkage) {
    shrunk <- lfcShrink(
      dds, res = raw, type = lfc_shrinkage_type,
      parallel = fit$parallel_enabled, BPPARAM = fit$bp
    )
    tab$log2FoldChange_shrunken <- shrunk$log2FoldChange
    tab$posterior_sd <- shrunk$lfcSE
  }
  
  # Preserve raw DESeq2 estimates and tests; store the shrunken effect and
  # the shrinkage method's returned uncertainty in separate columns.
  tab$shrinkage_applied <- use_lfc_shrinkage
  tab$tested <- !is.na(tab$pvalue)
  tab$padj_available <- !is.na(tab$padj)
  tab$test_status <- "pvalue_unavailable"
  tab$test_status[tab$tested] <- "independent_filter"
  tab$test_status[tab$padj_available] <- "adjusted_pvalue_available"
  
  # Combine per-model BH-adjusted p-values with a separate effect-size cutoff.
  # This cutoff does not formally test whether the effect exceeds the threshold.
  effect <- result_effect(tab)
  tab$significant_deg <- !is.na(tab$padj) & tab$padj < deseq_alpha &
    !is.na(effect) & abs(effect) > deseq_lfc_threshold
  tab
}

make_gene_universe <- function(rse, keep, tab) {
  # Track all input genes to distinguish filtered or untested genes from
  # tested, nonsignificant genes; do not perform enrichment here.
  universe <- annotate_results(
    data.frame(row.names = rownames(rse)), as.data.frame(rowData(rse))
  )
  idx <- match(universe$gene_id, tab$gene_id)
  universe$passed_count_filter <- keep
  universe$tested <- !is.na(idx) & tab$tested[idx] %in% TRUE
  universe$padj_available <- !is.na(idx) & tab$padj_available[idx] %in% TRUE
  universe$test_status <- ifelse(!keep, "low_count_filter", tab$test_status[idx])
  universe
}

# Figures --------------------------------------------------------------------

plot_sex_results <- function(tab, project_id, analysis, paths) {
  # Use the selected effect estimate and unmodified adjusted p-values for
  # classification; do not exclude genes with missing annotations.
  tab$effect <- result_effect(tab)
  tab$direction <- classify_deg_direction(tab, lfc = tab$effect)
  effect_label <- if (use_lfc_shrinkage) {
    "Shrunken log2 fold change"
  } else {
    "Log2 fold change"
  }
  comparison <- if (analysis == "interaction") {
    "Male effect minus female effect"
  } else {
    "Tumor vs normal"
  }
  title <- paste(project_id, analysis, comparison, sep = " | ")
  ma <- ggplot(tab[is.finite(tab$effect) & tab$baseMean > 0, ],
               aes(baseMean, effect, color = direction)) +
    geom_point(alpha = ma_point_alpha, size = ma_point_size) +
    geom_hline(yintercept = 0, linetype = "dashed") + scale_x_log10() +
    scale_color_manual(values = deg_direction_colors, drop = FALSE) +
    coord_cartesian(ylim = ma_plot_ylim) + theme_minimal() +
    labs(title = title, x = "Mean normalized count", y = effect_label, color = "Direction")
  volcano <- tab[is.finite(tab$effect) & !is.na(tab$padj), ]
  volcano$neg_log10_padj <- -log10(pmax(volcano$padj, .Machine$double.xmin))
  vp <- ggplot(volcano, aes(effect, neg_log10_padj, color = direction)) +
    geom_point(alpha = volcano_point_alpha, size = volcano_point_size) +
    geom_vline(xintercept = c(-deseq_lfc_threshold, deseq_lfc_threshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(deseq_alpha), linetype = "dashed") +
    scale_color_manual(values = deg_direction_colors, drop = FALSE) + theme_minimal() +
    labs(title = title, x = effect_label, y = "-log10 adjusted p-value", color = "Direction")
  # Save both figures even when the volcano plot has no eligible points.
  ggsave(paths$ma, ma, width = ma_plot_width, height = ma_plot_height, dpi = plot_dpi)
  ggsave(
    paths$volcano, vp,
    width = volcano_plot_width, height = volcano_plot_height, dpi = plot_dpi
  )
}

# Sample-level QC plots ------------------------------------------------------

# Select a reproducible, group-balanced subset for large distance heatmaps.
# PCA always includes all samples; only the distance heatmap is capped.
select_heatmap_samples <- function(sample_metadata, max_samples = 150L) {
  sample_ids <- rownames(sample_metadata)
  if (length(sample_ids) <= max_samples) {
    return(sample_ids)
  }
  
  groups <- interaction(
    sample_metadata$sex,
    sample_metadata$condition,
    drop = TRUE
  )
  grouped_ids <- split(sample_ids, groups, drop = TRUE)
  group_quota <- max(1L, floor(max_samples / length(grouped_ids)))
  
  selected <- unlist(lapply(grouped_ids, function(ids) {
    n <- min(length(ids), group_quota)
    ids[unique(as.integer(round(seq(1, length(ids), length.out = n))))]
  }), use.names = FALSE)
  
  if (length(selected) < max_samples) {
    remaining <- sample_ids[!sample_ids %in% selected]
    selected <- c(selected, head(remaining, max_samples - length(selected)))
  }
  
  sample_ids[sample_ids %in% selected]
}

# Plot VST-based PCA and sample distances for the fitted model's samples.
# These visualizations are descriptive QC, not inputs to DEG testing.
plot_sex_model_qc <- function(dds, project_id, analysis, paths) {
  vsd <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
  sample_metadata <- as.data.frame(SummarizedExperiment::colData(vsd))
  sample_metadata$condition <- factor(
    sample_metadata$condition, levels = c("normal", "tumor")
  )
  
  comparison <- paste(project_id, analysis, sep = " | ")
  intgroup <- if (analysis == "interaction") {
    c("sex", "condition")
  } else {
    "condition"
  }
  
  # Two variable genes are needed to draw a two-dimensional PCA.
  transformed <- SummarizedExperiment::assay(vsd)
  gene_variances <- matrixStats::rowVars(transformed)
  n_variable <- sum(is.finite(gene_variances) & gene_variances > 0)
  pca_created <- FALSE
  
  if (n_variable >= 2L && ncol(vsd) >= 3L) {
    pca_data <- DESeq2::plotPCA(
      vsd,
      intgroup = intgroup,
      ntop = min(500L, n_variable),
      returnData = TRUE
    )
    percent_variance <- round(100 * attr(pca_data, "percentVar"), 1)
    
    if (analysis == "interaction") {
      pca_plot <- ggplot2::ggplot(
        pca_data,
        ggplot2::aes(x = PC1, y = PC2, color = condition, shape = sex)
      )
    } else {
      pca_plot <- ggplot2::ggplot(
        pca_data,
        ggplot2::aes(x = PC1, y = PC2, color = condition)
      )
    }
    
    pca_plot <- pca_plot +
      ggplot2::geom_point(size = 2.5, alpha = 0.8) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste(comparison, "VST PCA"),
        x = paste0("PC1 (", percent_variance[[1]], "% variance)"),
        y = paste0("PC2 (", percent_variance[[2]], "% variance)")
      )
    
    ggplot2::ggsave(
      paths$pca,
      plot = pca_plot,
      width = pca_plot_width,
      height = pca_plot_height,
      dpi = plot_dpi
    )
    pca_created <- TRUE
  } else {
    warning("PCA unavailable: fewer than two variable genes or three samples")
  }
  
  # Limit large distance heatmaps to 150 samples while retaining each group.
  heatmap_ids <- select_heatmap_samples(sample_metadata)
  sample_distances <- stats::dist(t(transformed[, heatmap_ids, drop = FALSE]))
  distance_matrix <- as.matrix(sample_distances)
  annotation_variables <- if (analysis == "interaction") {
    c("sex", "condition")
  } else {
    "condition"
  }
  annotations <- sample_metadata[
    heatmap_ids, annotation_variables, drop = FALSE
  ]
  annotations[] <- lapply(annotations, factor)
  
  title <- sprintf(
    "%s | VST sample distances (%d/%d samples)",
    comparison,
    length(heatmap_ids),
    ncol(vsd)
  )
  
  pheatmap::pheatmap(
    distance_matrix,
    clustering_distance_rows = sample_distances,
    clustering_distance_cols = sample_distances,
    annotation_row = annotations,
    annotation_col = annotations,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = title,
    filename = paths$sample_distance,
    width = 9,
    height = 9
  )
  
  invisible(list(
    pca_created = pca_created,
    n_pca_samples = ncol(vsd),
    n_heatmap_samples = length(heatmap_ids),
    heatmap_sample_ids = heatmap_ids
  ))
}

# Reports and output ---------------------------------------------------------

model_report <- function(model, status, reason = NA_character_, warnings = character()) {
  formula <- if (is.null(model$formula)) NA_character_ else format_design_formula(model$formula)
  status <- data.frame(
    project_id = model$project_id,
    analysis_type = model$analysis,
    status = status,
    reason = reason,
    design_formula = formula,
    design_source = model$design_source,
    used_covariates = paste(model$used, collapse = ";"),
    dropped_covariates = paste(model$dropped, collapse = ";"),
    n_tested_genes = model$n_tested,
    n_significant_deg = model$n_significant,
    warnings = paste(unique(warnings), collapse = "; "),
    output_directory = model$output_directory
  )
  audit <- model$audit
  if (is.null(audit)) {
    audit <- data.frame(
      stage = "unavailable",
      n_samples = NA_integer_,
      n_normal = NA_integer_, n_tumor = NA_integer_,
      n_female_normal = NA_integer_, n_male_normal = NA_integer_,
      n_female_tumor = NA_integer_, n_male_tumor = NA_integer_,
      n_missing_sex = NA_integer_, n_missing_condition = NA_integer_,
      n_missing_required = NA_integer_,
      n_excluded_stage = NA_integer_, n_excluded_total = NA_integer_,
      n_genes_before = NA_integer_, n_genes_after = NA_integer_
    )
  }
  audit <- cbind(
    project_id = model$project_id, analysis_type = model$analysis, audit
  )
  list(status = status, audit = audit)
}

save_model_outputs <- function(model, fit, tab, universe, context, warnings) {
  prefix <- paste(model$project_id, model$analysis, sep = "_")
  
  comparison_label <- if (model$analysis == "interaction") {
    "Sex_Interaction"
  } else {
    paste0("Tumor_vs_Normal_", tools::toTitleCase(model$analysis))
  }
  
  paths <- list(
    results = file.path(
      context$table_dir,
      paste0(prefix, "_results.csv")
    ),
    raw_results = file.path(
      context$table_dir,
      paste0(model$project_id, "_DGE_Raw_Results_", comparison_label, ".csv")
    ),
    shrunken_results = file.path(
      context$table_dir,
      paste0(model$project_id, "_DGE_Results_", comparison_label, ".csv")
    ),
    universe = file.path(
      context$table_dir,
      paste0(prefix, "_gene_universe.csv")
    ),
    dds = file.path(
      context$data_dir,
      paste0(prefix, "_dds.rds")
    ),
    ma = file.path(
      context$figure_dir,
      paste0(prefix, "_MA.png")
    ),
    volcano = file.path(
      context$figure_dir,
      paste0(prefix, "_volcano.png")
    )
  )
  qc_paths <- list(
    pca = file.path(context$figure_dir, paste0(prefix, "_PCA.png")),
    sample_distance = file.path(
      context$figure_dir, paste0(prefix, "_sample_distance.png")
    )
  )

  # Remove previous QC figures so a failed QC plot cannot leave a stale file.
  unlink(unlist(qc_paths))

  # Export one complete table with raw statistics and shrunken effect sizes.
  # Keep genes with missing adjusted p-values and rank available padj first.
  # Preserve the complete sex-specific results table.
  ranked <- order(
    tab$padj,
    -abs(result_effect(tab)),
    tab$gene_id,
    na.last = TRUE
  )
  
  write.csv(
    tab[ranked, , drop = FALSE],
    paths$results,
    row.names = FALSE
  )
  
  # Match the original DESeq2 CSV export format.
  required_cols <- c(
    "baseMean", "log2FoldChange", "lfcSE",
    "stat", "pvalue", "padj", "gene_id"
  )
  
  if (!all(required_cols %in% names(tab))) {
    stop("Missing required DESeq2 result columns")
  }
  
  # Same filtering and sorting as the original analysis.
  filtered <- tab[!is.na(tab$padj), , drop = FALSE]
  filtered <- filtered[order(filtered$padj), , drop = FALSE]
  
  annotation_cols <- intersect(
    c(
      "gene_id", "gene_name", "gene_type",
      "seqnames", "start", "end", "strand"
    ),
    names(filtered)
  )
  
  # Raw results: original DESeq2 estimates and statistics.
  raw_cols <- c(
    "baseMean", "log2FoldChange", "lfcSE",
    "stat", "pvalue", "padj"
  )
  
  raw_table <- filtered[
    , c(raw_cols, annotation_cols), drop = FALSE
  ]
  
  write.csv(
    raw_table,
    paths$raw_results,
    row.names = FALSE
  )
  
  # Shrunken results: same structure as the original ashr output.
  shrunken_cols <- if (use_lfc_shrinkage) {
    setdiff(raw_cols, "stat")
  } else {
    raw_cols
  }
  
  shrunken_table <- filtered[
    , c(shrunken_cols, annotation_cols), drop = FALSE
  ]
  
  if (use_lfc_shrinkage) {
    if (!all(c("log2FoldChange_shrunken", "posterior_sd") %in% names(filtered))) {
      stop("Missing shrinkage result columns")
    }
    
    shrunken_table$log2FoldChange <- filtered$log2FoldChange_shrunken
    shrunken_table$lfcSE <- filtered$posterior_sd
  }
  
  write.csv(
    shrunken_table,
    paths$shrunken_results,
    row.names = FALSE
  )
  
  # Preserve the complete gene universe.
  write.csv(
    universe,
    paths$universe,
    row.names = FALSE
  )
  plot_sex_results(tab, model$project_id, model$analysis, paths)
  
  dds <- fit$dds
  # A QC failure must not invalidate completed DESeq2 inference.
  qc_info <- tryCatch(
    plot_sex_model_qc(dds, model$project_id, model$analysis, qc_paths),
    error = function(e) {
      warning("Sample-level QC plot generation failed: ", conditionMessage(e))
      list(qc_failed = TRUE)
    }
  )
  
  # Record only the QC files that were actually generated.
  qc_generated <- qc_paths[file.exists(unlist(qc_paths))]
  paths <- c(paths, qc_generated)
  
  metadata(dds)$sex_stratified <- list(
    project_id = model$project_id,
    analysis_type = model$analysis,
    formula = format_design_formula(model$formula),
    used_covariates = model$used,
    dropped_covariates = model$dropped,
    coefficient = if (model$analysis == "interaction") {
      "sexmale.conditiontumor"
    } else {
      "condition_tumor_vs_normal"
    },
    required_count_samples = model$required_samples,
    sample_ids = colnames(dds),
    paths = paths,
    qc = qc_info,
    warnings = unique(warnings),
    settings = context$settings,
    BPPARAM = class(fit$bp)[1],
    workers = BiocParallel::bpnworkers(fit$bp),
    session = sessionInfo()
  )
  saveRDS(dds, paths$dds)
  paths
}

# Individual analysis execution ----------------------------------------------

run_sex_model <- function(project_id, analysis, project, context) {
  model <- list(
    project_id = project_id, analysis = analysis,
    used = character(), dropped = character(),
    design_source = NA_character_,
    n_tested = NA_integer_, n_significant = NA_integer_,
    output_directory = context$table_dir
  )
  warnings <- character()
  
  # Handle skips and errors per analysis, preserving available audit records
  # while allowing the remaining analyses and projects to continue.
  tryCatch(withCallingHandlers({
    if (inherits(project, "error")) {
      stop(project)
    }
    model$used <- project$covariates
    model$design_source <- project$design_source

    model <- prepare_model_samples(model, project)
    model <- prepare_model_design(model)
    if (!is.null(model$skip_reason)) {
      skip_model(model$skip_reason)
    }
    model <- filter_model_genes(model, project$rse)
    if (!model$n_after) {
      skip_model("No genes pass the low-count filter")
    }
    
    fit <- fit_sex_model(model)
    tab <- extract_model_results(fit, analysis)
    universe <- make_gene_universe(project$rse, model$keep_genes, tab)
    model$n_tested <- sum(tab$tested)
    model$n_significant <- sum(tab$significant_deg)
    save_model_outputs(model, fit, tab, universe, context, warnings)
    model_report(model, "success", warnings = warnings)
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }),
  model_skip = function(e) model_report(model, "skipped", conditionMessage(e), warnings),
  error = function(e) model_report(model, "failed", conditionMessage(e), warnings))
}

create_run_context <- function() {
  # Use one fixed output directory per output type; reruns overwrite files.
  packages <- c(
    "DESeq2", "ashr", "SummarizedExperiment", "org.Hs.eg.db",
    "ggplot2", "pheatmap", "matrixStats"
  )
  settings <- list(
    min_normal_samples = min_normal_samples,
    min_tumor_samples = min_tumor_samples,
    min_count = min_count,
    min_count_samples = min_count_samples,
    alpha = deseq_alpha,
    lfc_threshold = deseq_lfc_threshold,
    shrinkage = use_lfc_shrinkage,
    shrinkage_type = lfc_shrinkage_type,
    r_version = R.version.string,
    packages = sapply(packages, function(p) as.character(packageVersion(p)))
  )
  context <- list(
    table_dir = file.path(results_table_dir, "05_sex_specific"),
    data_dir = file.path(data_processed_tcga_dir, "05_sex_specific"),
    figure_dir = file.path(results_figure_dir, "05_sex_specific"),
    settings = settings
  )
  for (path in context[c("table_dir", "data_dir", "figure_dir")]) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  context
}

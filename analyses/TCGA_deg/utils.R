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
select_covariates <- function(metadata, candidates = candidate_covariates) {
  metadata <- add_scaled_age_covariate(metadata)
  
  age_report <- evaluate_numeric_covariate(
    metadata,
    "age_at_diagnosis_years_scaled"
  )
  covariate_report <- dplyr::bind_rows(
    age_report,
    lapply(candidates, function(covariate) {
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

# DEG call as a factor for plots: the `deg` column written by step 05 (based
# on the raw log2FC), else computed from the table's own log2FC.
deg_factor <- function(results_df) {
  if ("deg" %in% colnames(results_df)) {
    factor(results_df$deg, levels = c("down", "not_significant", "up"))
  } else {
    classify_deg_direction(results_df)
  }
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
  
  plot_df$direction <- deg_factor(plot_df)
  
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
  
  plot_df$direction <- deg_factor(plot_df)
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

# Status rows -----------------------------------------------------------------

# One status row; columns not given are filled with NA when rows are combined.
status_row <- function(project_id, design_type, status, reason = NA_character_, ...) {
  data.frame(
    project_id = project_id,
    design_type = design_type,
    status = status,
    reason = reason,
    ...,
    stringsAsFactors = FALSE
  )
}

run_safely <- function(fun, project_id, design_type, ...) {
  tryCatch(
    fun(project_id, design_type, ...),
    error = function(e) status_row(project_id, design_type, "failed", conditionMessage(e))
  )
}

# Step 04: fit -------------------------------------------------------------------

# A fitted object can be reused if it is newer than both step-03 files of its
# design (so rerunning step 03 always forces a refit).
fit_is_current <- function(fit_file, dds_file, design_file) {
  file.exists(fit_file) &&
    file.mtime(fit_file) > file.mtime(dds_file) &&
    file.mtime(fit_file) > file.mtime(design_file)
}

# Fit the DESeq2 model of one design (size factors, dispersions, coefficients)
# and save the fitted object; no test is run here. With reuse_existing_fits,
# designs already fitted since step 03 are skipped, so an interrupted step 04
# resumes where it stopped.
fit_project_deseq2 <- function(project_id, design_type, create_plots) {
  dds_file <- tcga_project_dds_file(project_id, design_type)
  design_file <- tcga_project_design_file(project_id, design_type)
  fit_file <- tcga_project_fit_file(project_id, design_type)
  if (!file.exists(dds_file) || !file.exists(design_file)) {
    return(status_row(project_id, design_type, "skipped", "design or DDS file not found"))
  }

  design_info <- readRDS(design_file)
  if (reuse_existing_fits && fit_is_current(fit_file, dds_file, design_file)) {
    message("Reusing fitted ", design_type, " DESeq2 for ", project_id)
    return(status_row(
      project_id, design_type, "reused",
      "fitted object newer than its step-03 design",
      n_normal = design_info$n_normal,
      n_tumor = design_info$n_tumor,
      design_formula = format_design_formula(design_info$design_formula)
    ))
  }

  message("Fitting ", design_type, " DESeq2 for ", project_id)
  dds <- readRDS(dds_file)
  bpparam <- make_biocparallel_param()
  parallel_enabled <- use_biocparallel && BiocParallel::bpnworkers(bpparam) > 1

  pca_variables <- if (create_plots) {
    plot_project_pcas(dds, project_id, design_info)
  } else {
    character()
  }

  dds <- DESeq(dds, parallel = parallel_enabled, BPPARAM = bpparam)
  metadata(dds)$tcga_deg_fit <- list(
    design_type = design_type,
    coefficient = design_info$coefficient,
    pca_variables = pca_variables
  )
  # Write to a temporary file, then rename, so that an interrupted save never
  # leaves a truncated file that would look like a finished fit.
  tmp_file <- paste0(fit_file, ".tmp")
  saveRDS(dds, tmp_file)
  if (!file.rename(tmp_file, fit_file)) {
    stop("Could not move ", tmp_file, " to ", fit_file)
  }

  status_row(
    project_id, design_type, "completed",
    n_normal = design_info$n_normal,
    n_tumor = design_info$n_tumor,
    n_genes = nrow(dds),
    biocparallel_workers = BiocParallel::bpnworkers(bpparam),
    design_formula = format_design_formula(design_info$design_formula)
  )
}

# Step 05: test --------------------------------------------------------------

# DESeq2 results for the coefficient chosen in step 03 (the tumor effect, or
# for the interaction the male minus female tumor effect), with the test set
# by deseq_test in 00_config.R.
test_coefficient <- function(dds, coefficient) {
  if (!deseq_test %in% c("threshold", "zero")) {
    stop("deseq_test must be \"threshold\" or \"zero\"")
  }
  if (deseq_test == "threshold") {
    results(
      dds, name = coefficient, alpha = deseq_alpha,
      lfcThreshold = deseq_lfc_threshold, altHypothesis = "greaterAbs"
    )
  } else {
    results(dds, name = coefficient, alpha = deseq_alpha)
  }
}

# Test one fitted design, shrink log2FC, write raw and shrunken tables (with a
# `deg` column: up, down or not_significant, from the raw log2FC) and plots.
test_project_deseq2 <- function(project_id, design_type, create_plots) {
  message("Testing ", design_type, " DESeq2 for ", project_id)

  fit_file <- tcga_project_fit_file(project_id, design_type)
  design_file <- tcga_project_design_file(project_id, design_type)
  if (!file.exists(fit_file) || !file.exists(design_file)) {
    return(status_row(project_id, design_type, "skipped", "fitted object not found; run step 04"))
  }

  dds <- readRDS(fit_file)
  design_info <- readRDS(design_file)
  coefficient <- design_info$coefficient
  if (is.null(coefficient)) {
    stop("Design file has no coefficient; rerun step 03")
  }
  if (!coefficient %in% resultsNames(dds)) {
    stop("Coefficient ", coefficient, " not in model: ",
         paste(resultsNames(dds), collapse = ", "))
  }

  bpparam <- make_biocparallel_param()
  parallel_enabled <- use_biocparallel && BiocParallel::bpnworkers(bpparam) > 1

  raw_results <- test_coefficient(dds, coefficient)
  shrunken_results <- if (use_lfc_shrinkage) {
    lfcShrink(
      dds,
      coef = coefficient,
      res = raw_results,
      type = lfc_shrinkage_type,
      parallel = parallel_enabled,
      BPPARAM = bpparam,
      quiet = TRUE
    )
  } else {
    raw_results
  }

  row_data <- as.data.frame(rowData(dds))
  raw_results_df <- annotate_results(as.data.frame(raw_results), row_data)
  shrunken_results_df <- annotate_results(as.data.frame(shrunken_results), row_data)

  raw_results_df$deg <- as.character(classify_deg_direction(raw_results_df))
  shrunken_results_df$deg <- raw_results_df$deg[
    match(shrunken_results_df$gene_id, raw_results_df$gene_id)
  ]

  raw_results_df <- raw_results_df |>
    dplyr::filter(!is.na(padj)) |>
    dplyr::arrange(padj)
  shrunken_results_df <- shrunken_results_df |>
    dplyr::filter(!is.na(padj)) |>
    dplyr::arrange(padj)

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

  if (create_plots) {
    plot_project_ma(shrunken_results_df, project_id)
    plot_project_volcano(shrunken_results_df, project_id)
  }

  status_row(
    project_id, design_type, "completed",
    n_normal = design_info$n_normal,
    n_tumor = design_info$n_tumor,
    coefficient = coefficient,
    deseq_test = deseq_test,
    deseq_alpha = deseq_alpha,
    deseq_lfc_threshold = deseq_lfc_threshold,
    lfc_shrinkage_type = if (use_lfc_shrinkage) lfc_shrinkage_type else NA_character_,
    n_tested = nrow(raw_results_df),
    n_deg = sum(raw_results_df$deg != "not_significant"),
    n_deg_up = sum(raw_results_df$deg == "up"),
    n_deg_down = sum(raw_results_df$deg == "down"),
    covariates = paste(design_info$covariates, collapse = ";"),
    design_formula = format_design_formula(design_info$design_formula)
  )
}

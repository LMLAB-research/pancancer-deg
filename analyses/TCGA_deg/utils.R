#### Shared helpers ####

#### Helper functions ####

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


#### Helper functions ####

format_design_formula <- function(design_formula) {
  gsub("\\s+", " ", paste(deparse(design_formula), collapse = ""))
}

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

classify_deg_direction <- function(results_df) {
  direction <- rep("not_significant", nrow(results_df))
  significant <- !is.na(results_df$padj) &
    results_df$padj < deseq_alpha &
    !is.na(results_df$log2FoldChange) &
    abs(results_df$log2FoldChange) > deseq_lfc_threshold

  direction[significant & results_df$log2FoldChange > 0] <- "up"
  direction[significant & results_df$log2FoldChange < 0] <- "down"

  factor(direction, levels = c("down", "not_significant", "up"))
}

deg_direction_colors <- c(
  down = deg_down_color,
  not_significant = deg_neutral_color,
  up = deg_up_color
)

# Add gene annotation from rowData, filling missing symbols from org.Hs.eg.db.
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

  # Condition is always modeled with levels normal -> tumor.
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

  #### Save Results And QC Outputs ####

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


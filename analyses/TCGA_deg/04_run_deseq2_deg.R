#### Run DESeq2 Differential Expression ####
#
# This script runs projects marked eligible in the design plan, writes raw and
# shrunken DEG tables, and creates first-pass QC plots.

source(file.path("analyses", "TCGA_deg", "00_config.R"))

#### Helper Functions ####

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
  n_groups >= 2 && n_groups <= 12
}

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
    width = 7,
    height = 5,
    dpi = 300
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

  vst_data <- vst(dds, blind = FALSE)

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

plot_project_ma <- function(results_object, project_id) {
  png(tcga_ma_figure_file(project_id), width = 800, height = 600)
  on.exit(dev.off(), add = TRUE)

  plotMA(
    results_object,
    ylim = c(-5, 5),
    main = paste(project_id, "tumor vs normal")
  )
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

  volcano_plot <- EnhancedVolcano(
    plot_df,
    lab = plot_df$gene_name,
    x = "log2FoldChange",
    y = "padj",
    pCutoff = 0.05,
    FCcutoff = 1,
    title = paste(project_id, "tumor vs normal"),
    legendPosition = "bottom"
  )

  ggsave(
    tcga_volcano_figure_file(project_id),
    plot = volcano_plot,
    width = 8,
    height = 8,
    dpi = 300
  )

  TRUE
}

summarise_results <- function(project_id, design_info, results_df, raw_results_df) {
  data.frame(
    project_id = project_id,
    status = "completed",
    reason = NA_character_,
    n_results = nrow(results_df),
    n_raw_results = nrow(raw_results_df),
    n_padj_below_0_05 = sum(results_df$padj < 0.05, na.rm = TRUE),
    n_abs_lfc_above_1 = sum(abs(results_df$log2FoldChange) > 1, na.rm = TRUE),
    n_padj_below_0_05_and_abs_lfc_above_1 = sum(
      results_df$padj < 0.05 & abs(results_df$log2FoldChange) > 1,
      na.rm = TRUE
    ),
    covariates = paste(design_info$covariates, collapse = ";"),
    design_formula = format_design_formula(design_info$design_formula),
    stringsAsFactors = FALSE
  )
}

skip_status <- function(project_id, reason) {
  data.frame(
    project_id = project_id,
    status = "skipped",
    reason = reason,
    n_results = NA_integer_,
    n_raw_results = NA_integer_,
    n_padj_below_0_05 = NA_integer_,
    n_abs_lfc_above_1 = NA_integer_,
    n_padj_below_0_05_and_abs_lfc_above_1 = NA_integer_,
    covariates = NA_character_,
    design_formula = NA_character_,
    stringsAsFactors = FALSE
  )
}

failure_status <- function(project_id, error) {
  data.frame(
    project_id = project_id,
    status = "failed",
    reason = conditionMessage(error),
    n_results = NA_integer_,
    n_raw_results = NA_integer_,
    n_padj_below_0_05 = NA_integer_,
    n_abs_lfc_above_1 = NA_integer_,
    n_padj_below_0_05_and_abs_lfc_above_1 = NA_integer_,
    covariates = NA_character_,
    design_formula = NA_character_,
    stringsAsFactors = FALSE
  )
}

run_project_deseq2 <- function(project_id) {
  message("Running DESeq2 for ", project_id)

  if (!file.exists(tcga_project_dds_file(project_id))) {
    return(skip_status(project_id, "DDS file not found"))
  }

  if (!file.exists(tcga_project_design_file(project_id))) {
    return(skip_status(project_id, "design file not found"))
  }

  dds <- readRDS(tcga_project_dds_file(project_id))
  design_info <- readRDS(tcga_project_design_file(project_id))

  pca_variables <- plot_project_pcas(dds, project_id, design_info)

  dds <- DESeq(dds)

  # Condition is always modeled with levels normal -> tumor.
  contrast <- c("condition", "tumor", "normal")
  raw_results <- results(dds, contrast = contrast, alpha = 0.05)
  shrunken_results <- lfcShrink(
    dds,
    contrast = contrast,
    res = raw_results,
    type = "ashr"
  )

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
    tcga_deseq_raw_results_file(project_id),
    row.names = FALSE
  )
  write.csv(
    shrunken_results_df,
    tcga_deseq_results_file(project_id),
    row.names = FALSE
  )

  plot_project_ma(shrunken_results, project_id)
  volcano_created <- plot_project_volcano(shrunken_results_df, project_id)

  metadata(dds)$tcga_deg <- list(
    pca_variables = pca_variables,
    volcano_created = volcano_created,
    raw_results_file = tcga_deseq_raw_results_file(project_id),
    shrunken_results_file = tcga_deseq_results_file(project_id)
  )
  saveRDS(dds, tcga_project_dds_file(project_id))

  summarise_results(project_id, design_info, shrunken_results_df, raw_results_df)
}

run_project_deseq2_safe <- function(project_id) {
  tryCatch(
    run_project_deseq2(project_id),
    error = function(e) failure_status(project_id, e)
  )
}

#### Select Projects From Design Plan ####

design_plan <- readRDS(tcga_design_plan_rds)
eligible_projects <- design_plan |>
  dplyr::filter(eligible) |>
  dplyr::pull(project_id)
eligible_projects <- limit_projects(eligible_projects)

#### Run Each Eligible Project ####

run_status <- dplyr::bind_rows(lapply(eligible_projects, run_project_deseq2_safe))
write.csv(run_status, tcga_deseq_run_status_file, row.names = FALSE)

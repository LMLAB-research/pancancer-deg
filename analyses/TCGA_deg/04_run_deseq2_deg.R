#### Run DESeq2 Differential Expression ####
#
# This script runs only projects marked eligible in the design plan. It writes
# DEG tables and basic QC plots for each tumor-vs-normal comparison.

source(file.path("analyses", "TCGA_deg", "00_config.R"))

#### Helper Functions ####

# Add gene annotation from rowData, falling back to org.Hs.eg.db when needed.
annotate_results <- function(results_df, row_data) {
  results_df$gene_id <- rownames(results_df)

  if ("gene_id" %in% colnames(row_data)) {
    results_df <- dplyr::left_join(results_df, row_data, by = "gene_id")
  }

  if (!"gene_name" %in% colnames(results_df)) {
    results_df$ensembl_id <- gsub("\\..*", "", results_df$gene_id)
    results_df$gene_name <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = results_df$ensembl_id,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  }

  results_df
}

# PCA on VST data is the first-pass sample-level QC plot.
plot_project_pca <- function(dds, project_id) {
  vst_data <- vst(dds, blind = FALSE)
  pca_plot <- plotPCA(vst_data, intgroup = "condition") +
    theme_minimal() +
    labs(
      title = paste(project_id, "VST PCA"),
      color = "Condition"
    )

  ggsave(
    tcga_pca_figure_file(project_id),
    plot = pca_plot,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# MA plot checks global fold-change behavior after DESeq2 fitting.
plot_project_ma <- function(results_object, project_id) {
  png(tcga_ma_figure_file(project_id), width = 800, height = 600)
  plotMA(
    results_object,
    ylim = c(-5, 5),
    main = paste(project_id, "tumor vs normal")
  )
  dev.off()
}

# Volcano plot is an interpretable overview of significant DEGs.
plot_project_volcano <- function(results_df, project_id) {
  volcano_plot <- EnhancedVolcano(
    results_df,
    lab = results_df$gene_name,
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
}

#### Select Projects From Design Plan ####

design_plan <- readRDS(tcga_design_plan_rds)
eligible_projects <- design_plan |>
  dplyr::filter(eligible) |>
  dplyr::pull(project_id)
eligible_projects <- limit_projects(eligible_projects)

#### Run Each Eligible Project ####

for (project_id in eligible_projects) {
  message("Running DESeq2 for ", project_id)

  dds <- readRDS(tcga_project_dds_file(project_id))

  plot_project_pca(dds, project_id)

  dds <- DESeq(dds)

  # Condition is always modeled with levels normal -> tumor.
  results_object <- results(
    dds,
    contrast = c("condition", "tumor", "normal"),
    alpha = 0.05
  )

  results_df <- as.data.frame(results_object)
  results_df <- annotate_results(results_df, as.data.frame(rowData(dds)))
  results_df <- results_df |>
    dplyr::filter(!is.na(padj)) |>
    dplyr::arrange(padj)

  #### Save Results And QC Outputs ####

  write.csv(results_df, tcga_deseq_results_file(project_id), row.names = FALSE)

  plot_project_ma(results_object, project_id)
  plot_project_volcano(results_df, project_id)

  saveRDS(dds, tcga_project_dds_file(project_id))
}

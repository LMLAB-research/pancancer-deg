# ============================================================
# TCGA differential expression analysis (recount3, DESeq2)
# Example: CHOL (cholangiocarcinoma, bile duct cancer)
# ============================================================

# 0. Load libraries
library(recount3)
library(DESeq2)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)

# 1. Choose TCGA cancer type
cancer_type <- "CHOL"

# 2. Create recount3 cache folder if needed
dir.create(
  "C:/Users/danab/AppData/Local/R/cache/R/recount3",
  recursive = TRUE,
  showWarnings = FALSE
)

# 3. Find TCGA project in recount3
projects <- available_projects()

tcga_project <- subset(
  projects,
  project == cancer_type &
    project_home == "data_sources/tcga" &
    organism == "human"
)

# 4. Download gene-level RangedSummarizedExperiment object
rse <- create_rse(
  project_info = tcga_project,
  type = "gene",
  annotation = "gencode_v29"
)

# 5. Define tumor / normal condition
sample_type_col <- "tcga.gdc_cases.samples.sample_type"

colData(rse)$condition <- ifelse(
  colData(rse)[[sample_type_col]] == "Solid Tissue Normal",
  "normal",
  "tumor"
)

# TCGA sample types include categories such as:
# "Primary Tumor", "Metastatic", "Recurrent Tumor",
# and "Solid Tissue Normal".
#
# Here we define:
# - normal samples = "Solid Tissue Normal"
# - all other sample types = tumor
#
# For a stricter analysis, it may be better to keep only
# "Primary Tumor" and "Solid Tissue Normal" samples.

colData(rse)$condition <- factor(
  colData(rse)$condition,
  levels = c("normal", "tumor")
)

table(colData(rse)$condition)

# 6. Keep only tumor and normal samples
rse_sub <- rse[, colData(rse)$condition %in% c("normal", "tumor")]

# 7. Extract transformed counts and metadata
counts_mat <- transform_counts(rse_sub)
counts_mat <- round(counts_mat)

metadata <- as.data.frame(colData(rse_sub))

dim(counts_mat)
dim(metadata)

# 7b. Inspect input data
counts_mat[1:10, 1:5]
View(counts_mat[1:500, 1:10])
View(metadata)

table(metadata$condition)

gene_info <- as.data.frame(rowData(rse_sub))
View(gene_info)

library_sizes <- colSums(counts_mat)

barplot(
  library_sizes,
  las = 2,
  main = paste("Library sizes:", cancer_type),
  ylab = "Total counts per sample"
)

# 8. Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = metadata,
  design = ~ condition
)

# 9. Filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= 3

cat("Genes before filtering:", length(keep), "\n")

dds <- dds[keep, ]

cat("Genes after filtering:", nrow(dds), "\n")

# 10. Run differential expression analysis
dds <- DESeq(dds)

# 10b. QC plots
plotDispEsts(dds)

res <- results(
  dds,
  contrast = c("condition", "tumor", "normal")
)

plotMA(
  res,
  main = paste("MA plot:", cancer_type),
  ylim = c(-8, 8)
)

# 10c. VST normalization and PCA
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)

percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = condition
  )
) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = paste("PCA:", cancer_type),
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  )

pca_plot

ggsave(
  paste0("TCGA_", cancer_type, "_PCA.png"),
  plot = pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# 11. Convert DESeq2 results to data frame
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Remove genes without adjusted p-value
res_df <- res_df[!is.na(res_df$padj), ]

# 12. Add gene annotation
res_df$ensembl_id <- gsub("\\..*", "", res_df$gene_id)

res_df$symbol <- mapIds(
  org.Hs.eg.db,
  keys = res_df$ensembl_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Reorder columns
res_df <- res_df[, c(
  "gene_id", "ensembl_id", "symbol",
  "baseMean", "log2FoldChange", "lfcSE",
  "stat", "pvalue", "padj"
)]

# Sort by adjusted p-value
res_df <- res_df[order(res_df$padj), ]

# 13. Extract upregulated and downregulated genes
upregulated <- subset(
  res_df,
  padj < 0.05 & log2FoldChange > 1
)

downregulated <- subset(
  res_df,
  padj < 0.05 & log2FoldChange < -1
)

# 14. Show top genes
head(upregulated, 20)
head(downregulated, 20)

# 15. Save result tables
write.csv(
  res_df,
  paste0("TCGA_", cancer_type, "_DE_results_annotated.csv"),
  row.names = FALSE
)

write.csv(
  upregulated,
  paste0("TCGA_", cancer_type, "_upregulated_annotated.csv"),
  row.names = FALSE
)

write.csv(
  downregulated,
  paste0("TCGA_", cancer_type, "_downregulated_annotated.csv"),
  row.names = FALSE
)

# 16. Print summary
cat("Cancer type:", cancer_type, "\n")
cat("Samples:", ncol(counts_mat), "\n")
cat("Genes after filtering:", nrow(res_df), "\n")
cat("Upregulated genes:", nrow(upregulated), "\n")
cat("Downregulated genes:", nrow(downregulated), "\n")

# 17. Volcano plot
res_df$significance <- "Not significant"

res_df$significance[
  res_df$padj < 0.05 & res_df$log2FoldChange > 1
] <- "Upregulated"

res_df$significance[
  res_df$padj < 0.05 & res_df$log2FoldChange < -1
] <- "Downregulated"

volcano_plot <- ggplot(
  res_df,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = significance
  )
) +
  geom_point(alpha = 0.6, size = 1.2) +
  theme_minimal() +
  labs(
    title = paste("Volcano plot:", cancer_type),
    x = "log2 fold change: tumor vs normal",
    y = "-log10 adjusted p-value"
  )

volcano_plot


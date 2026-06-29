if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ashr")
library(ashr)

library(DESeq2)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano) # Optional, for a high-quality volcano plot
library(pheatmap)        # Optional, for sample-to-sample distance heatmaps

rse <- readRDS("TCGA-KIRC_DGE_ready_RSE.rds")
rse_complete <- rse[, colData(rse)$has_complete_metadata]
# 1b. CRITICAL FIX: Ensure no NAs in sample_type_f, and subset to our groups
valid_samples <- !is.na(colData(rse_complete)$sample_type_f) & 
  colData(rse_complete)$sample_type_f %in% c("Solid Tissue Normal", "Primary Tumor")

rse_complete <- rse_complete[, valid_samples]

message(sprintf("Proceeding with %d samples with complete metadata.", ncol(rse_complete)))
print(table(colData(rse_complete)$sample_type))


# Centering age to ensure GLM convergence
colData(rse_complete)$age_centered <- as.numeric(scale(colData(rse_complete)$age_at_diagnosis, scale = FALSE))

#  FIX: Use the proper factor columns
dds <- DESeqDataSet(rse_complete, design = ~ sex_at_birth_f + age_centered + sample_type_f)
# 3. Best Practice Quality Control & Filtering
# -----------------------------------------------------------------------------
# Let's keep genes with at least 10 counts in 
# at least as many samples as the smallest group size 
smallest_group_size <- min(table(colData(dds)$sample_type_f))
keep <- rowSums(counts(dds) >= 10) >= smallest_group_size
dds <- dds[keep,]
dds
message(sprintf("Filtering retained %d genes out of %d.", nrow(dds), length(keep)))


# 4. Exploratory Data Analysis (PCA)
# -----------------------------------------------------------------------------
# Use VST (Variance Stabilizing Transformation) for visualization/clustering.
# blind = FALSE ensures it accounts for our design variables.
vst_data <- vst(dds, blind = FALSE)

# Generate PCA Plot
pca_plot <- plotPCA(vst_data, intgroup = c("sample_type_f")) +
  theme_minimal() +
  labs(title = "PCA of TCGA-KIRC (VST transformed)", color = "Sample Type")
ggsave("TCGA-KIRC_PCA_plot.png", plot = pca_plot, width = 7, height = 5)

# 5. Run the DESeq2 Differential Expression Analysis
# -----------------------------------------------------------------------------
message("Running DESeq2 pipeline...")
dds <- DESeq(dds)


# 6. Extract Results & Apply Shrinkage (Crucial Best Practice)
# -----------------------------------------------------------------------------
# List the coefficients to ensure we extract the right contrast
print(resultsNames(dds))

# Extract the raw results for Tumor vs Normal
res_raw <- results(dds, 
                   contrast = c("sample_type_f", "Primary Tumor", "Solid Tissue Normal"),
                   alpha = 0.05) # Target FDR

# Apply lfcShrink (apeglm method is highly recommended)
# This reduces noise/false positives in low-expression genes with high variance
res_shrunk <- lfcShrink(dds, 
                        coef = "sample_type_f_Primary.Tumor_vs_Solid.Tissue.Normal", 
                       type = "ashr")

# 7. Annotate & Save Results
# -----------------------------------------------------------------------------
# Merge DESeq2 results with the gene annotations already present in your rowData
gene_info <- as.data.frame(rowData(dds))

res_df <- as.data.frame(res_shrunk) %>%
  mutate(gene_id = rownames(.)) %>%
  left_join(gene_info, by = "gene_id")


write.csv(res_df, "TCGA-KIRC_DGE_Results_Tumor_vs_Normal.csv", row.names = FALSE)
message("Saved DGE results to 'TCGA-KIRC_DGE_Results_Tumor_vs_Normal.csv'")

png("TCGA-KIRC_MA_plot.png", width = 800, height = 600)
plotMA(res_shrunk, ylim = c(-5, 5), main = "MA Plot (Shrunk LFC)")
dev.off()

# 8b. Volcano Plot using EnhancedVolcano
volcano <- EnhancedVolcano(res_df,
                           lab = res_df$gene_name,
                           x = 'log2FoldChange',
                           y = 'padj',
                           pCutoff = 0.05,
                           FCcutoff = 1.5,
                           title = 'Primary Tumor vs Solid Tissue Normal',
                           subtitle = 'Controlled for Sex and Age',
                           legendPosition = 'bottom')
ggsave("TCGA-KIRC_Volcano_plot.png", plot = volcano, width = 8, height = 8)

message("Pipeline complete!")
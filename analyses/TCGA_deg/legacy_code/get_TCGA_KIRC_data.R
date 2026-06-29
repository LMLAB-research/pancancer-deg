Sys.setenv(R_MAX_VSIZE = 32e9)

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(readr)

proj <- "TCGA-KIRC"

message("=======================================================")
message(" Processing: ", proj)
message("=======================================================")


# -----------------------------------------------------------------------------
# 1. Query, Download & Build RSE
# -----------------------------------------------------------------------------
rna_query <- GDCquery(
  project               = proj,
  data.category         = "Transcriptome Profiling",
  data.type             = "Gene Expression Quantification",
  workflow.type         = "STAR - Counts",
  experimental.strategy = "RNA-Seq"
)

GDCdownload(rna_query, method = "api", files.per.chunk = 20)

rse <- GDCprepare(rna_query, summarizedExperiment = TRUE)

message("RSE built — dims: ", nrow(rse), " genes x ", ncol(rse), " samples")


# -----------------------------------------------------------------------------
# 2. Slim down assays — keep only raw unstranded counts for DGE
# -----------------------------------------------------------------------------
assays(rse) <- list(counts = assay(rse, "unstranded"))


# -----------------------------------------------------------------------------
# 3. Add metadata factors to colData
# -----------------------------------------------------------------------------
cd <- as.data.frame(colData(rse))

# Factor for sample_type (reference level = Solid Tissue Normal)
colData(rse)$sample_type_f <- factor(
  cd$sample_type,
  levels = c("Solid Tissue Normal", "Primary Tumor", "Metastatic")
)

# Factor for sex_at_birth
colData(rse)$sex_at_birth_f <- factor(
  toupper(cd$sex_at_birth),
  levels = c("FEMALE", "MALE")
)

# Numeric age at diagnosis
colData(rse)$age_at_diagnosis <- suppressWarnings(as.numeric(cd$age_at_diagnosis))

# Completeness flag
colData(rse)$has_complete_metadata <-
  !is.na(colData(rse)$age_at_diagnosis) &
  !is.na(colData(rse)$sex_at_birth_f)

message(sprintf(
  "Samples — Primary Tumor: %d | Solid Tissue Normal: %d | complete metadata: %d",
  sum(cd$sample_type == "Primary Tumor",       na.rm = TRUE),
  sum(cd$sample_type == "Solid Tissue Normal", na.rm = TRUE),
  sum(colData(rse)$has_complete_metadata,      na.rm = TRUE)
))


# -----------------------------------------------------------------------------
# 4. Trim rowData to essential gene annotation columns
# -----------------------------------------------------------------------------
keep_row_cols <- intersect(
  c("gene_id", "gene_name", "gene_type", "seqnames", "start", "end", "strand"),
  colnames(rowData(rse))
)
rowData(rse) <- rowData(rse)[, keep_row_cols, drop = FALSE]


# -----------------------------------------------------------------------------
# 5. Pre-filter: remove genes with zero counts across all samples
# -----------------------------------------------------------------------------
keep_genes   <- rowSums(assay(rse, "counts")) > 0
rse_filtered <- rse[keep_genes, ]
message(sprintf(
  "Gene filtering: %d → %d genes (removed %d all-zero rows)",
  nrow(rse), nrow(rse_filtered), nrow(rse) - nrow(rse_filtered)
))
rse <- rse_filtered
rm(rse_filtered)


# -----------------------------------------------------------------------------
# 6. Export
# -----------------------------------------------------------------------------

# 6a. Main DGE-ready RSE object
saveRDS(rse, file = "TCGA-KIRC_DGE_ready_RSE.rds")
message("Saved: TCGA-KIRC_DGE_ready_RSE.rds")

# 6b. Flat metadata CSV
clean_meta_export <- as.data.frame(colData(rse)) %>%
  select(any_of(c(
    "barcode", "sample_type", "sex_at_birth", "age_at_diagnosis",
    "has_complete_metadata", "tumor_stage", "primary_diagnosis",
    "tissue_or_organ_of_origin"
  )))

write.csv(clean_meta_export, file = "TCGA-KIRC_metadata.csv", row.names = FALSE)
message("Saved: TCGA-KIRC_metadata.csv")

# 6c. Quick summary
message("\n--- Metadata summary ---")
print(table(SampleType = colData(rse)$sample_type,
            Sex        = colData(rse)$sex_at_birth_f))
message(sprintf("Age range: %.1f – %.1f years (median %.1f)",
                min(colData(rse)$age_at_diagnosis, na.rm = TRUE),
                max(colData(rse)$age_at_diagnosis, na.rm = TRUE),
                median(colData(rse)$age_at_diagnosis, na.rm = TRUE)
))

gc()
message("\nDone. RSE is ready for DESeq2 / edgeR.")
# Per-tumour scores for the RAS-induced-senescence (RIS) gene classes, ERK output, MYC
# targets, proliferation and a stromal proxy, in every local TCGA cohort; joined to
# RAS-hotspot / MYC-amplification calls. Exploratory; run from the pancancer-deg root.
eval(parse(".claude/worktrees/two-sided-fisher/.Rprofile")[[2]])  # sandbox proxy fix (msigdbr download)
suppressMessages({ library(SummarizedExperiment); library(DESeq2) })
OUT <- "results/exploratory/ras_vs_myc"
cls <- read.csv("results/exploratory/ras_vs_ois/strict/ois_gene_classes.csv")
drv <- read.csv(file.path(OUT, "driver_status.csv"))
drv$sample15 <- substr(drv$sample, 1, 15)

# Hallmark sets via msigdbr from the senescence_stress bulk library (read-only)
local({
  .libPaths(c("/Users/ro/Academic/10_Research/senescence_stress/modules/systems_analysis_arrest/bulk_RNA_seq/.Rlib", .libPaths()))
  h <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  assign("hs", split(h$gene_symbol, h$gs_name), envir = globalenv())
})
stroma <- c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A3", "PDGFRB", "FAP", "THY1", "SPARC", "LUM")
sets <- list(
  RIS_RASauto_UP = cls$gene[cls$OIS_direction == "UP" & cls$arrest_class == "RAS-autonomous"],
  RIS_arrest_UP = cls$gene[cls$OIS_direction == "UP" & cls$arrest_class == "arrest-coupled"],
  RIS_intermediate_UP = cls$gene[cls$OIS_direction == "UP" & cls$arrest_class == "intermediate"],
  RIS_NFKB_BRD4_UP = cls$gene[cls$OIS_direction == "UP" & cls$nfkb_brd4 == "dependent"],
  RIS_RASauto_DOWN = cls$gene[cls$OIS_direction == "DOWN" & cls$arrest_class == "RAS-autonomous"],
  ERK_output = c("DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "FOSL1", "PHLDA1", "SPRED1", "SPRED2", "DUSP4"),
  MYC_targets = unique(c(hs$HALLMARK_MYC_TARGETS_V1, hs$HALLMARK_MYC_TARGETS_V2)),
  Proliferation = unique(c(hs$HALLMARK_E2F_TARGETS, hs$HALLMARK_G2M_CHECKPOINT)),
  Stroma = stroma
)
# keep the class sets free of the stromal proxy genes and of each other's overlaps with MYC/prolif
for (k in grep("^RIS_", names(sets), value = TRUE)) sets[[k]] <- setdiff(sets[[k]], stroma)

files <- c(list.files("data/TCGA_deg/processed", pattern = "_DGE_ready_RSE.rds$", full.names = TRUE),
           # tumour-only cohorts downloaded for the driver analysis (07_download_paad_skcm.R)
           list.files("data/exploratory/tcga_extra", pattern = "_STAR_counts_RSE.rds$", full.names = TRUE))
all <- do.call(rbind, lapply(files, function(f) {
  r <- readRDS(f)
  # primary tumours; for SKCM (mostly metastases in TCGA) also metastatic samples, with
  # sample type kept as a covariate (stype) for the models
  keep_types <- if (grepl("SKCM", f)) c("TP", "TM") else "TP"
  is_tumour <- if ("condition" %in% colnames(colData(r))) r$condition == "tumor" else TRUE
  r <- r[rowData(r)$gene_type == "protein_coding", is_tumour & r$shortLetterCode %in% keep_types]
  if (!"unstranded" %in% assayNames(r) && !"counts" %in% assayNames(r)) stop("no counts assay in ", f)
  if (!"counts" %in% assayNames(r)) assay(r, "counts") <- assay(r, "unstranded")
  m <- assay(r, "counts"); storage.mode(m) <- "integer"
  m <- rowsum(m, rowData(r)$gene_name)
  m <- m[rowMeans(m >= 10) >= 0.2, ]
  v <- assay(vst(DESeqDataSetFromMatrix(m, data.frame(row.names = colnames(m), x = rep(1, ncol(m))), ~1), blind = TRUE))
  z <- t(scale(t(v)))
  sc <- sapply(sets, function(g) colMeans(z[intersect(g, rownames(z)), , drop = FALSE], na.rm = TRUE))
  data.frame(cohort = sub("_(DGE|STAR).*", "", sub("TCGA-", "", basename(f))), barcode = r$barcode,
             sample15 = substr(r$barcode, 1, 15), stype = r$shortLetterCode, sc, check.names = FALSE)
}))
all$cohort_cbio <- ifelse(all$cohort %in% c("COAD", "READ"), "COADREAD", all$cohort)
all <- merge(all, drv[, c("sample15", "ras_hotspot", "ras_mut", "myc_amp", "cna_profiled", "sequenced")],
             by = "sample15", all.x = TRUE)
write.csv(all, file.path(OUT, "tumour_scores.csv"), row.names = FALSE)
cat("tumours scored:", nrow(all), " with driver calls:", sum(!is.na(all$ras_mut)), "\n")
print(table(all$cohort, ifelse(is.na(all$ras_mut), "no_call", paste0("RAS=", all$ras_mut, " MYCamp=", all$myc_amp))))

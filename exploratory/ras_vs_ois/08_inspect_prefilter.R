# Read-only: inspect GSE245807 house objects before the LMNB1/G1 filter.
source("exploratory/paths.R")   # repo-relative inputs; run from the pancancer-deg repo root
add_lib(SC_LIB)
suppressMessages(library(Seurat))
H <- file.path(SS_DIR, "modules/scrna_seq/cs_dataset_survey")
t <- readRDS(file.path(H, "house_rerun/GSE245807/treated_clean_rebuilt_seurat_objs.rds"))
for (p in names(t)) {
  o <- t[[p]]; cat("==", p, ncol(o), "cells; assays", Assays(o), "\n")
  cat("  cols:", intersect(c("Phase", "LMNB1", "ctrl_like", "frac_CTRL1_nn", "nCount_RNA", "percent.mt"), colnames(o[[]])), "\n")
  if ("Phase" %in% colnames(o[[]])) print(table(o$Phase, LMNB1pos = if ("LMNB1" %in% colnames(o[[]])) o$LMNB1 > 0 else NA))
  cat("  barcodes:", head(colnames(o), 2), "\n")
}
a <- readRDS(file.path(H, "analysis/GSE245807/conditions_oneplate_seurat.rds"))
cat("== all-cell object", ncol(a), "cells; assays", Assays(a), "; barcodes:", head(colnames(a), 2), "\n")
cat("  cols:", head(colnames(a[[]]), 30), "\n")

# Read-only inspection of the house GSE245807 RIS final object (senescence_stress).
source("exploratory/paths.R")   # repo-relative inputs; run from the pancancer-deg repo root
add_lib(SC_LIB)
suppressMessages(library(Seurat))
f <- file.path(SS_DIR, "modules/scrna_seq/cs_dataset_survey/house_rerun/GSE245807/final/RIS_seurat.rds")
o <- readRDS(f)
print(o)
md <- o@meta.data
cat("meta columns:", colnames(md), "\n")
for (k in colnames(md)) {
  v <- md[[k]]
  if ((is.character(v) || is.factor(v)) && length(unique(v)) <= 15) { cat("\n", k, ":\n"); print(table(v)) }
}
cat("assays:", Assays(o), "\n")

# GSE173879 (Evans 2023, LF1 ER:HRAS clone 1): OIS cells kept in 4-OHT (OISON) vs OIS cells moved to
# DMSO before harvest (OISOFF, RAS withdrawn after senescence). House-final cells from the house
# one-plate preprocessing (all house filters applied; outputs in data/exploratory/house_GSE173879_OIS).
# Prediction: RAS-autonomous genes revert in OISOFF, arrest-coupled genes persist. UCell via scutils.
.libPaths(c("/Users/ro/Academic/10_Research/senescence_stress/modules/scrna_seq/GSE226225/final_analysis/.Rlib/main", .libPaths()))
suppressMessages({ library(Seurat); library(UCell) })
source("/Users/ro/code/academic/scutils/R/ucell_parallel.R")
D <- "data/exploratory/house_GSE173879_OIS/"
out <- "results/exploratory/ras_vs_ois/scrna_GSE173879"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cls <- read.csv("results/exploratory/ras_vs_ois/strict/ois_gene_classes.csv")
sel <- function(d, c) cls$gene[cls$OIS_direction == d & cls$arrest_class == c]
sets0 <- list(MAPK_output = c("DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "FOSL1", "EGR1"),
              RASauto_UP = sel("UP", "RAS-autonomous"), intermediate_UP = sel("UP", "intermediate"),
              arrest_UP = sel("UP", "arrest-coupled"),
              NFKB_BRD4_dep_UP = cls$gene[cls$OIS_direction == "UP" & cls$nfkb_brd4 == "dependent"],
              RASauto_DOWN = sel("DOWN", "RAS-autonomous"), arrest_DOWN = sel("DOWN", "arrest-coupled"))
source(".claude/worktrees/two-sided-fisher/exploratory/ras_vs_ois/pathway_sets.R"); sets0 <- c(sets0, pathway_sets())  # added 2026-11-07: pathway scores (TGF-beta, cGAS-STING, IFN, ...)
score <- function(o, lab) {
  DefaultAssay(o) <- "RNA"
  s <- lapply(sets0, intersect, rownames(o))
  o <- score_ucell_parallel(o, signatures = s, ncores = 1, assay = "RNA", slot = "counts")
  m <- o[[]][, grep("_UCell$", colnames(o[[]])), drop = FALSE]; colnames(m) <- sub("_UCell$", "", colnames(m))
  data.frame(group = lab, m)
}
fin <- list(OISON = readRDS(paste0(D, "treated_clean_rebuilt_seurat_objs.rds"))[["OISON_1"]],
            OISOFF = readRDS(paste0(D, "treated_clean_rebuilt_seurat_objs.rds"))[["OISOFF_1"]])
raw <- readRDS(paste0(D, "raw_seurat_objs.rds"))
ctrl <- subset(raw[[1]]$merged, cells = WhichCells(raw[[1]]$merged, expression = accession == "CTRL_1"))
rm(raw); gc()
sc <- rbind(score(ctrl, "Proliferating (CTRL)"), score(fin$OISON, "OIS, 4-OHT kept (house final)"),
            score(fin$OISOFF, "OIS, 4-OHT withdrawn (house final)"))
write.csv(sc, file.path(out, "withdrawal_cell_scores.csv"))
cat("cells:\n"); print(table(sc$group))
med <- t(sapply(split(sc[, -1], sc$group), function(x) sapply(x, median)))
print(round(med, 4))
on <- med["OIS, 4-OHT kept (house final)", ]; off <- med["OIS, 4-OHT withdrawn (house final)", ]; c0 <- med["Proliferating (CTRL)", ]
cat("\nFraction of the OIS change (vs proliferating) remaining after 4-OHT withdrawal:\n")
print(round((off - c0) / (on - c0), 2))
p <- sapply(names(sets0), function(k) wilcox.test(sc[[k]][sc$group == "OIS, 4-OHT withdrawn (house final)"],
                                                    sc[[k]][sc$group == "OIS, 4-OHT kept (house final)"])$p.value)
cat("\nWilcoxon p, withdrawn vs kept (BH):\n"); print(signif(p.adjust(p, "BH"), 2))

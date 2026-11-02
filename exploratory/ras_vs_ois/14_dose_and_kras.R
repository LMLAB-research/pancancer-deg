# RIS gene classes against RAS dose and early KRAS time course.
# (a) GSE222951 (Chan et al., "Titration of RAS alters senescent state"): RPE1-hTERT ER:HRAS-G12V, 4-OHT 6 days,
#     flow-sorted by RAS-reporter level S < M < L < XL (hashtagged monocultures, 10x). No untreated monoculture, so
#     this tests dose within RAS-on cells. Cell calling: CLR-normalised HTO, HTODemux singlets (positive.quantile
#     0.99), as in the house staging script.
# (b) GSE202126: HPNE (hTERT pancreatic duct) parental vs KRAS days 1, 2, 3, 5 (VASA-seq single cells, exon UMIs).
# QC: house default (Wechter) cutoffs only: genes >= 300, UMI 1,000-120,000, mitochondrial < 12 %. The house
# senescent-state filters (control-like kNN, LMNB1, G1) are NOT applied: this is a dose / time-course comparison,
# not a definition of a senescent state, and those filters would remove the cells of interest.
# Scores: per-cell UCell (scutils fork-based helper, serial) on the script-13 sets plus proliferation (capped at
# 1,400 genes for maxRank 1,500; the list has 671 genes, so the cap does not bind).
# Statistics: median per group; Spearman rho with group order; fraction of the day-5 change reached per day;
# Wilcoxon (scores) or Fisher (detection fractions) tests, XL vs S and day 5 vs parental, BH within dataset.
# Run from the pancancer-deg repo root with Rscript --vanilla.
.libPaths(c("/Users/ro/Academic/10_Research/senescence_stress/modules/scrna_seq/GSE226225/final_analysis/.Rlib/main", .libPaths()))
suppressMessages({ library(Seurat); library(UCell); library(Matrix) })
source("/Users/ro/code/academic/scutils/R/ucell_parallel.R")
set.seed(1)
out <- "results/exploratory/ras_vs_ois/dose_and_kras"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cls <- read.csv("results/exploratory/ras_vs_ois/strict/ois_gene_classes.csv")
prolif <- readLines("results/exploratory/two_sided_2026-11-01/prolif_genes.txt")
sel <- function(d, c) cls$gene[cls$OIS_direction == d & cls$arrest_class == c]
sets0 <- list(MAPK_output = c("DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "FOSL1", "EGR1"),
              Proliferation = prolif,
              RASauto_UP = sel("UP", "RAS-autonomous"), intermediate_UP = sel("UP", "intermediate"),
              arrest_UP = sel("UP", "arrest-coupled"),
              NFKB_BRD4_dep_UP = cls$gene[cls$OIS_direction == "UP" & cls$nfkb_brd4 == "dependent"],
              RASauto_DOWN = sel("DOWN", "RAS-autonomous"), arrest_DOWN = sel("DOWN", "arrest-coupled"))
FRAC <- c("MKI67_pos", "CDKN1A_pos")   # detection fractions; MKI67/CDKN1A log1p(CPM) are reported as medians

qc <- function(o, tag, umi_max = 120000) {
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
  n0 <- table(o$group)
  o <- subset(o, subset = nFeature_RNA >= 300 & nCount_RNA >= 1000 & nCount_RNA <= umi_max & percent.mt < 12)
  qcn <- data.frame(dataset = tag, group = names(n0), before_qc = as.integer(n0),
                    after_qc = as.integer(table(factor(o$group, levels = names(n0)))))
  list(o = o, n = qcn)
}

analyse <- function(o, lev, tag, ref, top) {
  s <- lapply(sets0, intersect, rownames(o)); s$Proliferation <- head(s$Proliferation, 1400)
  cat("\n==", tag, "genes found per set:\n"); print(sapply(s, length))
  o <- score_ucell_parallel(o, signatures = s, ncores = 1, assay = "RNA", slot = "counts")
  sc <- o[[]][, grep("_UCell$", colnames(o[[]])), drop = FALSE]; colnames(sc) <- sub("_UCell$", "", colnames(sc))
  cnt <- GetAssayData(o, assay = "RNA", layer = "counts")
  sc$MKI67_pos <- as.numeric(cnt["MKI67", ] > 0); sc$CDKN1A_pos <- as.numeric(cnt["CDKN1A", ] > 0)
  lib <- Matrix::colSums(cnt)
  sc$MKI67_logCPM <- log1p(cnt["MKI67", ] / lib * 1e6); sc$CDKN1A_logCPM <- log1p(cnt["CDKN1A", ] / lib * 1e6)
  sc$group <- factor(o$group, levels = lev)
  write.csv(sc, file.path(out, paste0(tag, "_cell_scores.csv")))
  vars <- setdiff(colnames(sc), "group")
  med <- sapply(vars, function(k) tapply(sc[[k]], sc$group, if (k %in% FRAC) mean else median))
  rho <- sapply(vars, function(k) cor(as.integer(sc$group), sc[[k]], method = "spearman"))
  p <- sapply(vars, function(k) {
    a <- sc[[k]][sc$group == top]; b <- sc[[k]][sc$group == ref]
    if (k %in% FRAC) fisher.test(table(factor(c(a, b), 0:1), rep(1:2, c(length(a), length(b)))))$p.value
    else wilcox.test(a, b)$p.value })
  tab <- data.frame(score = vars, t(med), check.names = FALSE)
  tab$spearman_rho <- rho; tab$delta <- med[top, ] - med[ref, ]
  tab$p <- p; tab$padj_BH <- p.adjust(p, "BH"); tab$test <- ifelse(vars %in% FRAC, "Fisher", "Wilcoxon")
  colnames(tab)[colnames(tab) == "delta"] <- paste0("delta_", top, "_vs_", ref)
  cat("\n==", tag, "cells per group:\n"); print(table(sc$group))
  cat("Median score per group (MKI67_pos/CDKN1A_pos = fraction of cells detected), Spearman rho with group order,",
      top, "vs", ref, "(BH):\n")
  print(format(tab[, -1], digits = 3), quote = FALSE)
  list(sc = sc, tab = tab, med = med)
}

# (a) GSE222951 RAS dose
R <- "data/exploratory/GSE222951/raw/"
m <- readMM(paste0(R, "GSM7948808_rna_matrix.mtx.gz"))
f <- read.delim(gzfile(paste0(R, "GSM7948808_rna_features.tsv.gz")), header = FALSE)
colnames(m) <- sub("-1$", "", readLines(gzfile(paste0(R, "GSM7948808_rna_barcodes.tsv.gz"))))   # HTO barcodes carry no "-1"
m <- as(m, "CsparseMatrix"); rownames(m) <- make.unique(f$V2)   # symbols; duplicates are made unique, as Read10X does
h <- readMM(paste0(R, "GSM7948809_hto_matrix.mtx.gz"))
rownames(h) <- readLines(gzfile(paste0(R, "GSM7948809_hto_features.tsv.gz")))
colnames(h) <- readLines(gzfile(paste0(R, "GSM7948809_hto_barcodes.tsv.gz")))
h <- as(h[rownames(h) != "unmapped", ], "CsparseMatrix")
cells <- intersect(colnames(m), colnames(h))
cat("GSE222951 barcodes: RNA", ncol(m), "HTO", ncol(h), "shared", length(cells), "\n")
o <- CreateSeuratObject(counts = m[, cells]); o[["HTO"]] <- CreateAssayObject(counts = h[, cells])
o <- NormalizeData(o, assay = "HTO", normalization.method = "CLR", verbose = FALSE)
o <- HTODemux(o, assay = "HTO", positive.quantile = 0.99, verbose = FALSE)
cat("HTODemux global classification:\n"); print(table(o$HTO_classification.global))
o <- subset(o, subset = HTO_classification.global == "Singlet")
cat("Singlets per hashtag:\n"); print(table(o$hash.ID))
o$group <- sub("^monoculture-([A-Z]+)-d6.*", "\\1", gsub("_", "-", as.character(o$hash.ID)))
dose_lev <- c("S", "M", "L", "XL")
o <- subset(o, cells = colnames(o)[o$group %in% dose_lev])
q1 <- qc(o, "GSE222951"); rm(o, m, h); gc(verbose = FALSE)
A <- analyse(q1$o, dose_lev, "GSE222951_RAS_dose", ref = "S", top = "XL")
write.csv(A$tab, file.path(out, "GSE222951_RAS_dose_summary.csv"), row.names = FALSE)
q1$o <- NULL; gc(verbose = FALSE)

# (b) GSE202126 HPNE KRAS time course
G <- "data/exploratory/GSE202126/"
rdm <- function(f) {
  x <- read.delim(gzfile(f), check.names = FALSE)
  as(rowsum(as.matrix(x[, -(1:2)]), x$gene_symbol), "CsparseMatrix")   # sum duplicated symbols
}
a <- rdm(paste0(G, "GSE202126_HPNE_exon_UMI_count_matrix.txt.gz"))
b <- rdm(paste0(G, "GSE202126_HPNE_KRAS_exon_UMI_count_matrix.txt.gz"))
g <- intersect(rownames(a), rownames(b)); mk <- cbind(a[g, ], b[g, ]); rm(a, b)
o <- CreateSeuratObject(counts = mk)
o$group <- ifelse(grepl("^HPNE_Kras_day", colnames(o)), sub("^HPNE_Kras_(day[0-9]+)_.*", "\\1", colnames(o)), "parental")
tc_lev <- c("parental", "day1", "day2", "day3", "day5")
cat("\nGSE202126 cells per group before QC:\n"); print(table(factor(o$group, tc_lev)))
cat("UMI / genes / mt% quantiles before QC:\n")
o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
print(sapply(o[[c("nCount_RNA", "nFeature_RNA", "percent.mt")]], quantile, c(.05, .25, .5, .75, .95)))
# Primary: house default QC as for the 10x data. Sensitivity: VASA-seq is plate-based full-length total RNA with a
# median of ~156k UMI per cell, so the droplet-oriented 120,000 UMI upper cap removes most cells; rerun without it.
tc <- function(q, tag) {
  B <- analyse(q$o, tc_lev, tag, ref = "parental", top = "day5")
  fr <- t(sapply(c("day1", "day2", "day3", "day5"), function(d) (B$med[d, ] - B$med["parental", ]) / (B$med["day5", ] - B$med["parental", ])))
  cat("\nFraction of the day-5 change (vs parental) reached at each day:\n"); print(round(fr, 2))
  tabB <- B$tab; for (d in rownames(fr)) tabB[[paste0("frac_of_day5_change_", d)]] <- fr[d, ]
  write.csv(tabB, file.path(out, paste0(tag, "_summary.csv")), row.names = FALSE)
}
q2 <- qc(o, "GSE202126")
tc(q2, "GSE202126_HPNE_KRAS_timecourse")
q3 <- qc(o, "GSE202126_no_UMI_cap", umi_max = Inf); rm(o, mk); gc(verbose = FALSE)
tc(q3, "GSE202126_HPNE_KRAS_timecourse_noUMIcap")

write.csv(rbind(q1$n, q2$n, q3$n), file.path(out, "cells_per_group_qc.csv"), row.names = FALSE)
writeLines(c("QC: house default (Wechter) cutoffs only: genes >= 300, UMI 1,000-120,000, mitochondrial < 12 %.",
             "House senescent-state filters (control-like kNN, LMNB1, G1) were NOT applied: dose / time-course comparison,",
             "not a senescent-state definition; those filters would remove the cells of interest.",
             "GSE222951: CLR-normalised HTO, HTODemux singlets (positive.quantile 0.99), S/M/L/XL monocultures only.",
             "GSE222951 has no untreated monoculture: it tests RAS dose within RAS-on cells.",
             "GSE202126 (VASA-seq, plate-based, median ~156k UMI/cell): primary run uses the house cutoffs; the",
             "120,000 UMI upper cap removes most cells, so a sensitivity run without the upper cap (*_noUMIcap) is also given.",
             "MKI67 and CDKN1A are detected in every VASA-seq cell, so their log1p(CPM) medians are also reported.",
             "Tests: Wilcoxon for UCell scores, Fisher for detection fractions; BH within dataset. set.seed(1); serial UCell."),
           file.path(out, "README_methods.txt"))
cat("\nCells per group (QC):\n"); print(rbind(q1$n, q2$n, q3$n))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))

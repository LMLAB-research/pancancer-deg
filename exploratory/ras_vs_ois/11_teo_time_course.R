# GSE115301 (Teo 2019) Smart-seq2 ER:RAS IMR90 time course: growing, day 2, day 4, senescent.
# Do RAS-autonomous genes rise before the arrest and arrest-coupled genes later?
# Uses the authors' published (downsampled, QC'd) cells as-is. House single-cell state filters
# (control-like, LMNB1, G1) are NOT applied: they would remove pre-arrest cells by design; this
# is a time course, not a senescent-state definition. Raw counts read READ-ONLY from
# senescence_stress; UCell via scutils. Run from the pancancer-deg repo root.
.libPaths(c("/Users/ro/Academic/10_Research/senescence_stress/modules/scrna_seq/GSE226225/final_analysis/.Rlib/main", .libPaths()))
suppressMessages({ library(Seurat); library(UCell); library(SummarizedExperiment) })
source("/Users/ro/code/academic/scutils/R/ucell_parallel.R")
f <- "/Users/ro/Academic/10_Research/senescence_stress/modules/scrna_seq/cs_dataset_survey/data/raw/GSE115301/GSE115301_Time_course_downsampled_HTSeq.txt.gz"
out <- "results/exploratory/ras_vs_ois/scrna_GSE115301"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
x <- read.delim(gzfile(f), row.names = 1, check.names = FALSE)
x <- x[grepl("^ENSG", rownames(x)), ]
ann <- as.data.frame(rowData(readRDS("data/exploratory/ras_vs_ois/SRP062871_gene_rse.rds")))
ann$id <- sub("\\..*", "", ann$gene_id)
sym <- ann$gene_name[match(sub("\\..*", "", rownames(x)), ann$id)]
keep <- !is.na(sym)
m <- rowsum(as.matrix(x[keep, ]), sym[keep])
tp <- sub("_.*", "", colnames(m))
cat("cells per time point:\n"); print(table(tp))
o <- CreateSeuratObject(counts = m, meta.data = data.frame(row.names = colnames(m), tp = tp))

cls <- read.csv("results/exploratory/ras_vs_ois/strict/ois_gene_classes.csv")
prolif <- readLines("results/exploratory/two_sided_2026-11-01/prolif_genes.txt")
sets <- list(
  MAPK_output = c("DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "FOSL1", "EGR1"),
  Proliferation = prolif,
  RASauto_UP = cls$gene[cls$OIS_direction == "UP" & cls$arrest_class == "RAS-autonomous"],
  intermediate_UP = cls$gene[cls$OIS_direction == "UP" & cls$arrest_class == "intermediate"],
  arrest_UP = cls$gene[cls$OIS_direction == "UP" & cls$arrest_class == "arrest-coupled"],
  NFKB_BRD4_dep_UP = cls$gene[cls$OIS_direction == "UP" & cls$nfkb_brd4 == "dependent"],
  RASauto_DOWN = cls$gene[cls$OIS_direction == "DOWN" & cls$arrest_class == "RAS-autonomous"],
  arrest_DOWN = cls$gene[cls$OIS_direction == "DOWN" & cls$arrest_class == "arrest-coupled"]
)
source(".claude/worktrees/two-sided-fisher/exploratory/ras_vs_ois/pathway_sets.R"); sets <- c(sets, pathway_sets())  # added 2026-11-07: pathway scores (TGF-beta, cGAS-STING, IFN, ...)
sets <- lapply(sets, function(g) intersect(g, rownames(o)))
sets$Proliferation <- head(sets$Proliferation, 1400)  # UCell maxRank 1500
o <- score_ucell_parallel(o, signatures = sets, ncores = 1, assay = "RNA", slot = "counts")
sc <- o[[]][, c("tp", grep("_UCell$", colnames(o[[]]), value = TRUE))]
colnames(sc) <- sub("_UCell$", "", colnames(sc))
sc$MKI67 <- as.numeric(GetAssayData(o, layer = "counts")["MKI67", ] > 0)
write.csv(sc, file.path(out, "time_course_cell_scores.csv"))
lev <- intersect(c("Growing", "day2", "day4", "Sen1", "Sen2", "Sen"), unique(sc$tp))
lev <- c(lev, setdiff(unique(sc$tp), lev))
med <- sapply(setdiff(colnames(sc), "tp"), function(k) tapply(sc[[k]], factor(sc$tp, levels = lev), median))
cat("\nMedian score per time point (MKI67 = fraction of cells detected):\n"); print(round(med, 3))
# fraction of the final (senescent) change already reached at each time point, per score
ref <- lev[1]; fin <- tail(lev[grepl("^Sen", lev)], 1)
if (length(fin)) {
  frac <- sweep(sweep(med, 2, med[ref, ], "-"), 2, med[fin, ] - med[ref, ], "/")
  cat("\nFraction of the", ref, "->", fin, "change reached:\n"); print(round(frac, 2))
}

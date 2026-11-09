# GSE245807 (Olan 2024, IMR90 ER:RAS): RAS-exposed cells that are still cycling vs RAS cells
# in G1, each against phase-matched growing cells. House filters kept (published QC, nuclear
# integrity, control-like kNN removal); the house LMNB1/G1 subset is NOT applied for this
# comparison (approved by Roberto 2026-11-02), but reproduced to check the rebuild.
# Reads senescence_stress objects and scutils READ-ONLY; writes only to pancancer-deg.
.libPaths(c("/Users/ro/Academic/10_Research/senescence_stress/modules/scrna_seq/GSE226225/final_analysis/.Rlib/main", .libPaths()))
suppressMessages({ library(Seurat); library(UCell) })
source("/Users/ro/code/academic/scutils/R/ucell_parallel.R")
H <- "/Users/ro/Academic/10_Research/senescence_stress/modules/scrna_seq/cs_dataset_survey/house_rerun/GSE245807"
out <- "results/exploratory/ras_vs_ois/scrna_GSE245807"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(1)

# House rebuild of a cell set (lines 328-341 of scRNA_seq_preprocess_final_oneplate_savefull_res_opt1.R)
rebuild <- function(obj_old) {
  counts <- GetAssayData(obj_old, assay = "RNA", layer = "counts")
  obj <- CreateSeuratObject(counts = counts, meta.data = obj_old@meta.data)
  obj$LMNB1 <- as.numeric(counts["LMNB1", , drop = TRUE])
  obj <- SCTransform(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  for (nb in c(24, 20, 16, 12, 8)) {
    res <- tryCatch(CellCycleScoring(obj, s.features = cc.genes.updated.2019$s.genes,
                                     g2m.features = cc.genes.updated.2019$g2m.genes, nbin = nb), error = function(e) NULL)
    if (!is.null(res)) { obj <- res; break }
  }
  stopifnot(!is.null(res))
  obj
}

raw <- readRDS(file.path(H, "raw_seurat_objs.rds"))
final_n <- c(RIS_1 = 773, RIS_2 = 652)
cells <- list()
for (b in names(raw)) {
  o <- rebuild(raw[[b]]$clean)
  n_house <- sum(o$LMNB1 == 0 & o$Phase == "G1")
  cat(sprintf("%s: %d cells after control-like removal; house LMNB1==0 & G1 subset = %d (house final %d)\n",
              b, ncol(o), n_house, final_n[[b]]))
  cells[[b]] <- o
}
ctrl_ids <- WhichCells(raw[[1]]$merged, expression = accession == "CTRL_1")
cells[["CTRL_1"]] <- rebuild(subset(raw[[1]]$merged, cells = ctrl_ids))
rm(raw); gc()

cls <- read.csv("results/exploratory/ras_vs_ois/strict/ois_gene_classes.csv")
sets <- list(
  MAPK_output = c("DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "FOSL1", "EGR1"),
  RASauto_UP = cls$gene[cls$OIS_direction == "UP" & cls$arrest_class == "RAS-autonomous"],
  intermediate_UP = cls$gene[cls$OIS_direction == "UP" & cls$arrest_class == "intermediate"],
  arrest_UP = cls$gene[cls$OIS_direction == "UP" & cls$arrest_class == "arrest-coupled"],
  RASauto_DOWN = cls$gene[cls$OIS_direction == "DOWN" & cls$arrest_class == "RAS-autonomous"],
  arrest_DOWN = cls$gene[cls$OIS_direction == "DOWN" & cls$arrest_class == "arrest-coupled"],
  NFKB_BRD4_dep_UP = cls$gene[cls$OIS_direction == "UP" & cls$nfkb_brd4 == "dependent"]
)
source(".claude/worktrees/two-sided-fisher/exploratory/ras_vs_ois/pathway_sets.R"); sets <- c(sets, pathway_sets())  # added 2026-11-07: pathway scores (TGF-beta, cGAS-STING, IFN, ...)
sc <- do.call(rbind, lapply(names(cells), function(b) {
  o <- cells[[b]]; DefaultAssay(o) <- "RNA"
  s <- lapply(sets, function(g) intersect(g, rownames(o)))
  o <- score_ucell_parallel(o, signatures = s, ncores = 1, assay = "RNA", slot = "counts")
  m <- o[[]][, grep("_UCell$", colnames(o[[]])), drop = FALSE]
  colnames(m) <- sub("_UCell$", "", colnames(m))
  data.frame(cell = colnames(o), plate = b, Phase = o$Phase, LMNB1 = o$LMNB1, m)
}))
sc$cond <- ifelse(sc$plate == "CTRL_1", "Growing", "RAS")
sc$cycling <- ifelse(sc$Phase == "G1", "G1", "cycling")
sc$group <- with(sc, ifelse(cond == "Growing", paste0("Growing_", cycling),
                     ifelse(cycling == "cycling", "RAS_cycling",
                     ifelse(LMNB1 == 0, "RAS_G1_LMNB1neg (house final)", "RAS_G1_LMNB1pos"))))
write.csv(sc, file.path(out, "cycling_vs_arrested_cell_scores.csv"), row.names = FALSE)
cat("\ncells per group:\n"); print(table(sc$group, sc$plate))

feat <- names(sets)
cmp <- list(c("RAS_cycling", "Growing_cycling"), c("RAS_G1_LMNB1neg (house final)", "Growing_G1"),
            c("RAS_G1_LMNB1pos", "Growing_G1"), c("RAS_cycling", "RAS_G1_LMNB1neg (house final)"))
res <- do.call(rbind, lapply(feat, function(f) do.call(rbind, lapply(cmp, function(cc) {
  a <- sc[[f]][sc$group == cc[1]]; b <- sc[[f]][sc$group == cc[2]]
  data.frame(score = f, comparison = paste(cc, collapse = " vs "), n1 = length(a), n2 = length(b),
             median1 = median(a), median2 = median(b), diff = median(a) - median(b),
             p = wilcox.test(a, b)$p.value)
}))))
res$padj_BH <- ave(res$p, ifelse(grepl("^HS_", res$score), "house", ifelse(grepl("^P_(GLYC|GAG|HEXOSAMINE)", res$score), "met", ifelse(grepl("^P_", res$score), "path", "core"))), FUN = function(p) p.adjust(p, "BH"))  # pathway scores (2026-11-07) are a separate BH family
write.csv(res, file.path(out, "cycling_vs_arrested_comparisons.csv"), row.names = FALSE)
# retention: RAS effect in cycling cells as a fraction of the RAS effect in house-final G1 cells
ret <- sapply(feat, function(f) {
  e_cyc <- res$diff[res$score == f & res$comparison == "RAS_cycling vs Growing_cycling"]
  e_g1 <- res$diff[res$score == f & res$comparison == "RAS_G1_LMNB1neg (house final) vs Growing_G1"]
  e_cyc / e_g1
})
cat("\nMedian UCell differences (RAS minus phase-matched growing) and retention in cycling RAS cells:\n")
print(data.frame(score = feat,
                 cycling = sapply(feat, function(f) res$diff[res$score == f & res$comparison == "RAS_cycling vs Growing_cycling"]),
                 G1_final = sapply(feat, function(f) res$diff[res$score == f & res$comparison == "RAS_G1_LMNB1neg (house final) vs Growing_G1"]),
                 retention = ret), digits = 3, row.names = FALSE)

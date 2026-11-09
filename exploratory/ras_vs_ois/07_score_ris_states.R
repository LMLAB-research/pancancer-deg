# Score the RAS-vs-OIS gene classes in house GSE245807 RIS final states (all house QC
# filters kept, G1 senescent cells only). State 3 has MAPK output back at growing level
# while arrested: RAS-autonomous genes should drop there, arrest-coupled genes should not.
# Reads senescence_stress objects and scutils READ-ONLY; writes only to pancancer-deg.
.libPaths(c("/Users/ro/Academic/10_Research/senescence_stress/modules/scrna_seq/GSE226225/final_analysis/.Rlib/main", .libPaths()))
suppressMessages({ library(Seurat); library(UCell) })
source("/Users/ro/code/academic/scutils/R/ucell_parallel.R")
CS <- "/Users/ro/Academic/10_Research/senescence_stress/modules/scrna_seq/cs_dataset_survey"
out <- "results/exploratory/ras_vs_ois/scrna_GSE245807"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

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
groups <- read.csv(file.path(CS, "notch_secondary/out/step2/cell_scores_GSE245807_RIS.csv.gz"))[, c("cell", "group", "plate")]

objs <- readRDS(file.path(CS, "house_rerun/GSE245807/final/RIS_seurat.rds"))
scores <- do.call(rbind, lapply(names(objs), function(p) {
  o <- objs[[p]]; DefaultAssay(o) <- "RNA"
  s <- lapply(sets, function(g) intersect(g, rownames(o)))
  sc <- score_ucell_parallel(o, signatures = s, ncores = 1, assay = "RNA", slot = "counts")
  m <- sc[[]][, grep("_UCell$", colnames(sc[[]])), drop = FALSE]
  data.frame(cell = colnames(o), plate = p, m, check.names = FALSE)
}))
# groups file uses "S_<barcode>" plus a plate column; join on plate + barcode
groups$cell <- sub("^S_", "", groups$cell)
scores <- merge(scores, groups[, c("cell", "plate", "group")], by = c("cell", "plate"))
write.csv(scores, file.path(out, "cell_class_scores.csv"), row.names = FALSE)

cat("cells per group:\n"); print(table(scores$group, scores$plate))
feat <- grep("_UCell$", colnames(scores), value = TRUE)
res <- do.call(rbind, lapply(feat, function(f) {
  do.call(rbind, lapply(c("state_1", "state_2"), function(ref) {
    a <- scores[[f]][scores$group == "state_3"]; b <- scores[[f]][scores$group == ref]
    data.frame(score = sub("_UCell$", "", f), comparison = paste0("state_3 vs ", ref),
               median_state3 = median(a), median_ref = median(b),
               log2_ratio = log2((median(a) + 1e-3) / (median(b) + 1e-3)),
               p = wilcox.test(a, b)$p.value)
  }))
}))
res$padj_BH <- ave(res$p, ifelse(grepl("^HS_", res$score), "house", ifelse(grepl("^P_(GLYC|GAG|HEXOSAMINE)", res$score), "met", ifelse(grepl("^P_", res$score), "path", "core"))), FUN = function(p) p.adjust(p, "BH"))  # pathway scores (2026-11-07) are a separate BH family
write.csv(res, file.path(out, "state3_vs_states12.csv"), row.names = FALSE)
print(res, digits = 3, row.names = FALSE)

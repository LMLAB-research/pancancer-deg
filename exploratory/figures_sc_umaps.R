# WHY: UMAPs for the single-cell tests in FINDINGS.md, drawn with the house plotting calls
#   (house plot_states_umaps.R: DimPlot clusters first, no labels, larger points; customSignaturePlot
#   for scores; save_p). Scores are the UCell class scores already saved by scripts 07, 09 and 13.
# USAGE: from the pancancer-deg repo root:
#   Rscript --vanilla .claude/worktrees/two-sided-fisher/exploratory/figures_sc_umaps.R
# OUT: exploratory/figures/umap_*.png (git-tracked). Reads senescence_stress objects READ-ONLY.
Sys.setenv(KMP_WARNINGS = "0")
SC <- "/Users/ro/Academic/10_Research/senescence_stress/modules/scrna_seq/GSE226225/final_analysis"
.libPaths(c(file.path(SC, ".Rlib/main"), .libPaths()))
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork); library(viridis) })
source(file.path(SC, "lib/scRNA_seq.R"))       # customSignaturePlot
source(file.path(SC, "lib/general_scripts.R"))  # save_p
FIG <- file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))), "figures")
dir.create(FIG, showWarnings = FALSE)
H <- "/Users/ro/Academic/10_Research/senescence_stress/modules/scrna_seq/cs_dataset_survey/house_rerun"
R <- "results/exploratory/ras_vs_ois"
score_cols <- c(MAPK_output = "MAPK output", RASauto_UP = "RAS-autonomous UP", arrest_UP = "Arrest-coupled UP",
                NFKB_BRD4_dep_UP = "NF-kB/BRD4-dependent UP", RASauto_DOWN = "RAS-autonomous DOWN",
                arrest_DOWN = "Arrest-coupled DOWN")
FD <- "results/exploratory/figure_data"; dir.create(FD, recursive = TRUE, showWarnings = FALSE)
save_meta <- function(o, tag, extra = NULL) {   # per-cell cluster, group and scores for the cluster tests
  m <- o[[]][, intersect(c("seurat_clusters", extra, names(score_cols)), colnames(o[[]])), drop = FALSE]
  write.csv(data.frame(cell = colnames(o), m), file.path(FD, paste0(tag, "_cells.csv")), row.names = FALSE)
}
umap_panel <- function(o, title, extra_group = NULL) {
  pt <- if (ncol(o) > 3000) 1 else if (ncol(o) > 1000) 1.6 else 3
  o$seurat_clusters <- factor(as.character(o$seurat_clusters),
                              levels = sort(unique(as.integer(as.character(o$seurat_clusters)))))
  p <- list(DimPlot(o, reduction = "umap", group.by = "seurat_clusters", pt.size = pt, label = FALSE) +
              ggtitle(sprintf("clusters (n = %d)", ncol(o))))
  if (!is.null(extra_group)) p <- c(p, list(DimPlot(o, reduction = "umap", group.by = extra_group, pt.size = pt, label = FALSE) + ggtitle(extra_group)))
  for (k in intersect(names(score_cols), colnames(o[[]])))
    p <- c(p, list(customSignaturePlot(o, feature = k, pt.size = pt, name = "scaled UCell", remove_axes = TRUE) +
                     ggtitle(score_cols[[k]]) + theme(text = element_text(size = 14), plot.title = element_text(size = 14))))
  wrap_plots(p, ncol = 4) + plot_annotation(title = title, theme = theme(plot.title = element_text(size = 18)))
}

# 1. GSE245807 house final RAS states, scores from script 07
l <- readRDS(file.path(H, "GSE245807/final/RIS_seurat.rds"))
o <- merge(l[[1]], l[-1], merge.dr = "umap")
sc <- read.csv(file.path(R, "scrna_GSE245807/cell_class_scores.csv"))
key_o <- paste(o$accession, sub("_[0-9]+$", "", colnames(o)))
for (k in names(score_cols)) if (paste0(k, "_UCell") %in% colnames(sc) || k %in% colnames(sc)) {
  v <- if (k %in% colnames(sc)) sc[[k]] else sc[[paste0(k, "_UCell")]]
  o[[k]] <- v[match(key_o, paste(sc$plate, sc$cell))]
}
cat("GSE245807 final: scores matched for", sum(!is.na(o$RASauto_UP)), "of", ncol(o), "cells\n")
save_meta(o, "GSE245807_final_states", "accession")
save_p(umap_panel(o, "GSE245807 RAS (house final states)", extra_group = "accession"),
       file_name = "umap_GSE245807_RAS_final_states", save_dir = paste0(FIG, "/"), p_width = 20, p_height = 10)

# 2. GSE173879 OIS arms (house one-plate run), scores from script 13
D <- "data/exploratory/house_GSE173879_OIS"
ws <- read.csv(file.path(R, "scrna_GSE173879/withdrawal_cell_scores.csv"), row.names = 1)
for (cond in c("OISON", "OISOFF")) {
  o <- readRDS(file.path(D, paste0(cond, "_seurat.rds")))[[1]]
  lab <- if (cond == "OISON") "OIS, 4-OHT kept (house final)" else "OIS, 4-OHT withdrawn (house final)"
  w <- ws[ws$group == lab, ]
  rn <- sub("\\.[0-9]+$", "", rownames(w))   # rbind made row names unique with .1 suffixes
  for (k in names(score_cols)) o[[k]] <- w[[k]][match(colnames(o), rn)]
  cat(cond, ": scores matched for", sum(!is.na(o$RASauto_UP)), "of", ncol(o), "cells\n")
  save_meta(o, paste0("GSE173879_", cond))
  save_p(umap_panel(o, paste("GSE173879", lab)), file_name = paste0("umap_GSE173879_", cond),
         save_dir = paste0(FIG, "/"), p_width = 20, p_height = 10)
}

# 3. GSE245807 RAS cells before the house LMNB1/G1 subset (script 09 groups), UMAP on the house rebuild PCA
cs <- read.csv(file.path(R, "scrna_GSE245807/cycling_vs_arrested_cell_scores.csv"))
raw <- readRDS(file.path(H, "GSE245807/raw_seurat_objs.rds"))
objs <- lapply(names(raw), function(b) { x <- raw[[b]]$clean; x$plate <- b; x })
o <- merge(objs[[1]], objs[-1]); rm(raw, objs); gc()
DefaultAssay(o) <- "RNA"; o[["SCT"]] <- NULL; o[["RNA"]] <- JoinLayers(o[["RNA"]]); o <- SCTransform(o, verbose = FALSE); o <- RunPCA(o, npcs = 50, verbose = FALSE)
set.seed(1); o <- RunUMAP(o, dims = 1:30, verbose = FALSE)
o <- FindNeighbors(o, dims = 1:30, verbose = FALSE); o <- FindClusters(o, resolution = 0.12, verbose = FALSE)
cs <- cs[cs$cond == "RAS", ]
key <- paste(o$plate, colnames(o)); keyc <- paste(cs$plate, cs$cell)
o$group <- cs$group[match(key, keyc)]
for (k in names(score_cols)) o[[k]] <- cs[[k]][match(key, keyc)]
cat("GSE245807 pre-filter RAS cells: groups matched for", sum(!is.na(o$group)), "of", ncol(o), "\n")
save_meta(o, "GSE245807_prefilter", "group")
save_p(umap_panel(o, "GSE245807 RAS cells before LMNB1/G1 subset", extra_group = "group"),
       file_name = "umap_GSE245807_RAS_cycling_vs_G1", save_dir = paste0(FIG, "/"), p_width = 20, p_height = 10)

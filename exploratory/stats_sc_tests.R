# WHY: statistics for the single-cell, time-course and dose comparisons and for the UMAP clusters,
#   using the house wilcox_dotplot_wrapper() (lib/scRNA_seq.R): Kruskal-Wallis per score (BH across
#   scores) plus all pairwise Wilcoxon tests (rstatix, BH across all tests in the call), violin + jitter,
#   significant brackets only. Also bootstrap 95% CIs for the two ratio read-outs in Fig. 3d/e.
# USAGE: from the pancancer-deg repo root:
#   Rscript --vanilla .claude/worktrees/two-sided-fisher/exploratory/stats_sc_tests.R
# OUT: exploratory/figures/stats_*.png; exploratory/summary_tables/06_stats_*_{kw,wilcox}.csv,
#   06_bootstrap_ratios.csv. Reads local per-cell score files (results/exploratory/...).
Sys.setenv(KMP_WARNINGS = "0")
SC <- "/Users/ro/Academic/10_Research/senescence_stress/modules/scrna_seq/GSE226225/final_analysis"
.libPaths(c(file.path(SC, ".Rlib/main"), .libPaths()))
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(rstatix); library(ggpubr); library(patchwork) })
source(file.path(SC, "lib/scRNA_seq.R"))       # wilcox_dotplot_wrapper
source(file.path(SC, "lib/general_scripts.R"))  # save_p
E <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))))
FIG <- paste0(file.path(E, "figures"), "/"); TAB <- file.path(E, "summary_tables")
R <- "results/exploratory"
lab <- c(MAPK_output = "MAPK output", RASauto_UP = "RAS-autonomous UP", arrest_UP = "Arrest-coupled UP",
         NFKB_BRD4_dep_UP = "NF-kB/BRD4-dependent UP", RASauto_DOWN = "RAS-autonomous DOWN",
         arrest_DOWN = "Arrest-coupled DOWN")

run <- function(df, group_col, levels, tag, title) {
  df <- df[df[[group_col]] %in% levels, c(group_col, names(lab))]
  df[[group_col]] <- factor(df[[group_col]], levels = levels)
  colnames(df) <- c("group", unname(lab))
  res <- wilcox_dotplot_wrapper(df, group_col = "group", score_cols = unname(lab), ylab = "UCell score")
  p <- res$p + ggtitle(title) + theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), strip.text = element_text(size = 9))
  save_p(p, file_name = paste0("stats_", tag), save_dir = FIG, p_width = 13, p_height = 9)
  w <- as.data.frame(res$wilcox); w <- w[, !vapply(w, is.list, logical(1))]
  write.csv(w, file.path(TAB, paste0("06_stats_", tag, "_wilcox.csv")), row.names = FALSE)
  write.csv(as.data.frame(res$kw), file.path(TAB, paste0("06_stats_", tag, "_kw.csv")), row.names = FALSE)
  cat(sprintf("%-26s n=%d  KW BH-significant scores: %d of %d; pairwise BH<0.05: %d of %d\n", tag, nrow(df),
              sum(res$kw$p.adj < 0.05), nrow(res$kw), sum(w$p.adj < 0.05), nrow(w)))
}

# Time course, dose, KRAS days
tc <- read.csv(file.path(R, "ras_vs_ois/scrna_GSE115301/time_course_cell_scores.csv"))
tc$grp <- ifelse(grepl("^Growing", tc$tp), "Growing", ifelse(grepl("^Sen", tc$tp), "Senescent", tc$tp))
run(tc, "grp", c("Growing", "day2", "day4", "Senescent"), "GSE115301_time_course", "ER:RAS IMR90 time course (GSE115301)")
ds <- read.csv(file.path(R, "ras_vs_ois/dose_and_kras/GSE222951_RAS_dose_cell_scores.csv"))
run(ds, "group", c("S", "M", "L", "XL"), "GSE222951_RAS_dose", "RAS dose, RPE1 (GSE222951)")
kr <- read.csv(file.path(R, "ras_vs_ois/dose_and_kras/GSE202126_HPNE_KRAS_timecourse_noUMIcap_cell_scores.csv"))
run(kr, "group", c("parental", "day1", "day2", "day3", "day5"), "GSE202126_KRAS_days", "KRAS time course, HPNE (GSE202126; no UMI cap)")

# RAS withdrawal and cycling cells
wd <- read.csv(file.path(R, "ras_vs_ois/scrna_GSE173879/withdrawal_cell_scores.csv"), row.names = 1)
wl <- c("Proliferating (CTRL)", "OIS, 4-OHT kept (house final)", "OIS, 4-OHT withdrawn (house final)")
run(wd, "group", wl, "GSE173879_withdrawal", "RAS withdrawn after OIS (GSE173879, house pipeline)")
cy <- read.csv(file.path(R, "ras_vs_ois/scrna_GSE245807/cycling_vs_arrested_cell_scores.csv"))
cl <- c("Growing_G1", "Growing_cycling", "RAS_cycling", "RAS_G1_LMNB1neg (house final)")
run(cy, "group", cl, "GSE245807_cycling", "GSE245807: RAS cycling vs G1, with phase-matched growing cells")

# UMAP clusters (cluster vs cluster, all pairs)
fd <- function(f) read.csv(file.path(R, "figure_data", f))
x <- fd("GSE245807_final_states_cells.csv"); run(x, "seurat_clusters", sort(unique(x$seurat_clusters)), "clusters_GSE245807_final_states", "GSE245807 house final RAS states: clusters")
x <- fd("GSE245807_prefilter_cells.csv"); run(x, "seurat_clusters", sort(unique(x$seurat_clusters)), "clusters_GSE245807_prefilter", "GSE245807 RAS cells before LMNB1/G1 subset: clusters")
for (cond in c("OISON", "OISOFF")) {
  x <- fd(paste0("GSE173879_", cond, "_cells.csv"))
  run(x, "seurat_clusters", sort(unique(x$seurat_clusters)), paste0("clusters_GSE173879_", cond), paste("GSE173879", cond, "house final: clusters"))
}

# Bootstrap 95% CIs for the ratio read-outs (resample cells within each group; seed per score)
boot_ratio <- function(num_a, num_b, den_a, den_b, B = 2000, seed = 1) {
  set.seed(seed)
  est <- (median(num_a) - median(num_b)) / (median(den_a) - median(den_b))
  bs <- replicate(B, (median(sample(num_a, replace = TRUE)) - median(sample(num_b, replace = TRUE))) /
                     (median(sample(den_a, replace = TRUE)) - median(sample(den_b, replace = TRUE))))
  c(estimate = est, lo = unname(quantile(bs, 0.025, na.rm = TRUE)), hi = unname(quantile(bs, 0.975, na.rm = TRUE)))
}
out <- list()
for (i in seq_along(lab)) {
  k <- names(lab)[i]
  g <- function(d, l) d[[k]][d$group == l]
  r1 <- boot_ratio(g(wd, wl[3]), g(wd, wl[1]), g(wd, wl[2]), g(wd, wl[1]), seed = 1000 + i)
  r2 <- boot_ratio(g(cy, "RAS_cycling"), g(cy, "Growing_cycling"), g(cy, "RAS_G1_LMNB1neg (house final)"), g(cy, "Growing_G1"), seed = 2000 + i)
  out[[length(out) + 1]] <- data.frame(readout = "GSE173879 fraction of OIS change remaining after withdrawal", score = lab[[k]], t(r1))
  out[[length(out) + 1]] <- data.frame(readout = "GSE245807 RAS effect in cycling / G1 cells", score = lab[[k]], t(r2))
}
bt <- do.call(rbind, out)
write.csv(bt, file.path(TAB, "06_bootstrap_ratios.csv"), row.names = FALSE)
print(bt, digits = 3, row.names = FALSE)

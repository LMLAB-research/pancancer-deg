# WHY: summary figures for FINDINGS.md, drawn from the git-tracked summary_tables/ (so they can be
#   regenerated without the local results). Palette: dataviz reference categorical slots 1-5, validated
#   (all checks pass; three slots < 3:1 contrast, so every series is direct-labelled).
# USAGE: Rscript --vanilla .claude/worktrees/two-sided-fisher/exploratory/figures_summary.R  (any directory)
# OUT: exploratory/figures/fig*.png
.libPaths(c("/Users/ro/Academic/10_Research/senescence_stress/modules/scrna_seq/GSE226225/final_analysis/.Rlib/main", .libPaths()))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
E <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))))
T <- file.path(E, "summary_tables"); FIG <- file.path(E, "figures"); dir.create(FIG, showWarnings = FALSE)
col <- c("RAS-autonomous" = "#2a78d6", "Arrest-coupled" = "#eb6834", "Intermediate" = "#1baf7a",
         "MAPK output" = "#eda100", "NF-kB/BRD4-dependent" = "#e87ba4")
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid <- "#e1e0d9"
th <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = grid, linewidth = 0.3),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA), text = element_text(colour = ink),
        axis.text = element_text(colour = ink2), strip.text = element_text(face = "bold", hjust = 0),
        legend.position = "none", plot.title.position = "plot")
save <- function(p, f, w, h) ggsave(file.path(FIG, f), p, width = w, height = h, dpi = 200, bg = "#fcfcfb")
rd <- function(f) read.csv(file.path(T, f), check.names = FALSE)
stars <- function(p) ifelse(p < 1e-4, "****", ifelse(p < 1e-3, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "n.s."))))

# Fig 1: overlap of each OIS gene class with TCGA tumour DEGs (same direction), per cancer
o <- rd("02_tcga_overlap_by_class_strict.csv")
o <- o[o$cls %in% c("RAS-autonomous", "intermediate", "arrest-coupled"), ]
o$class <- factor(c("RAS-autonomous" = "RAS-autonomous", intermediate = "Intermediate", "arrest-coupled" = "Arrest-coupled")[o$cls],
                  levels = c("RAS-autonomous", "Intermediate", "Arrest-coupled"))
o$lor <- pmax(log2(pmax(o$odds_ratio, 1e-3)), -4)
o$dir <- factor(ifelse(o$ois_dir == "UP", "OIS-up genes vs tumour-up genes", "OIS-down genes vs tumour-down genes"),
                levels = c("OIS-up genes vs tumour-up genes", "OIS-down genes vs tumour-down genes"))
o$sig <- o$p < 0.05
med <- aggregate(lor ~ class + dir, o, median)
# paired Wilcoxon signed-rank tests between classes across the 17 cancers (paired by cancer), BH over the 6 tests
pairs <- list(c("RAS-autonomous", "Arrest-coupled"), c("RAS-autonomous", "Intermediate"), c("Intermediate", "Arrest-coupled"))
pt <- do.call(rbind, lapply(levels(o$dir), function(d) do.call(rbind, lapply(seq_along(pairs), function(i) {
  a1 <- o[o$dir == d & o$class == pairs[[i]][1], ]; b1 <- o[o$dir == d & o$class == pairs[[i]][2], ]
  m <- merge(a1[, c("cancer", "lor")], b1[, c("cancer", "lor")], by = "cancer")
  data.frame(dir = d, group1 = pairs[[i]][1], group2 = pairs[[i]][2], n_cancers = nrow(m),
             median_diff = median(m$lor.x - m$lor.y), p = wilcox.test(m$lor.x, m$lor.y, paired = TRUE)$p.value, i = i)
}))))
pt$padj <- p.adjust(pt$p, "BH")
write.csv(pt[, setdiff(names(pt), "i")], file.path(T, "07_fig1_paired_class_tests.csv"), row.names = FALSE)
lv <- levels(o$class)
ytop <- tapply(o$lor, o$dir, max)
pt$x1 <- match(pt$group1, lv); pt$x2 <- match(pt$group2, lv)
pt$y <- ytop[as.character(pt$dir)] + 0.35 * pt$i
pt$dir <- factor(pt$dir, levels = levels(o$dir))
pt$label <- stars(pt$padj)
p1 <- ggplot(o, aes(class, lor, colour = class)) +
  geom_hline(yintercept = 0, colour = "#c3c2b7") +
  geom_jitter(aes(shape = sig), width = 0.15, size = 2.6, stroke = 0.8) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  geom_crossbar(data = med, aes(ymin = lor, ymax = lor), width = 0.5, colour = ink, linewidth = 0.4) +
  geom_segment(data = pt, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, colour = ink2, linewidth = 0.4) +
  geom_text(data = pt, aes(x = (x1 + x2) / 2, y = y + 0.08, label = label), inherit.aes = FALSE, size = 4.2, colour = ink2) +
  scale_colour_manual(values = col) + facet_wrap(~dir) +
  labs(x = NULL, y = "log2 odds ratio (17 TCGA cancers)",
       title = "Overlap of OIS gene classes with TCGA tumour-vs-normal DEGs (same direction)",
       subtitle = paste0("Each point is one cancer; filled = p < 0.05 (two-sided Fisher, uncorrected); bar = median.\n",
                         "Brackets: paired Wilcoxon signed-rank across the 17 cancers (BH over 6 tests): * < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001.")) + th
save(p1, "fig1_tcga_overlap_by_class.png", 11, 6)

# Fig 2: fraction of each class changed in the same direction, per oncogene arm
a <- rd("04_oncogene_arms_class_sharing.csv")
keep <- c("UP RAS-autonomous", "DOWN RAS-autonomous", "UP intermediate", "DOWN intermediate",
          "UP arrest-coupled", "DOWN arrest-coupled", "UP NF-kB/BRD4-dependent")
long <- do.call(rbind, lapply(keep, function(k) data.frame(arm = a$arm, class = k, frac = a[[paste0("frac_same_dir: ", k)]],
                                                          MKI67 = a$log2FC_MKI67)))
long$arm <- factor(paste0(long$arm, "  (MKI67 ", sprintf("%+.1f", long$MKI67), ")"),
                   levels = rev(unique(paste0(a$arm, "  (MKI67 ", sprintf("%+.1f", a$log2FC_MKI67), ")"))))
long$class <- factor(long$class, levels = keep)
p2 <- ggplot(long, aes(class, arm, fill = frac)) +
  geom_tile(colour = "#fcfcfb", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", frac), colour = frac > 0.55), size = 3.6) +
  scale_colour_manual(values = c(`TRUE` = "#ffffff", `FALSE` = ink), guide = "none") +
  scale_fill_gradient(low = "#cde2fb", high = "#104281", limits = c(0, 1), name = "fraction") +
  labs(x = NULL, y = NULL, title = "Which OIS gene classes each oncogene arm reproduces",
       subtitle = "Fraction of the class changed in the same direction (house DEG thresholds); MKI67 log2FC marks arrest") +
  th + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right", panel.grid = element_blank())
save(p2, "fig2_oncogene_arms.png", 11, 5.5)

# Fig 3: single-cell and time/dose tests (scores min-max scaled within each dataset and score)
sc_show <- c(MAPK_output = "MAPK output", RASauto_UP = "RAS-autonomous", arrest_DOWN = "Arrest-coupled",
             NFKB_BRD4_dep_UP = "NF-kB/BRD4-dependent")
scale01 <- function(x) (x - min(x)) / (max(x) - min(x))
stat_lab <- c(MAPK_output = "MAPK output", RASauto_UP = "RAS-autonomous UP", arrest_DOWN = "Arrest-coupled DOWN",
              NFKB_BRD4_dep_UP = "NF-kB/BRD4-dependent UP")
# marker: square = reference group; filled = BH < 0.05 vs reference (house pairwise Wilcoxon, stats_sc_tests.R); open = n.s.
add_sig <- function(df, tag) {
  w <- rd(paste0("06_stats_", tag, "_wilcox.csv")); w$score <- sub(" \\(kw.*$", "", w$var_new)
  ref <- levels(df$x)[1]
  df$mark <- vapply(seq_len(nrow(df)), function(i) {
    if (df$x[i] == ref) return("reference")
    sl <- stat_lab[[names(sc_show)[sc_show == df$score[i]]]]
    r <- w[w$score == sl & ((w$group1 == ref & w$group2 == df$x[i]) | (w$group2 == ref & w$group1 == df$x[i])), ]
    if (nrow(r) && r$p.adj[1] < 0.05) "BH < 0.05 vs reference" else "n.s. vs reference" }, "")
  df$mark <- factor(df$mark, levels = c("reference", "BH < 0.05 vs reference", "n.s. vs reference")); df
}
line_panel <- function(df, xlab, title, note = NULL) {
  df$score <- factor(df$score, levels = sc_show)
  ggplot(df, aes(x, y, colour = score, group = score)) +
    geom_line(linewidth = 1) + geom_point(aes(shape = mark), size = 3, stroke = 1, fill = "#fcfcfb") +
    scale_shape_manual(values = c("reference" = 15, "BH < 0.05 vs reference" = 16, "n.s. vs reference" = 21), name = NULL, drop = FALSE) +
    scale_colour_manual(values = col, name = NULL) +
    labs(x = xlab, y = "scaled median score (0-1)", title = title, subtitle = note) + th +
    theme(legend.position = "bottom")
}
# 3a GSE115301 time course
tc <- rd("03_GSE115301_time_course_medians.csv")
tc$x <- ifelse(grepl("^Growing", tc$group), "Growing", ifelse(grepl("^Sen", tc$group), "Senescent", tc$group))
tca <- aggregate(. ~ x, tc[, c("x", names(sc_show))], median)
tca$x <- factor(tca$x, levels = c("Growing", "day2", "day4", "Senescent"))
df <- do.call(rbind, lapply(names(sc_show), function(k) data.frame(x = tca$x, score = sc_show[[k]], y = scale01(tca[[k]]))))
df$x <- factor(df$x, levels = levels(tca$x)); df <- add_sig(df, "GSE115301_time_course")
p3a <- line_panel(df, NULL, "ER:RAS time course (GSE115301)", "Arrest-coupled = DOWN genes; falls as arrest sets in")
# 3b GSE222951 dose
ds <- rd("03_GSE222951_RAS_dose_summary.csv")
df <- do.call(rbind, lapply(names(sc_show), function(k) { r <- ds[ds$score == k, c("S", "M", "L", "XL")]
  data.frame(x = factor(c("S", "M", "L", "XL"), levels = c("S", "M", "L", "XL")), score = sc_show[[k]], y = scale01(unlist(r))) }))
df <- add_sig(df, "GSE222951_RAS_dose")
p3b <- line_panel(df, "RAS-reporter level", "RAS dose, RPE1 (GSE222951)")
# 3c GSE202126 KRAS days (no UMI cap run)
kr <- rd("03_GSE202126_HPNE_KRAS_timecourse_noUMIcap_summary.csv")
lv <- c("parental", "day1", "day2", "day3", "day5")
df <- do.call(rbind, lapply(names(sc_show), function(k) { r <- kr[kr$score == k, lv]
  data.frame(x = factor(lv, levels = lv), score = sc_show[[k]], y = scale01(unlist(r))) }))
df <- add_sig(df, "GSE202126_KRAS_days")
p3c <- line_panel(df, NULL, "KRAS time course, HPNE (GSE202126)", "Run without the 120k UMI cap (see DECISIONS.md)")
# 3d GSE173879 withdrawal: fraction of the OIS change remaining after 4-OHT withdrawal
w <- rd("03_GSE173879_ras_withdrawal_medians.csv")
g <- function(lab) unlist(w[w$group == lab, names(sc_show)])
on <- g("OIS, 4-OHT kept (house final)"); off <- g("OIS, 4-OHT withdrawn (house final)"); c0 <- g("Proliferating (CTRL)")
bt <- rd("06_bootstrap_ratios.csv")
ci <- function(readout) { b <- bt[bt$readout == readout, ]; b$score <- c("MAPK output" = "MAPK output", "RAS-autonomous UP" = "RAS-autonomous",
  "Arrest-coupled DOWN" = "Arrest-coupled", "NF-kB/BRD4-dependent UP" = "NF-kB/BRD4-dependent")[b$score]; b[!is.na(b$score), ] }
wd <- data.frame(score = factor(sc_show, levels = sc_show), frac = (off - c0) / (on - c0))
wd <- merge(wd, ci("GSE173879 fraction of OIS change remaining after withdrawal")[, c("score", "lo", "hi")], by = "score")
p3d <- ggplot(wd, aes(score, frac, fill = score)) + geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, colour = ink, linewidth = 0.4) +
  geom_hline(yintercept = c(0, 1), colour = "#c3c2b7", linetype = c(1, 2)) +
  geom_text(aes(label = sprintf("%.2f", frac), y = ifelse(frac >= 0, hi, lo), vjust = ifelse(frac >= 0, -0.5, 1.4)), size = 3.8, colour = ink) +
  scale_fill_manual(values = col) + scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(x = NULL, y = "fraction of OIS change remaining", title = "RAS withdrawn after OIS (GSE173879, house pipeline)",
       subtitle = "1 = unchanged by withdrawal; 0 = back to proliferating\nerror bars: 95% bootstrap CI (2,000 resamples of cells)") + th
# 3e GSE245807 cycling RAS cells: RAS effect retained vs phase-matched growing cells
cy <- rd("03_GSE245807_cycling_vs_arrested.csv")
ret <- sapply(names(sc_show), function(k) {
  a1 <- cy$diff[cy$score == k & cy$comparison == "RAS_cycling vs Growing_cycling"]
  b1 <- cy$diff[cy$score == k & cy$comparison == "RAS_G1_LMNB1neg (house final) vs Growing_G1"]; a1 / b1 })
cd <- data.frame(score = factor(sc_show, levels = sc_show), frac = ret)
cd <- merge(cd, ci("GSE245807 RAS effect in cycling / G1 cells")[, c("score", "lo", "hi")], by = "score")
p3e <- p3d; p3e$data <- cd; p3e <- p3e + labs(y = "RAS effect in cycling cells / in G1 cells", title = "RAS cells still cycling (GSE245807)",
                         subtitle = "Each vs phase-matched growing cells; house G1 filter not applied\nerror bars: 95% bootstrap CI (2,000 resamples of cells)")
p3 <- ((p3a | p3b | p3c) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.box = "vertical")) /
  (p3d | p3e | plot_spacer()) + plot_layout(heights = c(1.1, 1))
save(p3, "fig3_single_cell_time_dose_tests.png", 20, 11)

# Fig 4: TCGA activity models (stroma + ABSOLUTE purity), per cohort
am <- rd("05_tcga_activity_models_stroma_purity.csv")
tg <- c(RIS_RASauto_UP = "RAS-autonomous UP", RIS_NFKB_BRD4_UP = "NF-kB/BRD4-dependent UP",
        RIS_arrest_UP = "Arrest-coupled UP", RIS_RASauto_DOWN = "RAS-autonomous DOWN")
am <- am[am$target %in% names(tg), ]
lg <- do.call(rbind, lapply(c(ERK_output = "ERK output", MYC_targets = "MYC targets", Proliferation = "Proliferation"),
  function(nm) NULL))
lg <- do.call(rbind, lapply(names(c(ERK_output = 1, MYC_targets = 1, Proliferation = 1)), function(p)
  data.frame(cohort = am$cohort, target = factor(tg[am$target], levels = tg), predictor = gsub("_", " ", p),
             b = am[[paste0("b_", p)]], sig = am[[paste0("p_", p)]] < 0.05)))
mm <- aggregate(b ~ target + predictor, lg, median)
# across cohorts: one-sample Wilcoxon (coefficient != 0) per target x predictor; paired Wilcoxon ERK vs MYC per target; BH within each family
one <- do.call(rbind, lapply(split(lg, list(lg$target, lg$predictor), drop = TRUE), function(x)
  data.frame(target = x$target[1], predictor = x$predictor[1], n_cohorts = nrow(x), median_b = median(x$b), p = wilcox.test(x$b, mu = 0)$p.value)))
one$padj <- p.adjust(one$p, "BH")
pr <- do.call(rbind, lapply(split(lg, lg$target), function(x) { e <- x[x$predictor == "ERK output", ]; m <- x[x$predictor == "MYC targets", ]
  m <- m[match(e$cohort, m$cohort), ]
  data.frame(target = x$target[1], n_cohorts = nrow(e), median_ERK_minus_MYC = median(e$b - m$b), p = wilcox.test(e$b, m$b, paired = TRUE)$p.value) }))
pr$padj <- p.adjust(pr$p, "BH")
write.csv(one, file.path(T, "07_fig4_coefficient_vs_zero_tests.csv"), row.names = FALSE)
write.csv(pr, file.path(T, "07_fig4_ERK_vs_MYC_paired_tests.csv"), row.names = FALSE)
ymax <- max(lg$b); ymin <- min(lg$b)
one$y <- ymin - 0.1; one$lab <- stars(one$padj)
pr$y <- ymax + 0.12; pr$lab <- stars(pr$padj)
p4 <- ggplot(lg, aes(predictor, b)) + geom_hline(yintercept = 0, colour = "#c3c2b7") +
  geom_jitter(aes(shape = sig), width = 0.15, size = 2.2, colour = "#2a78d6") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  geom_crossbar(data = mm, aes(ymin = b, ymax = b), width = 0.5, colour = ink, linewidth = 0.4) +
  geom_text(data = one, aes(predictor, y, label = lab), size = 4.2, colour = ink2) +
  geom_segment(data = pr, aes(x = 1, xend = 2, y = y - 0.04, yend = y - 0.04), colour = ink2, linewidth = 0.4) +
  geom_text(data = pr, aes(x = 1.5, y = y, label = lab), size = 4.2, colour = ink2) +
  facet_wrap(~target, nrow = 1) +
  labs(x = NULL, y = "standardised coefficient", title = "TCGA: what each OIS gene class tracks across tumours (20 cohorts)",
       subtitle = paste0("score ~ ERK output + MYC targets + proliferation + stroma + ABSOLUTE purity, fitted within each cohort. Point = one cohort (filled: p < 0.05 in that cohort); bar = median.\n",
                          "Stars under a column: is that coefficient different from 0 across the 20 cohorts? Bracket on top: is the ERK coefficient different from the MYC coefficient?\n",
                          "Wilcoxon signed-rank across cohorts, BH; * < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001")) + th
save(p4, "fig4_tcga_activity_models.png", 15, 6)

# Fig 5: TCGA driver groups (stroma + purity): RAS-mutant or MYC-amplified vs neither
dm <- rd("05_tcga_driver_models_stroma_purity.csv")
dm <- dm[dm$score %in% c("ERK_output", "MYC_targets", "RIS_RASauto_UP", "RIS_NFKB_BRD4_UP", "RIS_arrest_UP", "RIS_RASauto_DOWN"), ]
dm$score <- factor(c(ERK_output = "ERK output", MYC_targets = "MYC targets", RIS_RASauto_UP = "RAS-autonomous UP",
                     RIS_NFKB_BRD4_UP = "NF-kB/BRD4 UP", RIS_arrest_UP = "Arrest-coupled UP", RIS_RASauto_DOWN = "RAS-autonomous DOWN")[dm$score],
                   levels = c("ERK output", "MYC targets", "RAS-autonomous UP", "NF-kB/BRD4 UP", "Arrest-coupled UP", "RAS-autonomous DOWN"))
dm$driver <- factor(ifelse(dm$driver == "RAS", "RAS-hotspot", "MYC-amplified"), levels = c("RAS-hotspot", "MYC-amplified"))
dm$lab <- ifelse(dm$padj < 0.05, stars(dm$padj), "")
lim <- max(abs(dm$effect_sd))
p5 <- ggplot(dm, aes(driver, cohort, fill = effect_sd)) +
  geom_tile(colour = "#fcfcfb", linewidth = 0.8) +
  geom_text(aes(label = lab, colour = abs(effect_sd) > 0.6 * lim), size = 4) +
  scale_colour_manual(values = c(`TRUE` = "#ffffff", `FALSE` = ink), guide = "none") +
  scale_fill_gradient2(low = "#184f95", mid = "#f0efec", high = "#e34948", midpoint = 0, limits = c(-lim, lim),
                       name = "effect\n(within-cohort SD)") +
  scale_y_discrete(limits = rev) +
  facet_wrap(~score, nrow = 1) +
  labs(x = NULL, y = NULL, title = "TCGA tumours with a RAS hotspot or MYC amplification, each compared with tumours that have neither",
       subtitle = paste0("Colour = difference in score (within-cohort SD), adjusted for stroma and ABSOLUTE purity; white/blank cell = not tested (< 10 tumours in a group).\n",
                         "Stars = BH-corrected significance over all cohort x score x driver tests: * < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001")) +
  th + theme(legend.position = "right", panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
save(p5, "fig5_tcga_driver_models.png", 15, 7.5)
cat("figures written to", FIG, "\n")

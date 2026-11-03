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
p1 <- ggplot(o, aes(class, lor, colour = class)) +
  geom_hline(yintercept = 0, colour = "#c3c2b7") +
  geom_jitter(aes(shape = sig), width = 0.15, size = 2.6, stroke = 0.8) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  geom_crossbar(data = med, aes(ymin = lor, ymax = lor), width = 0.5, colour = ink, linewidth = 0.4) +
  scale_colour_manual(values = col) + facet_wrap(~dir) +
  labs(x = NULL, y = "log2 odds ratio (17 TCGA cancers)",
       title = "Overlap of OIS gene classes with TCGA tumour-vs-normal DEGs (same direction)",
       subtitle = "Each point is one cancer; filled = p < 0.05 (two-sided Fisher); bar = median. Classes from script 04 (strict).") + th
save(p1, "fig1_tcga_overlap_by_class.png", 11, 5)

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
line_panel <- function(df, xlab, title, note = NULL) {
  df$score <- factor(df$score, levels = sc_show)
  lab <- df[df$x == levels(df$x)[nlevels(df$x)], ]
  ggplot(df, aes(x, y, colour = score, group = score)) +
    geom_line(linewidth = 1) + geom_point(size = 2.6) +
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
df$x <- factor(df$x, levels = levels(tca$x))
p3a <- line_panel(df, NULL, "ER:RAS time course (GSE115301)", "Arrest-coupled = DOWN genes; falls as arrest sets in")
# 3b GSE222951 dose
ds <- rd("03_GSE222951_RAS_dose_summary.csv")
df <- do.call(rbind, lapply(names(sc_show), function(k) { r <- ds[ds$score == k, c("S", "M", "L", "XL")]
  data.frame(x = factor(c("S", "M", "L", "XL"), levels = c("S", "M", "L", "XL")), score = sc_show[[k]], y = scale01(unlist(r))) }))
p3b <- line_panel(df, "RAS-reporter level", "RAS dose, RPE1 (GSE222951)")
# 3c GSE202126 KRAS days (no UMI cap run)
kr <- rd("03_GSE202126_HPNE_KRAS_timecourse_noUMIcap_summary.csv")
lv <- c("parental", "day1", "day2", "day3", "day5")
df <- do.call(rbind, lapply(names(sc_show), function(k) { r <- kr[kr$score == k, lv]
  data.frame(x = factor(lv, levels = lv), score = sc_show[[k]], y = scale01(unlist(r))) }))
p3c <- line_panel(df, NULL, "KRAS time course, HPNE (GSE202126)", "Run without the 120k UMI cap (see DECISIONS.md)")
# 3d GSE173879 withdrawal: fraction of the OIS change remaining after 4-OHT withdrawal
w <- rd("03_GSE173879_ras_withdrawal_medians.csv")
g <- function(lab) unlist(w[w$group == lab, names(sc_show)])
on <- g("OIS, 4-OHT kept (house final)"); off <- g("OIS, 4-OHT withdrawn (house final)"); c0 <- g("Proliferating (CTRL)")
wd <- data.frame(score = factor(sc_show, levels = sc_show), frac = (off - c0) / (on - c0))
p3d <- ggplot(wd, aes(score, frac, fill = score)) + geom_col(width = 0.6) +
  geom_hline(yintercept = c(0, 1), colour = "#c3c2b7", linetype = c(1, 2)) +
  geom_text(aes(label = sprintf("%.2f", frac), vjust = ifelse(frac >= 0, -0.4, 1.2)), size = 3.8, colour = ink) +
  scale_fill_manual(values = col) + scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(x = NULL, y = "fraction of OIS change remaining", title = "RAS withdrawn after OIS (GSE173879, house pipeline)",
       subtitle = "1 = unchanged by withdrawal; 0 = back to proliferating") + th
# 3e GSE245807 cycling RAS cells: RAS effect retained vs phase-matched growing cells
cy <- rd("03_GSE245807_cycling_vs_arrested.csv")
ret <- sapply(names(sc_show), function(k) {
  a1 <- cy$diff[cy$score == k & cy$comparison == "RAS_cycling vs Growing_cycling"]
  b1 <- cy$diff[cy$score == k & cy$comparison == "RAS_G1_LMNB1neg (house final) vs Growing_G1"]; a1 / b1 })
cd <- data.frame(score = factor(sc_show, levels = sc_show), frac = ret)
p3e <- p3d; p3e$data <- cd; p3e <- p3e + labs(y = "RAS effect in cycling cells / in G1 cells", title = "RAS cells still cycling (GSE245807)",
                         subtitle = "Each vs phase-matched growing cells; house G1 filter not applied")
p3 <- (p3a | p3b | p3c) / (p3d | p3e | plot_spacer()) + plot_layout(heights = c(1, 1))
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
p4 <- ggplot(lg, aes(predictor, b)) + geom_hline(yintercept = 0, colour = "#c3c2b7") +
  geom_jitter(aes(shape = sig), width = 0.15, size = 2.2, colour = "#2a78d6") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  geom_crossbar(data = mm, aes(ymin = b, ymax = b), width = 0.5, colour = ink, linewidth = 0.4) +
  facet_wrap(~target, nrow = 1) +
  labs(x = NULL, y = "standardised coefficient", title = "TCGA: what each OIS gene class tracks across tumours (20 cohorts)",
       subtitle = "score ~ ERK output + MYC targets + proliferation + stroma + ABSOLUTE purity, within cohort; filled = p < 0.05; bar = median") + th
save(p4, "fig4_tcga_activity_models.png", 15, 5)

# Fig 5: TCGA driver groups (stroma + purity): RAS-mutant or MYC-amplified vs neither
dm <- rd("05_tcga_driver_models_stroma_purity.csv")
dm <- dm[dm$score %in% c("ERK_output", "MYC_targets", "RIS_RASauto_UP", "RIS_NFKB_BRD4_UP", "RIS_arrest_UP", "RIS_RASauto_DOWN"), ]
dm$score <- factor(c(ERK_output = "ERK output", MYC_targets = "MYC targets", RIS_RASauto_UP = "RAS-autonomous UP",
                     RIS_NFKB_BRD4_UP = "NF-kB/BRD4 UP", RIS_arrest_UP = "Arrest-coupled UP", RIS_RASauto_DOWN = "RAS-autonomous DOWN")[dm$score],
                   levels = c("ERK output", "MYC targets", "RAS-autonomous UP", "NF-kB/BRD4 UP", "Arrest-coupled UP", "RAS-autonomous DOWN"))
dm$driver <- factor(paste0(dm$driver, ifelse(dm$driver == "RAS", "-hotspot tumours", "-amplified tumours"), " vs neither"))
p5 <- ggplot(dm, aes(effect_sd, cohort, colour = padj < 0.05)) + geom_vline(xintercept = 0, colour = "#c3c2b7") +
  geom_point(size = 2.6) + scale_colour_manual(values = c(`TRUE` = "#2a78d6", `FALSE` = "#898781")) +
  facet_grid(driver ~ score, scales = "free_y", space = "free_y") +
  labs(x = "effect (within-cohort SD)", y = NULL, title = "TCGA driver groups (stroma + ABSOLUTE purity adjusted)",
       subtitle = "blue = BH < 0.05; cohorts with >= 10 tumours per group") + th
save(p5, "fig5_tcga_driver_models.png", 17, 7)
cat("figures written to", FIG, "\n")

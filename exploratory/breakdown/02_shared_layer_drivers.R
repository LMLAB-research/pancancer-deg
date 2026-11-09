# The OIS genes that RS / SIPS share: are they the tumour-like part, and what drives them outside RAS? (exploratory, 2026-11-09)
# Gene sets (per direction), from ras_vs_ois/strict/ois_gene_classes.csv and arrest_degs_final_RERUN.csv:
#   RAS-autonomous, shared  = RAS-autonomous OIS DEGs also changed the same way in RS or SIPS (meta)
#   RAS-autonomous, OIS-only = RAS-autonomous OIS DEGs changed in neither RS nor SIPS
#   (and the same split of all OIS senescence-specific genes, i.e. not changed in CICQ / SSCQ)
#   reference: arrest-coupled OIS DEGs
# PART A: each set vs each cancer's tumour DEGs (house overlap_function, two-sided; background = OIS DEGs of that
#   direction tested in the cancer), BH across all tests.
# PART B: what reproduces each set outside RAS? For each perturbation contrast (vs its own proliferating control), house
#   overlap_function of each set with the contrast's same-direction DEGs; background = RAS-autonomous (or
#   senescence-specific) OIS DEGs of that direction tested in the contrast, so shared and OIS-only are compared directly.
#   Contrasts, grouped by what they switch on:
#     DNA damage, p53 wild-type: IR / doxorubicin / etoposide (GSE162175, IMR90), bleomycin (NHDF, NHLF), doxorubicin BJ 8 d,
#       IR GSE259295, hydroxyurea RPE1 p53-wt (GSE185512), IR fibroblasts 4 d / 20 d (ERP021140)
#     DNA damage, p53 lost: hydroxyurea RPE1 p53-/- (GSE185512), Adriamycin and H2O2 in p53-mutant LFS fibroblasts (GSE60340)
#     p53 activation without damage: nutlin, RG7388, HDM201
#     arrest without damage: palbociclib (IMR90, NHDF, NHLF), abemaciclib, quiescence (IMR90, TIG3), rapamycin
#     cytokines (house modules, up genes only): NF-kB (TNF + IL-1b 3 h), TNF 24 h, IL-1b 14 d
#     other oncogenes: NOTCH1-ICD (IMR90), AKT (BJ, IMR90 array), MYC and cyclin E (RPE1 p53-wt)
#   Dependence inside OIS (retention = knockdown LFC / control LFC, as 04_classify_and_project.py): shp65, shBrd4, JQ1,
#   shIL1R (day 10), shp53, shp16/p21; shared vs OIS-only compared by Wilcoxon on retention.
#   Not independent of the meta-analysis: IMR90 etoposide (SRP062872 DDIS) is in the meta SIPS group (shown, flagged).
# Run from the pancancer-deg repo root. OUT: results/exploratory/breakdown/shared_layer/
SA <- "/Users/ro/Academic/10_Research/senescence_stress/modules/systems_analysis_arrest"
BULK <- file.path(SA, "bulk_RNA_seq"); PD <- getwd()
OUT <- file.path(PD, "results/exploratory/breakdown/shared_layer"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(BULK, ".Rlib"), .libPaths())); Sys.setenv(BULK_RNASEQ_DIR = SA)
source(file.path(BULK, "R/config.R")); source(file.path(BULK, "R/functions.R"))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
LFC_T <- log2(1.5)
tab <- function(gene, lfc, sig) { d <- data.frame(gene = gene, lfc = lfc, dir = ifelse(lfc > 0, "up", "down"), sig = sig %in% "y")
  d[!is.na(d$gene) & d$gene != "" & !duplicated(d$gene), ] }

# ---- gene sets ----
cls <- read.csv(file.path(PD, "results/exploratory/ras_vs_ois/strict/ois_gene_classes.csv"))
ad <- read.csv(file.path(RERUN_DIR, "arrest_degs_final_RERUN.csv"))
md <- function(g) { x <- ad[ad$group_1 == g, ]; tab(x$gene, x$log2FoldChange, x$sig) }
RS <- md("Replicative_CS"); SI <- md("Stress_induced_CS"); CI <- md("Contact_inhibited_CQ"); SS <- md("Serum_starved_CQ"); OI <- md("Oncogene_induced_CS")
sd <- function(t, d) t$gene[t$sig & t$dir == d]
SETS <- list(); BGK <- list()
for (d in c("up", "down")) { D <- toupper(d)
  rssi <- union(sd(RS, d), sd(SI, d)); cq <- union(sd(CI, d), sd(SS, d))
  ra <- cls$gene[cls$OIS_direction == D & cls$arrest_class == "RAS-autonomous"]
  SETS[[paste("RAS-autonomous | shared with RS/SIPS", d)]] <- intersect(ra, rssi)
  SETS[[paste("RAS-autonomous | OIS-only", d)]] <- setdiff(ra, rssi)
  ss <- setdiff(sd(OI, d), cq)
  SETS[[paste("OIS senescence-specific | shared with RS/SIPS", d)]] <- intersect(ss, rssi)
  SETS[[paste("OIS senescence-specific | OIS-only", d)]] <- setdiff(ss, rssi)
  SETS[[paste("arrest-coupled (reference)", d)]] <- cls$gene[cls$OIS_direction == D & cls$arrest_class == "arrest-coupled"]
  BGK[[paste("RAS-autonomous", d)]] <- ra; BGK[[paste("OIS senescence-specific", d)]] <- ss
  BGK[[paste("arrest-coupled (reference)", d)]] <- sd(OI, d) }
SI_ <- data.frame(set = names(SETS), n = lengths(SETS)); print(SI_, row.names = FALSE)
fam <- function(s) sub(" \\| .*| (up|down)$", "", s); dirof <- function(s) sub(".* ", "", s)
G <- bind_rows(lapply(names(SETS), function(s) data.frame(set = s, family = fam(s), dir = dirof(s), gene = SETS[[s]])))
save_csv(G, file_name = "shared_layer_sets", path = OUT)

# ---- PART A: tumours ----
tf <- list.files(file.path(PD, "results/tables/TCGA_deg"), pattern = "_DGE_Raw_Results_Tumor_vs_Normal_Adjusted.csv$", full.names = TRUE)
TU <- lapply(setNames(tf, sub("^TCGA-", "", sub("_DGE.*", "", basename(tf)))), function(f) { d <- read.csv(f)
  d <- d[d$gene_type == "protein_coding" & !is.na(d$gene_name), ]; d <- d[!duplicated(d$gene_name), ]
  tab(d$gene_name, d$log2FoldChange, ifelse(!is.na(d$padj) & d$padj < .05 & abs(d$log2FoldChange) > LFC_T, "y", "n")) })
A <- bind_rows(lapply(c("up", "down"), function(d) { bg0 <- sd(OI, d)
  bind_rows(lapply(names(TU), function(ca) { t <- TU[[ca]]; bg <- intersect(bg0, t$gene)
    s <- G[G$dir == d & G$gene %in% bg, ]; tt <- t[t$sig & t$gene %in% bg, ]
    o <- suppressWarnings(overlap_function(df_1 = s, gene_col_1 = "gene", group_col_1 = "set", df_2 = tt, gene_col_2 = "gene",
                                           group_col_2 = "dir", background = bg))
    data.frame(cancer = ca, set = o$set, sen_dir = d, tumour_dir = o$dir, actual = o$actual, expected = o$expected,
               odds = o$odds, pval = o$pval, n_set = o$n1, background_n = length(bg)) })) }))
A$adj_BH <- p.adjust(A$pval, "BH"); A$relation <- ifelse(A$sen_dir == A$tumour_dir, "same", "opposite")
save_csv(A, file_name = "shared_layer_tcga_overlaps", path = OUT)
AS <- A |> group_by(set, relation) |> summarise(enriched = sum(adj_BH < .05 & odds > 1), depleted = sum(adj_BH < .05 & odds < 1),
                                                 med_log2OR = round(median(log2(pmax(odds, 1e-3))), 2), .groups = "drop")
# paired comparison across cancers: shared vs OIS-only log2 OR (same direction), Wilcoxon signed-rank
pc <- bind_rows(lapply(c("RAS-autonomous", "OIS senescence-specific"), function(f) bind_rows(lapply(c("up", "down"), function(d) {
  a <- A[A$set == paste(f, "| shared with RS/SIPS", d) & A$relation == "same", ]; b <- A[A$set == paste(f, "| OIS-only", d) & A$relation == "same", ]
  m <- merge(a, b, by = "cancer"); la <- log2(pmax(m$odds.x, 1e-3)); lb <- log2(pmax(m$odds.y, 1e-3))
  data.frame(family = f, dir = d, median_log2OR_shared = median(la), median_log2OR_OISonly = median(lb),
             n_shared_higher = sum(la > lb), n = nrow(m), p_wilcoxon = wilcox.test(la, lb, paired = TRUE)$p.value) }))))
pc$padj_BH <- p.adjust(pc$p_wilcoxon, "BH")
save_csv(AS, file_name = "shared_layer_tcga_summary", path = OUT); save_csv(pc, file_name = "shared_layer_tcga_paired", path = OUT)
options(width = 250)
cat("\n== PART A: n cancers (of 17) enriched/depleted (BH < 0.05), median log2 OR ==\n")
print(as.data.frame(AS |> mutate(cell = sprintf("%d/%d (%.2f)", enriched, depleted, med_log2OR)) |> dplyr::select(set, relation, cell) |>
        pivot_wider(names_from = relation, values_from = cell)), row.names = FALSE, right = FALSE)
cat("\nShared vs OIS-only, same-direction log2 OR paired across cancers:\n"); print(pc, row.names = FALSE, digits = 3)

# ---- PART B: perturbation contrasts ----
C <- file.path(PD, "results/exploratory/ras_vs_ois/contrast_degs")
cd <- function(n) { x <- read.csv(file.path(C, paste0(n, ".csv"))); tab(x$gene, x$log2FoldChange, x$sig) }
ng <- read.csv(file.path(RERUN_DIR, "nongenotoxic", "degs_nongenotoxic.csv"))
arm <- function(a) { x <- ng[ng$arm == a, ]; tab(x$gene, x$log2FoldChange, x$sig) }
lf <- read.csv(file.path(RERUN_DIR, "robustness", "damage_not_p53", "gse60340_degs.csv"))
lfs <- function(k) { x <- lf[lf$contrast == k, ]; tab(x$gene, x$log2FC, x$sig) }
irt <- read.csv("/Users/ro/Academic/10_Research/senescence_stress/work/results/auc/temporal/Fibroblast/degs.csv")
ir <- function(tp) { x <- irt[irt$group_1 == tp, ]; tab(x$gene, x$log2FoldChange, x$sig) }
MD <- file.path(BULK, "external_inputs", "stress_modules", "v3")
am <- read.csv(file.path(MD, "stress_modules.csv")); am <- am[am$version == "full", ]; atg <- read.csv(file.path(MD, "stress_module_tested_genes.csv"))
cm <- read.csv(file.path(MD, "chronic_modules.csv")); cm <- cm[cm$version == "full", ]; ctg <- read.csv(file.path(MD, "chronic_module_tested_genes.csv"))
modt <- function(mods, tg, m) { g <- unique(c(tg$gene[tg$module == m], mods$gene[mods$module == m]))
  data.frame(gene = g, lfc = ifelse(g %in% mods$gene[mods$module == m], 1, 0), dir = "up", sig = g %in% mods$gene[mods$module == m]) }
K <- list(
  "DNA damage, p53 wt" = list("IR (IMR90)" = arm("IR (IMR90)"), "Doxo (IMR90)" = arm("Doxo (IMR90)"), "Eto (IMR90)" = arm("Eto (IMR90)"),
    "Bleo (NHDF)" = arm("Bleo (NHDF)"), "Bleo (NHLF)" = arm("Bleo (NHLF)"), "Doxo 8d (BJ)" = arm("Doxo 8d (BJ)"),
    "IR (IMR90 GSE259295)" = arm("IR (IMR90 GSE259295)"), "HU RPE1 p53-wt" = cd("RPE1_p53wt_HU120h"),
    "IR fibroblast 4 d" = ir("4_days"), "IR fibroblast 20 d" = ir("20_days"), "Eto IMR90 SRP062872 (in meta)" = cd("IMR90_DDIS")),
  "DNA damage, p53 lost" = list("HU RPE1 p53-/-" = cd("RPE1_p53ko_HU120h"), "Adriamycin LFS p53-mut" = lfs("Adriamycin vs immortal"),
    "H2O2 LFS p53-mut" = lfs("H2O2 vs immortal")),
  "p53 on, no damage" = list("Nutlin (PAAF)" = arm("Nutlin 24h (PAAF)"), "RG7388 (IMR90)" = arm("RG7388 (IMR90)"), "HDM201 (IMR90)" = arm("HDM201 (IMR90)")),
  "arrest, no damage" = list("Palbo (IMR90)" = arm("Palbo (IMR90)"), "Palbo (NHDF)" = arm("Palbo (NHDF)"), "Palbo (NHLF)" = arm("Palbo (NHLF)"),
    "Abema 8d" = arm("Abema 8d"), "Quiescence (IMR90)" = arm("Quiescence (IMR90)"), "Quiescence TIG3" = cd("TIG3_quiescent_serum"),
    "Rapamycin (NHDF)" = arm("Rapamycin (NHDF)")),
  "cytokines (modules)" = list("NF-kB: TNF+IL-1b 3 h" = modt(am, atg, "TNF_3h+IL1B_3h"), "TNF 24 h" = modt(am, atg, "TNF_24h"),
    "IL-1b 14 d" = modt(cm, ctg, "IL1B_14d")),
  "other oncogenes" = list("NOTCH1-ICD (IMR90)" = cd("IMR90_N1ICD_NIS"), "AKT (BJ)" = cd("BJhT_AKT_d6"), "AKT (IMR90 array)" = cd("IMR90_array_AKT"),
    "MYC (RPE1 p53-wt)" = cd("RPE1_p53wt_Myc_dox120h"), "Cyclin E (RPE1 p53-wt)" = cd("RPE1_p53wt_CyclinE_dox120h")))
B <- bind_rows(lapply(names(K), function(grp) bind_rows(lapply(names(K[[grp]]), function(k) { x <- K[[grp]][[k]]
  bind_rows(lapply(names(BGK), function(bk) { d <- dirof(bk); f <- sub(" (up|down)$", "", bk)
    if (grp == "cytokines (modules)" && d == "down") return(NULL)
    bg <- intersect(BGK[[bk]], x$gene); s <- G[G$family == f & G$dir == d & G$gene %in% bg, ]
    xx <- x[x$sig & x$dir == d & x$gene %in% bg, ]; if (!nrow(xx) || !nrow(s)) return(NULL)
    o <- suppressWarnings(overlap_function(df_1 = s, gene_col_1 = "gene", group_col_1 = "set", df_2 = xx, gene_col_2 = "gene",
                                           group_col_2 = "dir", background = bg))
    data.frame(group = grp, contrast = k, family = f, dir = d, set = o$set, actual = o$actual, n_set = o$n1,
               frac = o$actual / o$n1, odds = o$odds, pval = o$pval, background_n = length(bg)) })) })) ))
B$adj_BH <- p.adjust(B$pval, "BH")
save_csv(B, file_name = "shared_layer_perturbations", path = OUT)
cat("\n== PART B: fraction of each set changed the same way in each contrast (* = BH < 0.05 vs the family background) ==\n")
pb <- B |> filter(family != "arrest-coupled (reference)") |>
  mutate(lab = sub("^.*\\| ", "", set), cell = sprintf("%.2f%s", frac, ifelse(adj_BH < .05, ifelse(odds > 1, "+", "-"), ""))) |>
  dplyr::select(group, contrast, family, dir, lab, cell) |> pivot_wider(names_from = c(family, dir, lab), values_from = cell)
print(as.data.frame(pb), row.names = FALSE, right = FALSE)

# ---- PART B2: dependence inside OIS (retention ratios) ----
ret <- function(num, den, d) { a <- cd(num); b <- cd(den); m <- merge(a, b, by = "gene"); m <- m[abs(m$lfc.y) >= LFC_T, ]
  setNames(m$lfc.x / m$lfc.y, m$gene) }
KD <- list("p65 (NF-kB) knockdown" = c("IMR90_RAS_shp65", "IMR90_RAS_shRen"), "BRD4 knockdown" = c("IMR90_RAS_shBrd4", "IMR90_RAS_shRen"),
           "JQ1 (BET inhibitor)" = c("IMR90_RAS_JQ1", "IMR90_RAS_DMSO"), "IL-1R knockdown (d10)" = c("IMR90hT_RAS_d10_shIL1R", "IMR90hT_RAS_d10_shScr"),
           "p53 knockdown" = c("IMR90_RAS_shp53", "IMR90_RAS_shRen"), "p16/p21 knockdown" = c("IMR90_RAS_shp16p21", "IMR90_RAS_shRen"))
R2 <- bind_rows(lapply(names(KD), function(k) { r <- ret(KD[[k]][1], KD[[k]][2])
  bind_rows(lapply(c("RAS-autonomous", "OIS senescence-specific"), function(f) bind_rows(lapply(c("up", "down"), function(d) {
    a <- r[intersect(names(r), SETS[[paste(f, "| shared with RS/SIPS", d)]])]; b <- r[intersect(names(r), SETS[[paste(f, "| OIS-only", d)]])]
    data.frame(knockdown = k, family = f, dir = d, n_shared = length(a), n_OISonly = length(b),
               median_retention_shared = median(a), median_retention_OISonly = median(b),
               frac_dependent_shared = mean(a <= 0.5), frac_dependent_OISonly = mean(b <= 0.5),
               p_wilcoxon = suppressWarnings(wilcox.test(a, b)$p.value)) })))) }))
R2$padj_BH <- p.adjust(R2$p_wilcoxon, "BH")
save_csv(R2, file_name = "shared_layer_ois_dependence", path = OUT)
cat("\n== PART B2: retention of the OIS change under each perturbation (1 = unchanged, <= 0.5 = dependent) ==\n")
print(R2 |> mutate(across(where(is.numeric), ~ signif(.x, 3))), row.names = FALSE)

# ---- figure: part B, fraction reproduced, shared vs OIS-only ----
F <- B |> filter(family != "arrest-coupled (reference)") |> mutate(lab = sub("^.*\\| ", "", set), contrast = factor(contrast, levels = rev(unique(contrast))))
p <- ggplot(F, aes(lab, contrast, fill = frac)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f%s", frac, ifelse(adj_BH < .05, ifelse(odds > 1, "+", "-"), ""))), size = 2.4) +
  scale_fill_gradient(low = "#f7f7f7", high = "#2c7bb6", name = "fraction changed\nthe same way") +
  facet_grid(group ~ family + dir, scales = "free_y", space = "free_y") + labs(x = NULL, y = NULL,
  caption = "+/- : BH < 0.05 over/under-represented relative to the family background (RAS-autonomous or OIS senescence-specific genes of that direction).") +
  theme_classic(base_size = 8) + theme(axis.text.x = element_text(angle = 30, hjust = 1), strip.text.y = element_text(angle = 0))
save_p(p, file_name = "fig_shared_layer_perturbations", save_dir = OUT, p_width = 11, p_height = 11)
cat("done\n")

# Breakdown of the CQ, RS and SIPS signatures into components, and which components tumours share (exploratory, 2026-11-09).
# The OIS split (RAS-autonomous vs arrest-coupled, ras_vs_ois/) has no equivalent for the other arrest types, which have
# no bypass / escape arms. Three splits that the local senescence_stress data support (all read-only):
#  (1) ARREST-SHARED vs SENESCENCE-SPECIFIC
#      meta:      RS / SIPS / OIS DEGs also changed the same way in CICQ or SSCQ ("arrest-shared") or in neither
#                 ("senescence-specific"); CICQ / SSCQ DEGs also changed the same way in RS, SIPS or OIS ("shared with CS")
#                 or in none ("CQ-specific"). Source: arrest_degs_final_RERUN.csv.
#      within-study (GSE162175, IMR90, one control): IR, doxorubicin, etoposide, palbociclib and RAS arms split by
#                 whether the Quiescence (0.5% FBS) arm changes the gene the same way. Source: degs_nongenotoxic.csv.
#  (2) p53 PROGRAMME vs REST
#      p53 programme = DEGs changed the same way in >= 2 of the 3 MDM2-inhibitor arms (nutlin PAAF, RG7388 and HDM201
#      IMR90: p53 activation without damage). Each meta group and each genotoxic arm split into p53-programme and rest.
#      Doxorubicin p53 dependence: DEGs in >= 2 of the 4 p53-wild-type doxorubicin arms, split by whether the p53-mutant
#      LFS Adriamycin contrast (GSE60340) changes them the same way (p53-independent) or not (p53-dependent).
#  (3) TIMING: IR time course (ERP021140, per cell type; 4 d vs 20 d) and doxorubicin in BJ (1 d vs 8 d):
#      early-persistent (both), early-only, late-only.
# Each component (direction-specific) vs each cancer's tumour-up and tumour-down DEGs: house overlap_function (two-sided
# Fisher), background = genes tested in the component's source(s) AND in that cancer. BH across all tests of a version.
# Versions: "full" and "noprolif" (Hallmark proliferation genes, the 671-gene list of DECISIONS B5, removed from
# components and backgrounds). Each component is also tested for house Hallmark + lysosomal sets (house_sets()),
# background = genes tested in its source(s); BH within component.
# Run from the pancancer-deg repo root. OUT: results/exploratory/breakdown/
source("exploratory/paths.R"); PD <- getwd()
OUT <- file.path(PD, "results/exploratory/breakdown"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source_house_bulk()   # house R/functions.R (copy in exploratory/inputs/house_code)
source("exploratory/ras_vs_ois/pathway_sets.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
LFC_T <- log2(1.5)
ARR <- c(Contact_inhibited_CQ = "CICQ", Serum_starved_CQ = "SSCQ", Replicative_CS = "RS", Stress_induced_CS = "SIPS", Oncogene_induced_CS = "OIS")

# ---- senescence inputs as gene x set tables: sig (TRUE/FALSE) and dir (up/down), tested = all rows ----
tab <- function(gene, lfc, sig) { d <- data.frame(gene = gene, dir = ifelse(lfc > 0, "up", "down"), sig = sig %in% "y")
  d[!is.na(d$gene) & d$gene != "" & !duplicated(d$gene), ] }
ad <- read.csv(SEN$arrest_rerun)
META <- lapply(setNames(names(ARR), ARR), function(g) { x <- ad[ad$group_1 == g, ]; tab(x$gene, x$log2FoldChange, x$sig) })
ng <- read.csv(SEN$nongenotoxic)
ARM <- lapply(split(ng, ng$arm), function(x) tab(x$gene, x$log2FoldChange, x$sig))
lf <- read.csv(SEN$gse60340)
LFS <- lapply(split(lf, lf$contrast), function(x) tab(x$gene, x$log2FC, x$sig))
IRT <- list()
for (ct in c("Fibroblast", "Keratinocyte", "Melanocyte")) {   # per-cell-type files (combined file is corrupt on disk)
  x <- read.csv(ir_time(ct))
  for (tp in c("4_days", "20_days")) { y <- x[x$group_1 == tp, ]; IRT[[paste(ct, tp)]] <- tab(y$gene, y$log2FoldChange, y$sig) } }
prolif <- readLines(PROLIF)

sig_dir <- function(t, d) t$gene[t$sig & t$dir == d]
COMP <- list()   # each: list(family, source, component, dir, genes, tested)
add <- function(family, source, component, d, genes, tested)
  COMP[[length(COMP) + 1]] <<- list(family = family, source = source, component = component, dir = d,
                                    genes = unique(genes), tested = unique(tested))
for (d in c("up", "down")) {
  # (1) meta arrest-shared vs senescence-specific
  cq <- union(sig_dir(META$CICQ, d), sig_dir(META$SSCQ, d)); cs <- Reduce(union, lapply(META[c("RS", "SIPS", "OIS")], sig_dir, d = d))
  tm <- Reduce(intersect, lapply(META, `[[`, "gene"))
  for (g in c("RS", "SIPS", "OIS")) { x <- sig_dir(META[[g]], d)
    add("1 arrest-shared vs specific (meta)", g, "arrest-shared (also in CQ)", d, intersect(x, cq), tm)
    add("1 arrest-shared vs specific (meta)", g, "senescence-specific (not in CQ)", d, setdiff(x, cq), tm) }
  for (g in c("CICQ", "SSCQ")) { x <- sig_dir(META[[g]], d)
    add("1 arrest-shared vs specific (meta)", g, "shared with CS", d, intersect(x, cs), tm)
    add("1 arrest-shared vs specific (meta)", g, "CQ-specific", d, setdiff(x, cs), tm) }
  # (1b) within GSE162175 vs its quiescence arm
  q <- ARM[["Quiescence (IMR90)"]]; qd <- sig_dir(q, d)
  for (a in c("IR (IMR90)", "Doxo (IMR90)", "Eto (IMR90)", "Palbo (IMR90)", "RAS (IMR90)")) { x <- sig_dir(ARM[[a]], d)
    tt <- intersect(ARM[[a]]$gene, q$gene)
    add("1b arrest-shared vs specific (GSE162175)", a, "arrest-shared (also in quiescence)", d, intersect(x, qd), tt)
    add("1b arrest-shared vs specific (GSE162175)", a, "senescence-specific (not in quiescence)", d, setdiff(x, qd), tt) }
  # (2) p53 programme from MDM2i arms
  mdm <- c("Nutlin 24h (PAAF)", "RG7388 (IMR90)", "HDM201 (IMR90)")
  p53 <- names(which(table(unlist(lapply(ARM[mdm], sig_dir, d = d))) >= 2))
  tmdm <- Reduce(intersect, lapply(ARM[mdm], `[[`, "gene"))
  for (g in c("RS", "SIPS", "OIS", "CICQ", "SSCQ")) { x <- sig_dir(META[[g]], d); tt <- intersect(META[[g]]$gene, tmdm)
    add("2 p53 programme (MDM2i)", g, "p53 programme", d, intersect(x, p53), tt)
    add("2 p53 programme (MDM2i)", g, "not p53 programme", d, setdiff(x, p53), tt) }
  for (a in c("IR (IMR90)", "Doxo (IMR90)", "Eto (IMR90)", "Bleo (NHDF)", "Doxo 8d (BJ)")) { x <- sig_dir(ARM[[a]], d); tt <- intersect(ARM[[a]]$gene, tmdm)
    add("2 p53 programme (MDM2i)", a, "p53 programme", d, intersect(x, p53), tt)
    add("2 p53 programme (MDM2i)", a, "not p53 programme", d, setdiff(x, p53), tt) }
  # (2b) doxorubicin: p53-wild-type arms vs p53-mutant LFS Adriamycin
  dx <- c("Doxo (IMR90)", "Doxo 1d (BJ)", "Doxo 8d (BJ)", "Doxo 24h (PAAF)")
  wt <- names(which(table(unlist(lapply(ARM[dx], sig_dir, d = d))) >= 2)); mt <- sig_dir(LFS[["Adriamycin vs immortal"]], d)
  tt <- intersect(Reduce(union, lapply(ARM[dx], `[[`, "gene")), LFS[["Adriamycin vs immortal"]]$gene)
  add("2b doxorubicin p53 dependence (LFS)", "Doxorubicin", "p53-dependent (lost in p53-mutant)", d, setdiff(intersect(wt, tt), mt), tt)
  add("2b doxorubicin p53 dependence (LFS)", "Doxorubicin", "p53-independent (kept in p53-mutant)", d, intersect(wt, mt), tt)
  # (3) timing
  for (ct in c("Fibroblast", "Keratinocyte", "Melanocyte")) { e <- IRT[[paste(ct, "4_days")]]; l <- IRT[[paste(ct, "20_days")]]
    tt <- intersect(e$gene, l$gene); ed <- sig_dir(e, d); ld <- sig_dir(l, d)
    add("3 timing", paste("IR", ct), "early-persistent (4 d and 20 d)", d, intersect(ed, ld), tt)
    add("3 timing", paste("IR", ct), "early-only (4 d)", d, setdiff(ed, ld), tt)
    add("3 timing", paste("IR", ct), "late-only (20 d)", d, setdiff(ld, ed), tt) }
  e <- ARM[["Doxo 1d (BJ)"]]; l <- ARM[["Doxo 8d (BJ)"]]; tt <- intersect(e$gene, l$gene); ed <- sig_dir(e, d); ld <- sig_dir(l, d)
  add("3 timing", "Doxo BJ", "early-persistent (1 d and 8 d)", d, intersect(ed, ld), tt)
  add("3 timing", "Doxo BJ", "early-only (1 d)", d, setdiff(ed, ld), tt)
  add("3 timing", "Doxo BJ", "late-only (8 d)", d, setdiff(ld, ed), tt)
}
CI <- bind_rows(lapply(seq_along(COMP), function(i) with(COMP[[i]], data.frame(id = i, family, source, component, dir,
  n_genes = length(genes), n_tested = length(tested), frac_prolif = round(mean(genes %in% prolif), 3)))))
save_csv(CI, file_name = "components", path = OUT)
save_csv(bind_rows(lapply(seq_along(COMP), function(i) with(COMP[[i]], if (length(genes)) data.frame(id = i, family, source, component, dir, gene = genes)))),
         file_name = "component_genes", path = OUT)
cat("components:", nrow(CI), "\n"); print(CI[, -1], row.names = FALSE)

# ---- tumour DEGs ----
tf <- tcga_deg_files()
TU <- lapply(setNames(tf, sub("^TCGA-", "", sub("_DGE.*", "", basename(tf)))), function(f) { d <- read.csv(f)
  d <- d[d$gene_type == "protein_coding" & !is.na(d$gene_name), ]; d <- d[!duplicated(d$gene_name), ]
  data.frame(gene = d$gene_name, dir = ifelse(d$log2FoldChange > 0, "up", "down"),
             sig = !is.na(d$padj) & d$padj < .05 & abs(d$log2FoldChange) > LFC_T) })

# ---- component x cancer overlaps (house overlap_function), full and noprolif ----
run_overlaps <- function(version) {
  drop <- if (version == "noprolif") prolif else character(0)
  bind_rows(lapply(seq_along(COMP), function(i) { cp <- COMP[[i]]; g <- setdiff(cp$genes, drop)
    if (length(g) < 5) return(NULL)
    bind_rows(lapply(names(TU), function(ca) { t <- TU[[ca]]; bg <- setdiff(intersect(cp$tested, t$gene), drop)
      tt <- t[t$sig & t$gene %in% bg, ]; if (length(intersect(g, bg)) < 5) return(NULL)
      o <- suppressWarnings(overlap_function(df_1 = data.frame(gene = intersect(g, bg), comp = "x"), gene_col_1 = "gene", group_col_1 = "comp",
                                             df_2 = tt, gene_col_2 = "gene", group_col_2 = "dir", background = bg))
      data.frame(version = version, id = i, family = cp$family, source = cp$source, component = cp$component, sen_dir = cp$dir,
                 cancer = ca, tumour_dir = o$dir, actual = o$actual, expected = o$expected, odds = o$odds, pval = o$pval,
                 n_comp = o$n1, n_tumour = o$n2, background_n = length(bg)) })) })) }
R <- bind_rows(run_overlaps("full"), run_overlaps("noprolif"))
R$adj_BH <- ave(R$pval, R$version, FUN = function(p) p.adjust(p, "BH"))
R$relation <- ifelse(R$sen_dir == R$tumour_dir, "same", "opposite")
save_csv(R, file_name = "component_tcga_overlaps", path = OUT)
SM <- R |> group_by(version, id, family, source, component, sen_dir, relation) |>
  summarise(n_comp = max(n_comp), enriched = sum(adj_BH < .05 & odds > 1), depleted = sum(adj_BH < .05 & odds < 1),
            med_log2OR = round(median(log2(pmax(odds, 1e-3))), 2), .groups = "drop")
save_csv(SM, file_name = "component_tcga_summary", path = OUT)
W <- SM |> mutate(cell = sprintf("%d/%d (%.2f)", enriched, depleted, med_log2OR)) |>
  dplyr::select(version, family, source, component, sen_dir, n_comp, relation, cell) |>
  pivot_wider(names_from = c(version, relation), values_from = cell)
options(width = 250); cat("\nn cancers (of 17) enriched/depleted, BH < 0.05 (median log2 OR)\n"); print(as.data.frame(W), row.names = FALSE, right = FALSE)

# ---- house Hallmark + lysosomal content of each component ----
hs <- house_sets(); pw <- bind_rows(lapply(names(hs), function(k) data.frame(pathway = k, gene = hs[[k]])))
E <- bind_rows(lapply(seq_along(COMP), function(i) { cp <- COMP[[i]]; if (length(cp$genes) < 5) return(NULL)
  o <- suppressWarnings(overlap_function(df_1 = data.frame(gene = cp$genes, comp = "x"), gene_col_1 = "gene", group_col_1 = "comp",
                                         df_2 = pw[pw$gene %in% cp$tested, ], gene_col_2 = "gene", group_col_2 = "pathway", background = cp$tested))
  data.frame(id = i, family = cp$family, source = cp$source, component = cp$component, dir = cp$dir, pathway = o$pathway,
             actual = o$actual, expected = o$expected, odds = o$odds, pval = o$pval, adj_BH = p.adjust(o$pval, "BH")) }))
save_csv(E, file_name = "component_house_pathways", path = OUT)
cat("\nTop house pathways per component (BH < 0.05, over-represented, top 4 by p):\n")
print(as.data.frame(E |> filter(adj_BH < .05, odds > 1) |> group_by(id, family, source, component, dir) |> arrange(pval) |>
        summarise(top = paste(sub("^HS_", "", head(pathway, 4)), collapse = ", "), .groups = "drop") |> dplyr::select(-id)), right = FALSE, row.names = FALSE)

# ---- figure: same / opposite direction, full version ----
F <- SM |> filter(version == "full") |>
  mutate(row = paste0(source, " | ", component, " | ", sen_dir, " (", n_comp, ")"), lab = paste0(enriched, "/", depleted))
F$row <- factor(F$row, levels = rev(unique(F$row[order(F$family, F$source, F$sen_dir, F$component)])))
p <- ggplot(F, aes(relation, row, fill = pmax(pmin(med_log2OR, 2), -2))) + geom_tile(colour = "white") + geom_text(aes(label = lab), size = 2.4) +
  scale_fill_gradient2(low = "#c0392b", mid = "#f7f7f7", high = "#2c7bb6", name = "median log2 OR\n(17 cancers)") +
  facet_grid(family ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Direction in tumours relative to the component", y = NULL,
       caption = "Label: n cancers enriched / depleted (two-sided Fisher, BH < 0.05). Row: source | component | direction (n genes).") +
  theme_classic(base_size = 8) + theme(strip.text.y = element_text(angle = 0, size = 7))
save_p(p, file_name = "fig_component_tcga_overlaps", save_dir = OUT, p_width = 8.5, p_height = 16)
cat("done\n")

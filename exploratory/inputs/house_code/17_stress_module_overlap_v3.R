# ---- V3 COPY (2026-10-16) -- this file is a copy of 17_stress_module_overlap.R (left unchanged, as are its outputs). ----
# V3 CHANGES (Roberto, 2026-10-16): Rotenone, BPTES and Chloramphenicol (GSE241261, 6 h) re-included as acute stress
#   modules (derive_modules_house_v3.R / export_house_for_bulk_v3.R). Only the inputs and outputs move:
#   inputs  external_inputs/stress_modules/v3/ (canonical v2 files in external_inputs/stress_modules/ untouched);
#   outputs rerun_outputs/stress_modules_v3/. Figures also show the three re-included modules. Tests, corrections and functions are unchanged.
# USAGE: Rscript meta_analysis/17_stress_module_overlap_v3.R   (from bulk_RNA_seq/; 17_v3 before 18_v3 and 19_v3)
#
# 17_stress_module_overlap.R
#
# Over-representation of (1) the 49 MSigDB Hallmark sets + the lysosomal list (50 sets, meta_analysis/03) and (2) the
# fibroblast ACUTE stress modules among the arrest DEGs (meta-analysis, 5 conditions) and the ERP021140 irradiation
# time-course DEGs (3 cell types x 3 time points). Rewritten 2026-10-15 (previous version: *.PRE_HOUSEFN) so that every
# test is the house R/functions.R::overlap_function (GeneOverlap 2x2 table, two-sided fisher.test, odds; raw pval and Holm
# 'adj' over the call) and every figure uses summarise_overlaps / summarise_overlaps_facet / save_p / save_csv.
#
# MODULES. Derived with the house DEG path (find_degs_new / DESeq2, padj < 0.05 and |log2FC| > log2(1.5), Ensembl 100;
#   limma with the same cut-offs for arrays) in scrna_seq/cs_dataset_survey/stress_modules/scripts/derive_modules_house.R,
#   exported by export_house_for_bulk.R to external_inputs/stress_modules/ (registered in env/external_inputs.tsv).
#   Since 2026-10-15 (v2; previous inputs/outputs *.PRE_V2) the inputs come from derive_modules_house_v2.R /
#   export_house_for_bulk_v2.R: "DNA damage" = GSE328392 lung fibroblasts 15 Gy 8 h (was AG1522 2 Gy 12 h); "ISR" (DNP)
#   passed the artefact test and is used without its artefact-like genes. This script is unchanged apart from this note.
#   "full" = all induced genes; "unique" = genes induced by no other included module. TNF 24 h is a readout only.
# TESTS.
#   Hallmark family: called exactly as meta_analysis/03 (arrest) and temporal_analysis/12 (time course): one
#     overlap_function call per DEG analysis, background = all genes tested in that DEG analysis. Raw pval therefore equals
#     03 / 12; 'adj' (Holm over the call) equals their adj.
#   Acute-module family: one overlap_function call per module and version (the background is module-specific): background =
#     genes tested in the DEG analysis AND tested in the module's source contrast(s).
#   Correction: adj_BH = BH on overlap_function's raw pval within each family per DEG analysis (owner's decision; modules
#     also per version, readout separately). Holm kept twice: adj_holm_call = overlap_function's 'adj' (Holm over its own
#     call: the whole Hallmark family, or one module x all DEG groups) and adj_holm_family = Holm over the same rows as BH.
#   Expression-matched permutation (no house equivalent; robustness check only): B = 20,000 random gene sets of the
#     module's size from the same background, each module gene replaced by one from the same baseMean decile of that DEG
#     contrast (genes not tested in that contrast count as baseMean 0); p_perm = (1 + k)/(B + 1), k = draws with overlap
#     >= observed; raw p only, no correction (house convention for simulation nulls), floor 5e-5.
# USAGE: Rscript meta_analysis/17_stress_module_overlap.R      (from bulk_RNA_seq/; PLOT_ONLY=1 redraws figures from the CSVs)
# OUT:   rerun_outputs/stress_modules_v3/{hallmark_overlap_arrest.csv, hallmark_overlap_temporal.csv,
#        stress_module_overlap_arrest.csv, stress_module_overlap_temporal.csv, fig_overlap_arrest.png,
#        fig_overlap_temporal.png}

source("R/config.R")
source("R/functions.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork) })
OUT <- file.path(RERUN_DIR, "stress_modules_v3"); dir.create(OUT, showWarnings = FALSE)
B <- 20000L
PLOT_MODS <- c("UPR", "Heat shock", "NF-kB", "DNA damage", "ISR", "Rotenone (complex I)", "BPTES (glutaminase)", "Chloramphenicol (mito translation)")   # v3: + re-included
ARR <- c(Contact_inhibited_CQ = "CICQ", Serum_starved_CQ = "SSCQ", Replicative_CS = "RS", Stress_induced_CS = "SIPS", Oncogene_induced_CS = "OIS")

# ---- inputs ----
hall <- read.csv(file.path(RERUN_DIR, "stress_response_pathways_RERUN.csv"))          # meta_analysis/03: 49 Hallmark + lysosomal
SM <- read.csv(file.path(EXTERNAL_DIR, "stress_modules", "v3", "stress_modules.csv"))
TG <- read.csv(file.path(EXTERNAL_DIR, "stress_modules", "v3", "stress_module_tested_genes.csv")); TG <- split(TG$gene, TG$module)
ad <- read.csv(file.path(RERUN_DIR, "arrest_degs_final_RERUN.csv"))
td <- read.csv(file.path(RERUN_DIR, "ERP021140", "all_time_degs.csv"))

# ---- family 1: Hallmark + lysosomal, called as 03 (arrest) and 12 (time course) ----
hallmark_arrest <- function() {
  sig <- ad[ad$sig == "y", ]
  o <- suppressWarnings(overlap_function(df_1 = hall, df_2 = sig, gene_col_1 = "genes", gene_col_2 = "gene",
                                         group_col_2 = c("dir_accession"), group_col_1 = c("pathway"),
                                         carry_col_2 = c("group_1", "direction_1"), background = ad$gene))
  data.frame(analysis = "arrest", cell_type = "all", group = ARR[as.character(o$group_1)], direction = o$direction_1, o) }
hallmark_temporal <- function() {
  t <- td; t$cell_type_time_point <- paste0(t$cell_type, "_", t$group_1); sig <- t[t$sig == "y", ]
  sig$cell_type_time_point_dir <- paste0(sig$cell_type_time_point, "_", sig$direction_1)
  o <- suppressWarnings(overlap_function(df_1 = hall, df_2 = sig, gene_col_1 = "genes", gene_col_2 = "ensembl",
                                         group_col_2 = "cell_type_time_point_dir", group_col_1 = c("pathway"),
                                         carry_col_2 = c("cell_type", "group_1", "direction_1"), background = unique(t$ensembl)))
  data.frame(analysis = "temporal", cell_type = o$cell_type, group = o$group_1, direction = o$direction_1, o) }
add_corr <- function(d, fam) {   # BH and Holm within each family (rows sharing the value of 'fam')
  d$adj_holm_call <- d$adj; d$adj_BH <- NA; d$adj_holm_family <- NA
  for (f in unique(fam)) { k <- fam == f; d$adj_BH[k] <- p.adjust(d$pval[k], "BH"); d$adj_holm_family[k] <- p.adjust(d$pval[k], "holm") }
  d }

# ---- family 2 (and 3 in 18): stress modules, one overlap_function call per module x version ----
# D: gene, baseMean, sig, direction, group, cell_type (one row per gene x contrast). MODS: version, module, module_label, gene.
module_family <- function(D, MODS, TESTED, analysis, label_cols = c("module_label")) {
  D <- D[!is.na(D$gene) & D$gene != "", ]; D$contrast <- paste(D$cell_type, D$group, sep = "|")
  D$contrast_dir <- paste(D$contrast, D$direction, sep = "|"); Dsig <- D[D$sig == "y", ]
  bg_analysis <- unique(D$gene); out <- list()
  for (md in unique(MODS$module)) { U <- intersect(bg_analysis, TESTED[[md]])
    for (v in c("full", "unique")) { g <- intersect(MODS$gene[MODS$module == md & MODS$version == v], U); if (length(g) < 5) next
      M <- data.frame(module = md, gene = g)
      o <- suppressWarnings(overlap_function(df_1 = M, df_2 = Dsig, gene_col_1 = "gene", gene_col_2 = "gene", group_col_1 = "module",
                                             group_col_2 = "contrast_dir", carry_col_2 = c("cell_type", "group", "direction", "contrast"), background = U))
      # expression-matched permutation per DEG contrast (draws shared by up and down, as before)
      o$perm_mean <- NA; o$p_perm <- NA
      for (cn in unique(o$contrast)) { b <- D[D$contrast == cn & D$gene %in% U, ]; b <- b[!duplicated(b$gene), ]
        bm <- setNames(rep(0, length(U)), U); bm[b$gene] <- b$baseMean
        dec <- setNames(cut(rank(bm, ties.method = "first"), 10, labels = FALSE), U); byd <- split(names(dec), dec); gd <- split(g, dec[g]); set.seed(1)
        draws <- replicate(B, unlist(lapply(names(gd), function(d) sample(byd[[d]], length(gd[[d]])))), simplify = FALSE)
        for (i in which(o$contrast == cn)) { degs <- Dsig$gene[Dsig$contrast_dir == o$contrast_dir[i] & Dsig$gene %in% U]
          inn <- vapply(draws, function(x) sum(x %in% degs), 0L); o$perm_mean[i] <- round(mean(inn), 2); o$p_perm[i] <- (1 + sum(inn >= o$actual[i])) / (B + 1) } }
      lab <- MODS[MODS$module == md, ][1, , drop = FALSE]
      out[[length(out) + 1]] <- data.frame(analysis = analysis, cell_type = o$cell_type, group = o$group, direction = o$direction, module = md,
        lab[, label_cols, drop = FALSE], version = v, set_size = o$n1, n_deg = o$n2, universe = o$background_n, overlap = o$actual, expected = round(o$expected, 2),
        odds = o$odds, pval = o$pval, adj = o$adj, genes = o$int, perm_mean = o$perm_mean, p_perm = o$p_perm, row.names = NULL) }
    cat(analysis, md, "done\n"); flush.console() }
  do.call(rbind, out) }

ad_m <- data.frame(gene = ad$gene, baseMean = ad$baseMean, sig = ad$sig, direction = ad$direction_1, group = ARR[ad$group_1], cell_type = "all")
td_m <- data.frame(gene = td$ensembl, baseMean = td$baseMean, sig = td$sig, direction = td$direction_1, group = td$group_1, cell_type = td$cell_type)
SM$readout_only <- as.logical(SM$readout_only)

if (Sys.getenv("PLOT_ONLY") == "1") { rd <- function(f) read.csv(file.path(OUT, paste0(f, ".csv")))
  A <- list(H = rd("hallmark_overlap_arrest"), S = rd("stress_module_overlap_arrest")); T <- list(H = rd("hallmark_overlap_temporal"), S = rd("stress_module_overlap_temporal"))
} else {
  A <- list(H = hallmark_arrest()); T <- list(H = hallmark_temporal())
  A$H <- add_corr(A$H, rep(1, nrow(A$H))); T$H <- add_corr(T$H, rep(1, nrow(T$H)))
  for (nm in c("A", "T")) { L <- get(nm)
    S <- module_family(if (nm == "A") ad_m else td_m, SM, TG, if (nm == "A") "arrest" else "temporal", label_cols = c("module_label", "readout_only"))
    L$S <- add_corr(S, paste(S$version, S$readout_only)); assign(nm, L) }
  save_csv(A$H, file_name = "hallmark_overlap_arrest", path = OUT); save_csv(A$S, file_name = "stress_module_overlap_arrest", path = OUT)
  save_csv(T$H, file_name = "hallmark_overlap_temporal", path = OUT); save_csv(T$S, file_name = "stress_module_overlap_temporal", path = OUT)
}

# ---- figures: house tile style (summarise_overlaps); Hallmark panel = the pathways of Figure 3a ----
POI <- c("HALLMARK TNFA SIGNALING VIA NFKB", "HALLMARK P53 PATHWAY", "HALLMARK MYC TARGETS", "HALLMARK MTORC1 SIGNALING",
         "HALLMARK MITOTIC SPINDLE", "HALLMARK INTERFERON GAMMA RESPONSE", "HALLMARK INTERFERON ALPHA RESPONSE",
         "HALLMARK INFLAMMATORY RESPONSE", "HALLMARK IL6 JAK STAT3 SIGNALING", "HALLMARK HYPOXIA", "HALLMARK G2M CHECKPOINT",
         "HALLMARK E2F TARGETS", "HALLMARK APOPTOSIS", "HALLMARK DNA REPAIR", "Lysosomal Genes")
nice_pw <- function(x) { x <- sub("^HALLMARK ", "", x); x <- tools::toTitleCase(tolower(x)); x <- sub("Nfkb", "NF-kB", x); x <- sub("Tnfa", "TNFA", x); x <- sub("Myc", "MYC", x)
  x <- sub("Mtorc1", "MTORC1", x); x <- sub("P53", "p53", x); x <- sub("Il6 Jak Stat3", "IL6 JAK STAT3", x); x <- sub("Dna", "DNA", x)
  x <- sub("G2m", "G2M", x); x <- sub("E2f", "E2F", x); x }
ylev <- c("SSCQ D", "CICQ D", "OIS D", "SIPS D", "RS D", "SSCQ U", "CICQ U", "OIS U", "SIPS U", "RS U")
prep_arrest <- function(d, xcol) { d$acc <- factor(paste(d$group, ifelse(d$direction == "up", "U", "D")), levels = ylev); d$x <- d[[xcol]]; d }
h <- prep_arrest(A$H[A$H$pathway %in% POI, ], "pathway"); h$x <- nice_pw(h$x)
s <- prep_arrest(A$S[A$S$version == "full" & A$S$module_label %in% PLOT_MODS, ], "module_label"); s$x <- factor(s$x, levels = PLOT_MODS)
p1 <- summarise_overlaps(h, x = "x", y = "acc", pval_col = "adj_BH", xlab = "MSigDB Hallmark and lysosomal gene sets", ylab = "Arrest-DEGs")
p2 <- summarise_overlaps(s, x = "x", y = "acc", pval_col = "adj_BH", xlab = "Fibroblast stress modules", ylab = NULL) + theme(axis.text.y = element_blank(), legend.position = "none")
save_p(p1 + p2 + plot_layout(widths = c(15, 5), guides = "collect"), file_name = "fig_overlap_arrest", save_dir = OUT, p_width = 9, p_height = 4.5)

tlev <- c("20_days D", "10_days D", "4_days D", "20_days U", "10_days U", "4_days U")
prep_t <- function(d, xcol) { d$summary_y <- factor(paste(d$group, ifelse(d$direction == "up", "U", "D")), levels = tlev)
  levels(d$summary_y) <- sub("_days", " Days", levels(d$summary_y)); d$cell_type <- tools::toTitleCase(d$cell_type); d$x <- d[[xcol]]; d }
th <- prep_t(T$H[T$H$pathway %in% POI, ], "pathway"); th$x <- nice_pw(th$x)
ts <- prep_t(T$S[T$S$version == "full" & T$S$module_label %in% PLOT_MODS, ], "module_label"); ts$x <- factor(ts$x, levels = PLOT_MODS)
q1 <- summarise_overlaps_facet(db = th, x = "x", y = "summary_y", pval_col = "adj_BH", xlab = "MSigDB Hallmark and lysosomal gene sets",
                               ylab = "Temporal DEGs", add_line = TRUE, legend_title = NULL, label_cq = TRUE, facet_col = "cell_type", angle = 45)
q2 <- summarise_overlaps_facet(db = ts, x = "x", y = "summary_y", pval_col = "adj_BH", xlab = "Fibroblast stress modules",
                               ylab = "Temporal DEGs", add_line = TRUE, legend_title = NULL, label_cq = TRUE, facet_col = "cell_type", angle = 45) + ylab(NULL) + theme(axis.text.y = element_blank())
save_p(q1 + theme(legend.position = "none") + q2 + plot_layout(widths = c(15, 5)), file_name = "fig_overlap_temporal", save_dir = OUT, p_width = 10, p_height = 7)

# ---- console summary: acute modules, odds (* = adj_BH < 0.05), both versions and directions ----
for (L in list(A$S, T$S)) for (v in c("full", "unique")) for (dr in c("up", "down")) {
  x <- L[L$version == v & L$direction == dr & !L$readout_only, ]; x$cell <- sprintf("%.2f%s", x$odds, ifelse(x$adj_BH < 0.05, "*", ""))
  cat("\n==", L$analysis[1], v, dr, "==\n")
  print(reshape(x[, c("cell_type", "group", "module_label", "cell")], idvar = c("cell_type", "group"), timevar = "module_label", direction = "wide"), row.names = FALSE) }
cat("done\n")

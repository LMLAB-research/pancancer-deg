# Metabolic pathway over-representation in TCGA tumour DEGs, compared with senescence (exploratory, 2026-11-07).
#
# COPY of the house senescence_stress/modules/systems_analysis_arrest/bulk_RNA_seq/meta_analysis/19_metabolic_overlap_v3.R
# (read-only; not modified). Same pathways (SI Table 37, keep == TRUE: 14 classes, 346 rows), same test (house
# R/functions.R::overlap_function: GeneOverlap + two-sided Fisher), same backgrounds (PRIMARY "SI37_pathway_genes" =
# SI37 genes AND tested; SENSITIVITY "SI18_enzymes_sensitivity" = SI18 enzymes AND tested), same correction (add_corr
# from 17_stress_module_overlap_v3.R: BH within family x bg_type) and the house summarise_overlaps tile figure.
# WHAT CHANGED (inputs only):
#   (1) TCGA tumour-vs-normal DEGs (pancancer-deg pipeline, covariate-adjusted design; protein-coding; sig = padj < 0.05
#       and |log2FC| > log2(1.5), the pipeline thresholds, equal to the house 'sig'). One overlap_function call per cancer,
#       background = SI37 / SI18 genes tested in that cancer (the module-style per-source background of 17/18/19).
#       Family for BH = all cancers x directions x classes, per bg_type (as the arrest family in 19).
#   (2) Per-inducer senescence arms (rerun_outputs/nongenotoxic/degs_nongenotoxic.csv; 21 arms), one call per arm
#       (background = genes tested in that arm), family = all arms, per bg_type.
#   (3) The senescence arrest results are READ from the house output (stress_modules_v3/metabolic_SI37/
#       metabolic_overlap_arrest.csv), not recomputed.
#   (4) New comparison (no house equivalent; house overlap_function only): within metabolic genes, do senescence DEGs
#       and tumour DEGs change in the same or opposite direction? df_1 = senescence arrest DEGs restricted to SI37 genes
#       (grouped by group x direction), df_2 = tumour DEGs (cancer x direction), background = SI37 genes tested in both.
#       BH within family = all senescence-group x cancer x direction pairs.
# USAGE (from the pancancer-deg repo root):
#   Rscript --vanilla exploratory/metabolism/01_metabolic_overlap_tcga.R
# OUT: results/exploratory/metabolism/
source("exploratory/paths.R"); PD <- getwd()
NOPROLIF <- Sys.getenv("NOPROLIF") == "1"   # 2026-11-07: Hallmark proliferation genes removed from pathways and backgrounds
OUT <- file.path(PD, "results/exploratory/metabolism", if (NOPROLIF) "noprolif" else ""); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source_house_bulk()   # house R/functions.R (copy in exploratory/inputs/house_code)
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork); library(readxl) })
house_add_corr()   # add_corr() and ARR from the house 17_stress_module_overlap_v3.R (copy in inputs/house_code)

# ---- house inputs (as 19_v3) ----
met_bg_tab <- read_si("si18"); pathways_genes <- read_si("si37")   # SI Tables 18 / 37 (extracted from si_tables.xlsx; inputs/senescence_stress)
metabolism_background <- unique(met_bg_tab$gene)
metabolism_genes <- pathways_genes[pathways_genes$keep %in% TRUE, ]
stopifnot(length(metabolism_background) == 2894, nrow(metabolism_genes) == 346)
if (NOPROLIF) {   # same 671-gene list as the main-pipeline proliferation sensitivity (DECISIONS B5)
  prolif <- readLines(PROLIF)
  cat("NOPROLIF: removing", sum(unique(metabolism_genes$gene) %in% prolif), "pathway genes and",
      sum(metabolism_background %in% prolif), "background genes\n")
  metabolism_genes <- metabolism_genes[!metabolism_genes$gene %in% prolif, ]; metabolism_background <- setdiff(metabolism_background, prolif)
}
BGS <- function(tested) list(SI37_pathway_genes = intersect(unique(metabolism_genes$gene), tested),
                             SI18_enzymes_sensitivity = intersect(metabolism_background, tested))
ad <- read.csv(SEN$arrest_rerun)
if (!NOPROLIF) {
  A <- read.csv(SEN$metabolic_arrest)
} else {   # recomputed with the house call of 19_v3 metabolic_arrest(), on the reduced pathways / backgrounds
  sig_ad <- ad[ad$sig == "y", ]
  A <- do.call(rbind, lapply(names(BGS(ad$gene)), function(bt) { bg <- BGS(ad$gene)[[bt]]
    o <- suppressWarnings(overlap_function(df_1 = sig_ad, gene_col_1 = "gene", group_col_1 = "dir_accession",
                                           carry_col_1 = c("group_1", "direction_1"), df_2 = metabolism_genes, gene_col_2 = "gene",
                                           group_col_2 = "class", background = bg))
    data.frame(analysis = "arrest", bg_type = bt, background_n = length(bg), cell_type = "all",
               group = ARR[as.character(o$group_1)], direction = o$direction_1, o) }))
  A <- add_corr(A, A$bg_type); save_csv(A, file_name = "metabolic_overlap_arrest_noprolif", path = OUT)
}

# ---- (1) TCGA tumour DEGs ----
LFC_T <- log2(1.5)
tf <- tcga_deg_files()
tum <- do.call(rbind, lapply(tf, function(f) {
  d <- read.csv(f); d <- d[d$gene_type == "protein_coding" & !is.na(d$gene_name) & d$gene_name != "", ]
  d <- d[!duplicated(d$gene_name), ]
  d$cancer <- sub("^TCGA-", "", sub("_DGE.*", "", basename(f)))
  d$sig <- ifelse(!is.na(d$padj) & d$padj < 0.05 & abs(d$log2FoldChange) > LFC_T, "y", "n")
  d$direction_1 <- ifelse(d$log2FoldChange > 0, "up", "down")
  d$dir_accession <- paste(d$cancer, d$direction_1, sep = "_")
  d[, c("gene_name", "log2FoldChange", "padj", "sig", "direction_1", "dir_accession", "cancer")]
}))
names(tum)[1] <- "gene"
cat("TCGA cancers:", length(unique(tum$cancer)), "\n")
per_source <- function(D, src_col, analysis) bind_rows(lapply(unique(D[[src_col]]), function(s) {
  d <- D[D[[src_col]] == s, ]; tested <- unique(d$gene); sig <- d[d$sig == "y", ]
  bind_rows(lapply(names(BGS(tested)), function(bt) { bg <- BGS(tested)[[bt]]
    o <- suppressWarnings(overlap_function(df_1 = sig, gene_col_1 = "gene", group_col_1 = "dir_accession",
                                           carry_col_1 = "direction_1", df_2 = metabolism_genes, gene_col_2 = "gene",
                                           group_col_2 = "class", background = bg))
    data.frame(analysis = analysis, bg_type = bt, background_n = length(bg), group = s, direction = o$direction_1, o) })) }))
TC <- per_source(tum, "cancer", "tcga"); TC <- add_corr(TC, TC$bg_type)
save_csv(TC, file_name = "metabolic_overlap_tcga", path = OUT)

# ---- (2) per-inducer senescence arms ----
ng <- read.csv(SEN$nongenotoxic)
ng$dir_accession <- paste(ng$arm, ng$direction_1, sep = "_")
NG <- per_source(ng, "arm", "inducer_arm"); NG <- add_corr(NG, NG$bg_type)
NG$class_arm <- ng$class[match(NG$group, ng$arm)]
save_csv(NG, file_name = "metabolic_overlap_inducer_arms", path = OUT)

# ---- (4) senescence vs tumour direction within metabolic genes ----
sen <- ad[ad$sig == "y" & ad$gene %in% metabolism_genes$gene, ]
sen$grp_dir <- paste(ARR[as.character(sen$group_1)], sen$direction_1)
CMP <- bind_rows(lapply(unique(tum$cancer), function(ca) {
  d <- tum[tum$cancer == ca, ]; bg <- intersect(intersect(unique(metabolism_genes$gene), unique(ad$gene)), unique(d$gene))
  o <- suppressWarnings(overlap_function(df_1 = sen, gene_col_1 = "gene", group_col_1 = "grp_dir",
                                         df_2 = d[d$sig == "y", ], gene_col_2 = "gene", group_col_2 = "direction_1", background = bg))
  data.frame(cancer = ca, background_n = length(bg), o) }))
CMP <- add_corr(CMP, rep("all", nrow(CMP)))
CMP$sen_group <- sub(" .*", "", CMP$grp_dir); CMP$sen_dir <- sub(".* ", "", CMP$grp_dir)
CMP$relation <- ifelse(CMP$sen_dir == CMP$direction_1, "same direction", "opposite direction")
save_csv(CMP, file_name = "metabolic_genes_senescence_vs_tumour", path = OUT)

# per-gene table: metabolic class, direction in each senescence group, n cancers up / down
gt <- metabolism_genes |> distinct(gene, class) |>
  left_join(tum |> group_by(gene) |> summarise(n_cancer_up = sum(sig == "y" & direction_1 == "up"),
                                               n_cancer_down = sum(sig == "y" & direction_1 == "down"),
                                               median_lfc_tumour = median(log2FoldChange), .groups = "drop"), by = "gene")
for (g in names(ARR)) { x <- ad[ad$group_1 == g, ]; gt[[ARR[[g]]]] <- with(x[match(gt$gene, x$gene), ], ifelse(sig == "y", ifelse(direction_1 == "up", "up", "down"), ifelse(is.na(sig), NA, "ns"))) }
save_csv(gt, file_name = "metabolic_genes_direction_table", path = OUT)

# ---- figures (house summarise_overlaps tile; label_cq = FALSE because tumour rows are not CS/CQ) ----
nice_cls <- function(x) trimws(sub("Oxidative phosphorylation", "OXPHOS", x))
Ap <- A[A$bg_type == "SI37_pathway_genes", ]
Ap$y <- paste0("Senescence ", Ap$group, " ", ifelse(Ap$direction == "up", "U", "D"))
Tp <- TC[TC$bg_type == "SI37_pathway_genes", ]; Tp$y <- paste0(Tp$group, " ", ifelse(Tp$direction == "up", "U", "D"))
both <- bind_rows(Ap[, c("y", "class", "odds", "adj_BH", "direction")], Tp[, c("y", "class", "odds", "adj_BH", "direction")])
both$x <- nice_cls(both$class)
lv <- c(rev(sort(unique(Tp$y[Tp$direction == "down"]))), paste0("Senescence ", c("SSCQ", "CICQ", "OIS", "SIPS", "RS"), " D"),
        rev(sort(unique(Tp$y[Tp$direction == "up"]))), paste0("Senescence ", c("SSCQ", "CICQ", "OIS", "SIPS", "RS"), " U"))
both$y <- factor(both$y, levels = lv)
p <- summarise_overlaps(both, x = "x", y = "y", pval_col = "adj_BH", xlab = "Metabolic pathways (SI 37; background: SI 37 pathway genes tested)",
                        ylab = "DEGs (U = up, D = down)", add_line = FALSE, label_cq = FALSE)
save_p(p, file_name = "fig_metabolic_senescence_vs_tcga", save_dir = OUT, p_width = 9, p_height = 12)
Np <- NG[NG$bg_type == "SI37_pathway_genes", ]; Np$x <- nice_cls(Np$class)
Np$y <- factor(paste0(Np$group, " ", ifelse(Np$direction == "up", "U", "D")))
q <- summarise_overlaps(Np, x = "x", y = "y", pval_col = "adj_BH", xlab = "Metabolic pathways (SI 37)", ylab = "Inducer arm DEGs",
                        add_line = FALSE, label_cq = FALSE)
save_p(q, file_name = "fig_metabolic_inducer_arms", save_dir = OUT, p_width = 9, p_height = 11)
Cp <- CMP; Cp$x <- Cp$cancer; Cp$y <- factor(paste0(Cp$grp_dir, " vs tumour ", Cp$direction_1))
r <- summarise_overlaps(Cp, x = "x", y = "y", pval_col = "adj_BH", xlab = "Cancer", ylab = "Senescence metabolic DEGs vs tumour DEGs",
                        add_line = FALSE, label_cq = FALSE)
save_p(r, file_name = "fig_metabolic_genes_senescence_vs_tumour", save_dir = OUT, p_width = 10, p_height = 9)

# ---- console summaries ----
cat("\n== TCGA x metabolic class (primary background): n cancers sig over / under (BH < 0.05), median log2 OR ==\n")
s <- Tp |> group_by(class, direction) |> summarise(over = sum(adj_BH < 0.05 & odds > 1), under = sum(adj_BH < 0.05 & odds < 1),
                                                  med_log2OR = round(median(log2(pmax(odds, 1e-3))), 2), .groups = "drop")
print(as.data.frame(s), row.names = FALSE)
cat("\n== Senescence arrest x metabolic class (", if (NOPROLIF) "recomputed, proliferation genes removed" else "house output", ", primary background): odds (* BH < 0.05) ==\n", sep = "")
Ap$cell <- sprintf("%.2f%s", Ap$odds, ifelse(Ap$adj_BH < 0.05, "*", ""))
print(reshape(Ap[, c("group", "direction", "class", "cell")], idvar = c("group", "direction"), timevar = "class", direction = "wide"), row.names = FALSE)
cat("\n== Inducer arms: sig classes (BH < 0.05, primary background) ==\n")
print(Np |> filter(adj_BH < 0.05) |> transmute(arm = group, class_arm, direction, class, actual, odds = round(odds, 2), adj_BH = signif(adj_BH, 2)) |> arrange(class, direction), row.names = FALSE)
cat("\n== Metabolic genes: senescence vs tumour direction; n cancers sig over / under (BH < 0.05) ==\n")
print(as.data.frame(CMP |> group_by(grp_dir, relation, tumour = direction_1) |>
        summarise(over = sum(adj_BH < 0.05 & odds > 1), under = sum(adj_BH < 0.05 & odds < 1),
                  med_log2OR = round(median(log2(pmax(odds, 1e-3))), 2), .groups = "drop")), row.names = FALSE)
cat("done\n")

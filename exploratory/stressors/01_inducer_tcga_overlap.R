# Per-stressor senescence DEGs vs TCGA tumour DEGs (exploratory, 2026-11-07).
# Question: the meta-analysis groups pool inducers (SIPS = IR + etoposide + bleomycin + H2O2). Does each stressor
# overlap tumour DEGs the same way (opposite direction), and does any share tumour direction as OIS does?
# Sources (senescence_stress outputs, read-only; IR time course from the per-cell-type temporal/<cell>/degs.csv files):
#   meta    arrest_degs_final_RERUN.csv, the 5 meta-analysis groups (reference; the main pipeline uses the published table)
#   arm     rerun_outputs/nongenotoxic/degs_nongenotoxic.csv: 21 inducer arms (genotoxic, CDK4/6i, MDM2i, quiescence, OIS)
#   IR_time work/results/auc/temporal/<cell type>/degs.csv: ERP021140 IR 10 Gy, 3 cell types x 4/10/20 days
#   GSE60340 rerun_outputs/robustness/damage_not_p53/gse60340_degs.csv: p53-mutant LFS fibroblasts, 5 contrasts
# Tumour DEGs: pancancer-deg pipeline (covariate-adjusted), protein-coding, padj < 0.05 and |log2FC| > log2(1.5).
# Test: house overlap_function (GeneOverlap + two-sided Fisher), one call per source set x cancer, df_1 = source DEGs
# by direction, df_2 = tumour DEGs by direction, background = genes tested in both. BH across all tests.
# Run from the pancancer-deg repo root. OUT: results/exploratory/stressors/
SA <- "/Users/ro/Academic/10_Research/senescence_stress/modules/systems_analysis_arrest"
BULK <- file.path(SA, "bulk_RNA_seq"); PD <- getwd()
OUT <- file.path(PD, "results/exploratory/stressors"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(BULK, ".Rlib"), .libPaths())); Sys.setenv(BULK_RNASEQ_DIR = SA)
source(file.path(BULK, "R/config.R")); source(file.path(BULK, "R/functions.R"))
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
LFC_T <- log2(1.5)
std <- function(d, src, set, gene, lfc, padj, sig = NULL, cls = NA) {
  s <- if (is.null(sig)) ifelse(!is.na(d[[padj]]) & d[[padj]] < .05 & abs(d[[lfc]]) > LFC_T, "y", "n") else d[[sig]]
  data.frame(source = src, set = set, class = cls, gene = d[[gene]], sig = s, dir = ifelse(d[[lfc]] > 0, "up", "down"))
}
ARR <- c(Contact_inhibited_CQ = "CICQ", Serum_starved_CQ = "SSCQ", Replicative_CS = "RS", Stress_induced_CS = "SIPS", Oncogene_induced_CS = "OIS")
ad <- read.csv(file.path(RERUN_DIR, "arrest_degs_final_RERUN.csv"))
ng <- read.csv(file.path(RERUN_DIR, "nongenotoxic", "degs_nongenotoxic.csv"))
# per-cell-type files: the combined temporal/all_time_degs.csv is corrupt on disk (33 MB, almost all NUL bytes; 2026-11-07)
ir <- do.call(rbind, lapply(c("Fibroblast", "Keratinocyte", "Melanocyte"), function(ct) {
  x <- read.csv(file.path("/Users/ro/Academic/10_Research/senescence_stress/work/results/auc/temporal", ct, "degs.csv"))
  x$cell_type <- ct; x[, c("gene", "log2FoldChange", "padj", "sig", "group_1", "cell_type")] }))
lf <- read.csv(file.path(RERUN_DIR, "robustness", "damage_not_p53", "gse60340_degs.csv"))
S <- bind_rows(
  bind_rows(lapply(names(ARR), function(g) std(ad[ad$group_1 == g, ], "meta", ARR[[g]], "gene", "log2FoldChange", "padj", "sig", "meta-analysis group"))),
  bind_rows(lapply(unique(ng$arm), function(a) { x <- ng[ng$arm == a, ]; std(x, "arm", a, "gene", "log2FoldChange", "padj", "sig", x$class[1]) })),
  bind_rows(lapply(split(ir, paste(ir$cell_type, ir$group_1)), function(x) std(x, "IR_time", paste("IR", x$cell_type[1], sub("_days", "d", x$group_1[1])), "gene", "log2FoldChange", "padj", "sig", "genotoxic (time course)"))),
  bind_rows(lapply(unique(lf$contrast), function(k) std(lf[lf$contrast == k, ], "GSE60340", paste("LFS", k), "gene", "log2FC", "padj", "sig", "p53-mutant")))
)
# 2026-11-09: house stress-response modules (external_inputs/stress_modules/v3; "full" version). Module genes = DEGs
# (direction = sign of max_log2FC); background = the module's tested genes. Acute (GSE241261 etc.) and chronic /
# in-senescent (antimycin A, 2-DG, glucocorticoid, 21% O2, IL-1b 14 d, heat shock, UPR). Separate BH family.
MD <- file.path(BULK, "external_inputs", "stress_modules", "v3")
mod_src <- function(mods, tested, fam) bind_rows(lapply(unique(mods$module), function(m) {
  x <- mods[mods$module == m, ]; x <- x[!duplicated(x$gene), ]; tg <- unique(c(tested$gene[tested$module == m], x$gene))
  data.frame(source = "module", set = paste0(x$module_label[1], " [", fam, "]"), class = fam, gene = tg,
             sig = ifelse(tg %in% x$gene, "y", "n"), dir = ifelse(x$max_log2FC[match(tg, x$gene)] > 0, "up", "down") %||% "up")
}))
`%||%` <- function(a, b) ifelse(is.na(a), b, a)
am <- read.csv(file.path(MD, "stress_modules.csv")); am <- am[am$version == "full", ]
cm <- read.csv(file.path(MD, "chronic_modules.csv")); cm <- cm[cm$version == "full", ]
S <- bind_rows(S, mod_src(am, read.csv(file.path(MD, "stress_module_tested_genes.csv")), "acute"),
               mod_src(cm, read.csv(file.path(MD, "chronic_module_tested_genes.csv")), "chronic"))
S <- S[!is.na(S$gene) & S$gene != "", ]; S <- S[!duplicated(S[, c("set", "gene")]), ]
tf <- list.files(file.path(PD, "results/tables/TCGA_deg"), pattern = "_DGE_Raw_Results_Tumor_vs_Normal_Adjusted.csv$", full.names = TRUE)
TU <- bind_rows(lapply(tf, function(f) { d <- read.csv(f); d <- d[d$gene_type == "protein_coding" & !is.na(d$gene_name), ]
  d <- d[!duplicated(d$gene_name), ]; x <- std(d, "tcga", sub("^TCGA-", "", sub("_DGE.*", "", basename(f))), "gene_name", "log2FoldChange", "padj"); x }))
cat("senescence sets:", length(unique(S$set)), " cancers:", length(unique(TU$set)), "\n")
print(S |> group_by(source, set) |> summarise(tested = n(), up = sum(sig == "y" & dir == "up"), down = sum(sig == "y" & dir == "down"), .groups = "drop") |> as.data.frame())

res <- bind_rows(lapply(unique(S$set), function(s) { a <- S[S$set == s, ]
  bind_rows(lapply(unique(TU$set), function(ca) { b <- TU[TU$set == ca, ]
    bg <- intersect(a$gene, b$gene)
    o <- suppressWarnings(overlap_function(df_1 = a[a$sig == "y", ], gene_col_1 = "gene", group_col_1 = "dir",
                                           df_2 = b[b$sig == "y", ], gene_col_2 = "gene", group_col_2 = "dir", background = bg))
    names(o)[names(o) == "dir"] <- "sen_dir"; names(o)[names(o) == "dir.1"] <- "tumour_dir"
    data.frame(source = a$source[1], set = s, class = a$class[1], cancer = ca, o) })) }))
res$adj_BH <- ave(res$pval, res$source == "module", FUN = function(p) p.adjust(p, "BH"))  # modules: own BH family (2026-11-09)
res$relation <- paste0("sen ", res$sen_dir, " / tumour ", res$tumour_dir)
save_csv(res, file_name = "stressor_tcga_overlaps", path = OUT)

sm <- res |> group_by(source, class, set, relation) |>
  summarise(enriched = sum(adj_BH < .05 & odds > 1), depleted = sum(adj_BH < .05 & odds < 1),
            med_log2OR = median(log2(pmax(odds, 1e-3))), .groups = "drop")
save_csv(sm, file_name = "stressor_tcga_summary", path = OUT)
w <- sm |> mutate(cell = sprintf("%d/%d", enriched, depleted)) |> dplyr::select(source, class, set, relation, cell) |>
  tidyr::pivot_wider(names_from = relation, values_from = cell)
cat("\nn cancers (of 17) enriched/depleted, BH < 0.05\n"); print(as.data.frame(w), right = FALSE, row.names = FALSE)

# summary figure: median log2 OR across cancers, label = n cancers enriched / depleted (no house equivalent for this summary)
sm$lab <- sprintf("%d/%d", sm$enriched, sm$depleted)
sm$set <- factor(sm$set, levels = rev(unique(sm$set[order(sm$source, sm$class, sm$set)])))
p <- ggplot(sm, aes(relation, set, fill = pmax(pmin(med_log2OR, 2), -2))) + geom_tile(colour = "white") +
  geom_text(aes(label = lab), size = 2.6) +
  scale_fill_gradient2(low = "#c0392b", mid = "#f7f7f7", high = "#2c7bb6", name = "median log2 OR\n(17 cancers)") +
  facet_grid(source ~ ., scales = "free_y", space = "free_y") + labs(x = NULL, y = NULL,
  caption = "Label: n cancers enriched / depleted (two-sided Fisher, BH < 0.05). Blue = enrichment, red = depletion.") +
  theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 30, hjust = 1), strip.text.y = element_text(angle = 0))
save_p(p, file_name = "fig_stressor_tcga_overlaps", save_dir = OUT, p_width = 7.5, p_height = 11)
cat("done\n")

# House Hallmark + lysosomal over-representation in TCGA tumour DEGs and per-inducer arms, compared with senescence
# (exploratory, 2026-11-09).
# COPY of family 1 (hallmark_arrest) of senescence_stress/.../meta_analysis/17_stress_module_overlap_v3.R (read-only):
# same sets (rerun_outputs/stress_response_pathways_RERUN.csv: 49 Hallmark + "Lysosomal Genes"), same call
# (house overlap_function, df_1 = pathways by pathway, df_2 = significant DEGs by dir_accession, background = all genes
# tested), same correction (add_corr from 17: BH within family). WHAT CHANGED (inputs only):
#   TCGA tumour-vs-normal DEGs (protein-coding; padj < 0.05, |log2FC| > log2(1.5)), one call per cancer (background =
#   genes tested in that cancer), family = all cancers; per-inducer arms (degs_nongenotoxic.csv), one call per arm,
#   family = all arms. Senescence arrest results are READ from the house output (stress_modules_v3/hallmark_overlap_arrest.csv).
# Run from the pancancer-deg repo root. OUT: results/exploratory/hallmarks/
SA <- "/Users/ro/Academic/10_Research/senescence_stress/modules/systems_analysis_arrest"
BULK <- file.path(SA, "bulk_RNA_seq"); PD <- getwd()
OUT <- file.path(PD, "results/exploratory/hallmarks"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(BULK, ".Rlib"), .libPaths())); Sys.setenv(BULK_RNASEQ_DIR = SA)
source(file.path(BULK, "R/config.R")); source(file.path(BULK, "R/functions.R"))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
ex <- parse(file.path(BULK, "meta_analysis/17_stress_module_overlap_v3.R"))
fn <- sapply(ex, function(e) is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% c("add_corr", "ARR"))
for (e in ex[fn]) eval(e)
hall <- read.csv(file.path(RERUN_DIR, "stress_response_pathways_RERUN.csv"))
A <- read.csv(file.path(RERUN_DIR, "stress_modules_v3", "hallmark_overlap_arrest.csv"))
LFC_T <- log2(1.5)
tf <- list.files(file.path(PD, "results/tables/TCGA_deg"), pattern = "_DGE_Raw_Results_Tumor_vs_Normal_Adjusted.csv$", full.names = TRUE)
tum <- bind_rows(lapply(tf, function(f) { d <- read.csv(f); d <- d[d$gene_type == "protein_coding" & !is.na(d$gene_name) & d$gene_name != "", ]
  d <- d[!duplicated(d$gene_name), ]; ca <- sub("^TCGA-", "", sub("_DGE.*", "", basename(f)))
  data.frame(gene = d$gene_name, sig = ifelse(!is.na(d$padj) & d$padj < .05 & abs(d$log2FoldChange) > LFC_T, "y", "n"),
             direction_1 = ifelse(d$log2FoldChange > 0, "up", "down"), group_1 = ca,
             dir_accession = paste(ca, ifelse(d$log2FoldChange > 0, "up", "down"), sep = "_")) }))
ng <- read.csv(file.path(RERUN_DIR, "nongenotoxic", "degs_nongenotoxic.csv")); ng$group_1 <- ng$arm
ng$dir_accession <- paste(ng$arm, ng$direction_1, sep = "_")
per_source <- function(D, analysis) bind_rows(lapply(unique(D$group_1), function(s) { d <- D[D$group_1 == s, ]
  o <- suppressWarnings(overlap_function(df_1 = hall, df_2 = d[d$sig == "y", ], gene_col_1 = "genes", gene_col_2 = "gene",
                                         group_col_2 = "dir_accession", group_col_1 = "pathway",
                                         carry_col_2 = c("group_1", "direction_1"), background = unique(d$gene)))
  data.frame(analysis = analysis, group = o$group_1, direction = o$direction_1, o) }))
TC <- per_source(tum, "tcga"); TC <- add_corr(TC, rep("tcga", nrow(TC))); save_csv(TC, file_name = "hallmark_overlap_tcga", path = OUT)
NG <- per_source(ng, "inducer_arm"); NG <- add_corr(NG, rep("arm", nrow(NG))); save_csv(NG, file_name = "hallmark_overlap_inducer_arms", path = OUT)

# comparison table: per pathway, senescence groups enriched (BH < 0.05) by direction, and n cancers over / under
nice <- function(x) sub("^HALLMARK ", "", x)
sen <- A |> filter(adj_BH < 0.05, odds > 1) |> group_by(pathway, direction) |> summarise(sen = paste(group, collapse = " "), .groups = "drop") |>
  pivot_wider(names_from = direction, values_from = sen, names_prefix = "senescence_enriched_")
tcg <- TC |> group_by(pathway, direction) |> summarise(over = sum(adj_BH < .05 & odds > 1), under = sum(adj_BH < .05 & odds < 1), .groups = "drop") |>
  mutate(cell = paste0(over, "/", under)) |> dplyr::select(pathway, direction, cell) |>
  pivot_wider(names_from = direction, values_from = cell, names_prefix = "tumour_over_under_")
cmp <- full_join(sen, tcg, by = "pathway") |> mutate(pathway = nice(pathway)) |> arrange(pathway)
save_csv(cmp, file_name = "hallmark_senescence_vs_tumour_table", path = OUT)
options(width = 250); print(as.data.frame(cmp), row.names = FALSE, right = FALSE)

# figure (house summarise_overlaps tile; label_cq = FALSE)
Ap <- A; Ap$y <- paste0("Senescence ", Ap$group, " ", ifelse(Ap$direction == "up", "U", "D"))
Tp <- TC; Tp$y <- paste0(Tp$group, " ", ifelse(Tp$direction == "up", "U", "D"))
both <- bind_rows(Ap[, c("y", "pathway", "odds", "adj_BH", "direction")], Tp[, c("y", "pathway", "odds", "adj_BH", "direction")])
both$x <- nice(both$pathway)
lv <- c(rev(sort(unique(Tp$y[Tp$direction == "down"]))), paste0("Senescence ", c("SSCQ", "CICQ", "OIS", "SIPS", "RS"), " D"),
        rev(sort(unique(Tp$y[Tp$direction == "up"]))), paste0("Senescence ", c("SSCQ", "CICQ", "OIS", "SIPS", "RS"), " U"))
both$y <- factor(both$y, levels = lv)
p <- summarise_overlaps(both, x = "x", y = "y", pval_col = "adj_BH", xlab = "House Hallmark + lysosomal sets", ylab = "DEGs (U = up, D = down)",
                        add_line = FALSE, label_cq = FALSE)
save_p(p, file_name = "fig_hallmark_senescence_vs_tcga", save_dir = OUT, p_width = 16, p_height = 12)
cat("done\n")

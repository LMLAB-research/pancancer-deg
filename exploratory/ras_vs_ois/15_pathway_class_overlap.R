# Which RIS gene class holds the genes of each pathway (TGF-beta, cGAS-STING, interferon, NF-kB, ...)?
# House overlap_function() (GeneOverlap + two-sided Fisher; adj = Holm, the house default), run per OIS
# direction with the classified OIS DEGs of that direction as background, so each test asks whether a
# class holds more or fewer of the pathway's OIS genes than the other classes. All 50 Hallmark sets plus the
# two cGAS-STING sets. Run from the pancancer-deg repo root.
source(".claude/worktrees/two-sided-fisher/exploratory/ras_vs_ois/house_env.R")
source(".claude/worktrees/two-sided-fisher/exploratory/ras_vs_ois/pathway_sets.R")
suppressMessages({ library(dplyr); library(tidyr); library(GeneOverlap) })
out <- "results/exploratory/ras_vs_ois/pathways"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cls <- read.csv("results/exploratory/ras_vs_ois/strict/ois_gene_classes.csv")
cls <- cls[cls$arrest_class != "unassessable", ]
# extra grouping: NF-kB/BRD4-dependent OIS-up genes, as in the other analyses
nf <- cls[cls$OIS_direction == "UP" & cls$nfkb_brd4 == "dependent", ]; nf$arrest_class <- "NF-kB/BRD4-dependent"
cl2 <- rbind(cls, nf)

.libPaths(c("/Users/ro/Academic/10_Research/senescence_stress/modules/systems_analysis_arrest/bulk_RNA_seq/.Rlib", .libPaths()))
h <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
ps <- pathway_sets()
pw <- rbind(data.frame(pathway = sub("^HALLMARK_", "H_", h$gs_name), gene = h$gene_symbol),
            data.frame(pathway = "cGAS_STING_KEGG", gene = ps$P_cGAS_STING_KEGG),
            data.frame(pathway = "STING_Reactome", gene = ps$P_STING_Reactome)) |> distinct()

res <- do.call(rbind, lapply(c("UP", "DOWN"), function(d) {
  bg <- unique(cls$gene[cls$OIS_direction == d])
  r <- overlap_function(df_1 = cl2[cl2$OIS_direction == d, ], df_2 = pw[pw$gene %in% bg, ],
                        gene_col_1 = "gene", gene_col_2 = "gene",
                        group_col_1 = "arrest_class", group_col_2 = "pathway", background = bg)
  r$OIS_direction <- d; r
}))
res$padj_BH <- ave(res$pval, res$OIS_direction, res$arrest_class, FUN = function(p) p.adjust(p, "BH"))
write.csv(res, file.path(out, "class_pathway_overlap.csv"), row.names = FALSE)

focus <- c("H_TGF_BETA_SIGNALING", "cGAS_STING_KEGG", "STING_Reactome", "H_INTERFERON_ALPHA_RESPONSE",
           "H_INTERFERON_GAMMA_RESPONSE", "H_TNFA_SIGNALING_VIA_NFKB", "H_IL6_JAK_STAT3_SIGNALING",
           "H_INFLAMMATORY_RESPONSE", "H_P53_PATHWAY", "H_NOTCH_SIGNALING", "H_KRAS_SIGNALING_UP",
           "H_KRAS_SIGNALING_DN", "H_EPITHELIAL_MESENCHYMAL_TRANSITION", "H_E2F_TARGETS")
cat("Focus pathways: overlap / pathway genes among OIS DEGs of that direction, odds ratio, BH within class\n")
f <- res |> filter(pathway %in% focus) |>
  mutate(cell = sprintf("%d/%d OR %.2f%s", as.integer(actual), as.integer(n2), odds,
                        ifelse(padj_BH < 0.05, "*", ""))) |>
  dplyr::select(OIS_direction, pathway, arrest_class, cell) |>
  pivot_wider(names_from = arrest_class, values_from = cell)
print(as.data.frame(f), right = FALSE)
cat("\nAll Hallmarks with BH < 0.05 in any class:\n")
print(res |> filter(padj_BH < 0.05) |> arrange(OIS_direction, arrest_class, pval) |>
        transmute(OIS_direction, arrest_class, pathway, actual, expected = round(expected, 1), odds = round(odds, 2),
                  padj_BH = signif(padj_BH, 2)), row.names = FALSE)
# gene lists for the cGAS-STING / TGF-beta / IFN pathways, per class
gl <- cls |> dplyr::select(gene, OIS_direction, arrest_class, nfkb_brd4) |>
  inner_join(pw |> filter(pathway %in% focus[1:5]), by = "gene", relationship = "many-to-many") |>
  arrange(pathway, OIS_direction, arrest_class, gene)
write.csv(gl, file.path(out, "pathway_genes_by_class.csv"), row.names = FALSE)
cat("\nTGF-beta / cGAS-STING / IFN genes among OIS DEGs, by class:\n")
print(gl |> group_by(pathway, OIS_direction, arrest_class) |> summarise(genes = paste(gene, collapse = " "), .groups = "drop") |> as.data.frame(), right = FALSE)

# ---- 2026-11-09: same test with the HOUSE pathway sets (49 Hallmark + Lysosomal Genes; house_sets()) ----
hs <- house_sets()
pwh <- do.call(rbind, lapply(names(hs), function(k) data.frame(pathway = k, gene = hs[[k]])))
resh <- do.call(rbind, lapply(c("UP", "DOWN"), function(d) {
  bg <- unique(cls$gene[cls$OIS_direction == d])
  r <- overlap_function(df_1 = cl2[cl2$OIS_direction == d, ], df_2 = pwh[pwh$gene %in% bg, ],
                        gene_col_1 = "gene", gene_col_2 = "gene",
                        group_col_1 = "arrest_class", group_col_2 = "pathway", background = bg)
  r$OIS_direction <- d; r
}))
resh$padj_BH <- ave(resh$pval, resh$OIS_direction, resh$arrest_class, FUN = function(p) p.adjust(p, "BH"))
write.csv(resh, file.path(out, "class_pathway_overlap_house.csv"), row.names = FALSE)
cat("\n== HOUSE sets: BH < 0.05 in any class ==\n")
print(resh |> filter(padj_BH < 0.05) |> arrange(OIS_direction, arrest_class, pval) |>
        transmute(OIS_direction, arrest_class, pathway, actual, expected = round(expected, 1), odds = round(odds, 2),
                  padj_BH = signif(padj_BH, 2)), row.names = FALSE)
cat("\n== HOUSE Lysosomal Genes by class ==\n")
print(resh |> filter(pathway == "HS_LYSOSOMAL_GENES") |> transmute(OIS_direction, arrest_class, actual, n2, odds = round(odds, 2), padj_BH = signif(padj_BH, 2)), row.names = FALSE)

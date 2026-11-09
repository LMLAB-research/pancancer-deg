# Metabolic gene sets for scoring (single-cell tests and TCGA tumour scores), written to
# results/exploratory/metabolism/metabolic_score_sets.csv and read by ras_vs_ois/pathway_sets.R. Exploratory, 2026-11-07.
#   P_GLYC_SI37     SI Table 37 (keep == TRUE) class "Glycolysis" (the house curated enzyme set)
#   P_GLYC_HALLMARK Hallmark glycolysis (sensitivity; broader, includes non-enzyme genes)
#   P_GAG_SYNTH     SI Table 37 class "Hexosamine pathway" genes that make glycosaminoglycan chains (chondroitin /
#                   dermatan / heparan sulfate: linker, polymerases, sulfotransferases, epimerases). The rest of the
#                   class (hexosamine-pathway proper and sugar interconversion) is P_HEXOSAMINE_CORE.
#   The GAG / core split is a gene-name rule (Claude, 2026-11-07), printed below for checking.
# Run from the pancancer-deg repo root.
SA <- "/Users/ro/Academic/10_Research/senescence_stress/modules/systems_analysis_arrest"
.libPaths(c(file.path(SA, "bulk_RNA_seq/.Rlib"), .libPaths()))
suppressPackageStartupMessages(library(readxl))
si <- as.data.frame(read_excel(file.path(SA, "Final/final/si_tables.xlsx"), sheet = "SI Table 37", skip = 1, na = c("", "NA")))
si <- si[si$keep %in% TRUE, ]
hex <- unique(si$gene[si$class == "Hexosamine pathway"])
gag_rx <- "^(CHST|CHPF|CHSY|CSGALNACT|B3GAT[13]|B4GALT7|B3GALT6|XYLT|EXT|NDST|HS[236]ST|GLCE|UST|DSE|DSEL)"
h <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
sets <- list(P_GLYC_SI37 = unique(si$gene[si$class == "Glycolysis"]),
             P_GLYC_HALLMARK = unique(h$gene_symbol[h$gs_name == "HALLMARK_GLYCOLYSIS"]),
             P_GAG_SYNTH = sort(hex[grepl(gag_rx, hex)]),
             P_HEXOSAMINE_CORE = sort(hex[!grepl(gag_rx, hex)]))
for (k in names(sets)) cat(k, length(sets[[k]]), ":", paste(head(sets[[k]], 60), collapse = " "), "\n")
out <- "results/exploratory/metabolism"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, lapply(names(sets), function(k) data.frame(set = k, gene = sets[[k]]))),
          file.path(out, "metabolic_score_sets.csv"), row.names = FALSE)

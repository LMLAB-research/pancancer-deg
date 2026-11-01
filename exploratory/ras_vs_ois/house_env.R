# Load Roberto's house bulk RNA-seq functions READ-ONLY from senescence_stress
# (functions.R and its package library are only read, never modified).
SS_BULK <- "/Users/ro/Academic/10_Research/senescence_stress/modules/systems_analysis_arrest"
.libPaths(c(file.path(SS_BULK, "bulk_RNA_seq/.Rlib"), .libPaths()))
source(file.path(SS_BULK, "HANDOFF_2026-08-31/bulk_RNA_seq/R/functions.R"))

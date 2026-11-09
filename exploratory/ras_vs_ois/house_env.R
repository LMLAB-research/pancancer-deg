# Load Roberto's house bulk RNA-seq functions (senescence_stress HANDOFF_2026-08-31/bulk_RNA_seq/R/functions.R; copy in
# exploratory/inputs/house_code, md5 in inputs/MANIFEST.tsv). Run from the pancancer-deg repo root.
if (!exists("inp")) source("exploratory/paths.R")
source_house_bulk_handoff()

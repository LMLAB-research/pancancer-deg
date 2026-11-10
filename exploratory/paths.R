# Paths for the exploratory scripts (2026-11-09). Source from the pancancer-deg repo root:
#   source("exploratory/paths.R")
# - Inputs from other projects (senescence_stress tables, house code, purity, proliferation list) are read from the
#   repo copies in exploratory/inputs/ (see inputs/MANIFEST.tsv for the originals).
# - Intermediates the pipeline regenerates (TCGA DEG tables, within-study contrast DEGs) are read from results/ when
#   they exist there, otherwise from the copies in exploratory/inputs/.
# - Only the multi-GB single-cell objects and raw data that are not in the repo use SENESCENCE_STRESS_DIR
#   (default: Roberto's machine); see inputs/DATA.md.
# - R package libraries of the senescence_stress project are added to .libPaths() only if they exist.
EXP <- "exploratory"
INP <- file.path(EXP, "inputs")
if (!dir.exists(INP)) stop("Run from the pancancer-deg repo root (exploratory/inputs not found from ", getwd(), ")")
inp <- function(...) file.path(INP, ...)
SS_DIR <- Sys.getenv("SENESCENCE_STRESS_DIR", "/Users/ro/Academic/10_Research/senescence_stress")
add_lib <- function(p) if (dir.exists(p)) .libPaths(c(p, .libPaths()))
BULK_LIB <- file.path(SS_DIR, "modules/systems_analysis_arrest/bulk_RNA_seq/.Rlib")
SC_LIB <- file.path(SS_DIR, "modules/scrna_seq/GSE226225/final_analysis/.Rlib/main")

# house code (copies in inputs/house_code)
source_house_bulk <- function() { add_lib(BULK_LIB); source(inp("house_code", "bulk_functions.R")) }
source_house_bulk_handoff <- function() { add_lib(BULK_LIB); source(inp("house_code", "HANDOFF_2026-08-31_bulk_functions.R")) }
source_ucell <- function() source(inp("house_code", "scutils_ucell_parallel.R"))
house_add_corr <- function(env = parent.frame()) {   # add_corr() and ARR from the house 17_stress_module_overlap_v3.R
  ex <- parse(inp("house_code", "17_stress_module_overlap_v3.R"))
  fn <- sapply(ex, function(e) is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% c("add_corr", "ARR"))
  for (e in ex[fn]) eval(e, envir = env)
}

# senescence_stress tables
SEN <- list(
  arrest_rerun = inp("senescence_stress", "arrest_degs_final_RERUN.csv.gz"),
  arrest_published = inp("senescence_stress", "arrest_degs_final_published.csv.gz"),
  nongenotoxic = inp("senescence_stress", "degs_nongenotoxic.csv.gz"),
  gse60340 = inp("senescence_stress", "gse60340_degs.csv.gz"),
  hallmark_sets = inp("senescence_stress", "stress_response_pathways_RERUN.csv"),
  hallmark_arrest = inp("senescence_stress", "hallmark_overlap_arrest.csv"),
  metabolic_arrest = inp("senescence_stress", "metabolic_overlap_arrest_SI37.csv"),
  si18 = inp("senescence_stress", "SI_Table_18.csv"),
  si37 = inp("senescence_stress", "SI_Table_37.csv"),
  ris_groups = inp("senescence_stress", "cell_scores_GSE245807_RIS_groups.csv.gz"))
ir_time <- function(cell_type) inp("senescence_stress", paste0("temporal_IR_", cell_type, "_degs.csv.gz"))
stress_module <- function(name) inp("senescence_stress", paste0("stress_modules_v3_", name, ".csv.gz"))   # stress_modules, chronic_modules, *_tested_genes
read_si <- function(key) { d <- read.csv(SEN[[key]], check.names = FALSE); if ("keep" %in% names(d)) d$keep <- as.logical(d$keep); d }
PROLIF <- inp("other", "prolif_genes.txt")
PURITY_FILE <- inp("other", "TCGA_mastercalls.abs_tables_JSedit.fixed.txt.gz")

# regenerable intermediates: results/ if present, else the repo copy
tcga_deg_files <- function() {
  f <- list.files("results/tables/TCGA_deg", pattern = "_DGE_Raw_Results_Tumor_vs_Normal_Adjusted.csv$", full.names = TRUE)
  if (length(f)) f else list.files(inp("tcga_degs"), pattern = "_DGE_Raw_Results_Tumor_vs_Normal_Adjusted.csv.gz$", full.names = TRUE)
}
contrast_deg <- function(name) {
  f <- file.path("results/exploratory/ras_vs_ois/contrast_degs", paste0(name, ".csv"))
  if (file.exists(f)) f else inp("contrast_degs", paste0(name, ".csv.gz"))
}
overlap_results <- function() {   # two-sided overlap pipeline results: 2026-11-01 run if present, else run_both.sh output
  for (d in c("results/exploratory/two_sided_2026-11-01", Sys.getenv("RUN_BOTH_OUT", "results/exploratory/two_sided_rerun")))
    if (dir.exists(d)) return(d)
  stop("no overlap results: run bash exploratory/run_both.sh first")
}

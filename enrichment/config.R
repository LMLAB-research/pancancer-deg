signatures <- c("OIS", "RS", "SIPS", "CICQ", "SSCQ")
padj_threshold <- 0.05
lfc_threshold <- log2(1.5)
bh_cutoff <- 0.05
minimum_input_genes <- 5L
minimum_gene_set_size <- 10L
maximum_gene_set_size <- 500L
go_ontology <- "BP"
cancer_file_pattern <- "^TCGA-[A-Za-z0-9-]+_DGE_Raw_Results_Tumor_vs_Normal_Adjusted\\.csv$"
base_dir <- Sys.getenv("ENRICH_INPUT_ROOT", file.path(script_dir, "data"))
cancer_input_dir <- file.path(base_dir, "TCGA_deg")
senescence_input_dir <- file.path(base_dir, "meta_analysis")
overlap_root <- Sys.getenv("ENRICH_OVERLAP_ROOT", file.path(base_dir, "overlaps"))
legacy_ois_root <- Sys.getenv("ENRICH_LEGACY_OIS_ROOT", "")
output_root <- Sys.getenv("ENRICH_OUTPUT_ROOT", file.path(script_dir, "results"))
cache_root <- file.path(script_dir, "cache")
# The legacy OIS table is a full tested-gene table; opt in only after checking yours.
legacy_ois_full_table <- FALSE
if (file.exists(file.path(script_dir, "config.local.R"))) source(file.path(script_dir, "config.local.R"))

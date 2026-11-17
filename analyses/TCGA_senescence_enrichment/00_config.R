#### Paths ####
# Scripts are run from the repository root.

analysis_name <- "TCGA_senescence_enrichment"

# Settings shared with the overlap module (TCGA design, signatures, output
# folder) are read from its config, so the two modules always agree.
overlap_config <- new.env()
sys.source(file.path("analyses", "TCGA_senescence_overlap", "00_config.R"), envir = overlap_config)
cancer_design_type <- overlap_config$cancer_design_type
signatures <- overlap_config$expected_signatures
overlap_table_dir <- overlap_config$results_dir
tcga_table_dir <- overlap_config$cancer_input_dir

results_dir <- file.path("results", "tables", analysis_name)
figures_dir <- file.path("results", "figures", analysis_name)
gene_set_dir <- file.path("data", analysis_name, "gene_sets")

#### Gene sets ####

# All gene sets come from MSigDB (msigdbr), human gene symbols. Each entry:
# MSigDB collection and subcollection (NULL = whole collection).
gene_set_collections <- list(
  HALLMARK = list(collection = "H", subcollection = NULL),
  GO_BP = list(collection = "C5", subcollection = "GO:BP"),
  REACTOME = list(collection = "C2", subcollection = "CP:REACTOME"),
  KEGG = list(collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
  CGP = list(collection = "C2", subcollection = "CGP")
)

# Gene sets are tested if they have this many genes in the analysis universe.
min_gene_set_size <- 10
max_gene_set_size <- 500

#### Enrichment settings ####

# ORA: one-sided hypergeometric test; gene lists with fewer eligible genes are
# skipped. BH correction per gene list and collection.
ora_min_input_genes <- 5
fdr_threshold <- 0.05

# GSEA (fgsea, preranked): tumor genes ranked by shrunken log2 fold change of
# the TCGA design above, signature genes by their log2 fold change.
gsea_seed <- 1

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gene_set_dir, recursive = TRUE, showWarnings = FALSE)

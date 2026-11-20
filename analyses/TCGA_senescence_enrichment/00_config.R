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

# Data and results roots (see analyses/paths.R).
source(file.path("analyses", "paths.R"), local = TRUE)

results_dir <- file.path(results_root, "tables", analysis_name)
figures_dir <- file.path(results_root, "figures", analysis_name)
gene_set_dir <- file.path(data_root, analysis_name, "gene_sets")

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
# Gene weights in the enrichment score (fgsea gseaParam). 0 = unweighted:
# only gene ranks count, so rankings dominated by large changes of one sign
# (e.g. kidney tumors, where the largest changes are losses) do not distort
# the null distribution, and NES are comparable between tumors and signatures.
# 1 = classic weighted GSEA.
gsea_weight <- 0
# Permutations for gene sets that fgseaMultilevel cannot score (unbalanced
# rankings); p-values then have a floor of 1 / gsea_permutations.
gsea_permutations <- 10000
# Parallel workers for fgsea (it uses all cores by default).
gsea_workers <- 4
# Leading-edge genes are written for gene sets with padj below this value.
gsea_leading_edge_padj <- fdr_threshold

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gene_set_dir, recursive = TRUE, showWarnings = FALSE)

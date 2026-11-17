# Shared helpers of the enrichment module.

#### Rankings ####

# Named vector of log2 fold changes, one per gene symbol, without missing values.
make_ranking <- function(symbols, values) {
  ok <- !is.na(symbols) & nzchar(symbols) & is.finite(values)
  ranking <- stats::setNames(values[ok], toupper(symbols[ok]))
  ranking[!duplicated(names(ranking))]
}

cancer_ranking <- function(cancer_project) {
  genes <- read.csv(
    file.path(overlap_table_dir, "cancer_gene_tables", paste0(cancer_project, "_genes.csv")),
    stringsAsFactors = FALSE
  )
  shrunken_file <- file.path(
    tcga_table_dir,
    paste0(cancer_project, "_DGE_Results_Tumor_vs_Normal_", tools::toTitleCase(cancer_design_type), ".csv")
  )
  shrunken <- read.csv(shrunken_file, stringsAsFactors = FALSE)
  make_ranking(genes$gene_symbol, shrunken$log2FoldChange[match(genes$gene_id, shrunken$gene_id)])
}

signature_ranking <- function(signature) {
  genes <- read.csv(
    file.path(overlap_table_dir, signature, paste0(signature, "_signature_genes.csv")),
    stringsAsFactors = FALSE
  )
  make_ranking(genes$gene_symbol, genes$log2FoldChange)
}

# Gene sets of one collection as a named list of symbol vectors.
read_gene_sets <- function(collection) {
  table <- readRDS(file.path(gene_set_dir, paste0(collection, ".rds")))
  lapply(split(table$gene_symbol, table$gene_set), unique)
}

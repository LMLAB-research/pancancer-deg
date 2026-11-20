#### Build Gene-Set Collections (step 01) ####
#
# Downloads the gene sets of each collection in gene_set_collections from
# MSigDB (msigdbr) and saves them as one table per collection (gene set name,
# MSigDB ID, gene symbol), plus a version table. Later steps read only these
# files, so all analyses of a run use the same gene-set versions.

source(file.path("analyses", "TCGA_senescence_enrichment", "00_config.R"))

if (!requireNamespace("msigdbr", quietly = TRUE)) {
  stop("Install msigdbr (see README, R setup).")
}

# msigdbr >= 10 uses collection/subcollection; older versions category/subcategory.
fetch_msigdb <- function(collection, subcollection) {
  arguments <- names(formals(msigdbr::msigdbr))
  sets <- if ("collection" %in% arguments) {
    msigdbr::msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection)
  } else {
    msigdbr::msigdbr(species = "Homo sapiens", category = collection, subcategory = subcollection)
  }
  sets <- as.data.frame(sets)
  data.frame(
    gene_set = sets$gs_name,
    gene_set_id = sets$gs_id,
    gene_symbol = toupper(sets$gene_symbol),
    stringsAsFactors = FALSE
  )
}

versions <- dplyr::bind_rows(lapply(names(gene_set_collections), function(name) {
  spec <- gene_set_collections[[name]]
  message("Fetching ", name)
  sets <- unique(fetch_msigdb(spec$collection, spec$subcollection))
  sets <- sets[!is.na(sets$gene_symbol) & nzchar(sets$gene_symbol), , drop = FALSE]
  saveRDS(sets, file.path(gene_set_dir, paste0(name, ".rds")))
  data.frame(
    collection_name = name,
    msigdb_collection = spec$collection,
    msigdb_subcollection = if (is.null(spec$subcollection)) NA_character_ else spec$subcollection,
    n_gene_sets = length(unique(sets$gene_set)),
    n_genes = length(unique(sets$gene_symbol)),
    msigdbr_version = as.character(utils::packageVersion("msigdbr")),
    retrieved = as.character(Sys.Date()),
    stringsAsFactors = FALSE
  )
}))

write.csv(versions, file.path(results_dir, "01_gene_set_versions.csv"), row.names = FALSE)
print(versions, row.names = FALSE)

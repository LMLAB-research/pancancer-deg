#### Preranked Gene-Set Enrichment Analysis (step 03) ####
#
# GSEA (fgseaMultilevel) on rankings of all tested genes, without DEG cutoffs:
# - cancer: genes of each tumor ranked by shrunken log2 fold change (ashr) of
#   the TCGA design used by the overlap module (cancer_design_type), over the
#   genes of its prepared cancer table (protein-coding, one row per symbol);
# - signature: genes of each signature ranked by log2 fold change, over its
#   prepared gene table (after symbol mapping).
# Unweighted by default (gsea_weight = 0, see 00_config.R). NES > 0: gene set
# concentrated among up-regulated genes; NES < 0: among down-regulated genes. BH (fgsea padj) per ranking and collection. Cancer and
# signature NES are compared per gene set in step 04.

source(file.path("analyses", "TCGA_senescence_enrichment", "00_config.R"))

source(file.path("analyses", "TCGA_senescence_enrichment", "utils.R"))

if (!requireNamespace("fgsea", quietly = TRUE)) {
  stop("Install fgsea (see README, R setup).")
}

#### Rankings ####

rankings <- list()
cancer_projects <- sub(
  "_genes\\.csv$", "",
  sort(list.files(file.path(overlap_table_dir, "cancer_gene_tables"), "_genes\\.csv$"))
)
for (cancer_project in cancer_projects) {
  rankings[[length(rankings) + 1L]] <- list(
    kind = "cancer", cancer_project = cancer_project, signature = NA_character_,
    ranking = cancer_ranking(cancer_project)
  )
}
for (signature in signatures) {
  rankings[[length(rankings) + 1L]] <- list(
    kind = "signature", cancer_project = NA_character_, signature = signature,
    ranking = signature_ranking(signature)
  )
}
message(length(rankings), " rankings")

#### GSEA ####

bpparam <- if (gsea_workers > 1 && .Platform$OS.type == "unix") {
  BiocParallel::MulticoreParam(workers = gsea_workers)
} else {
  BiocParallel::SerialParam()
}

# fgseaMultilevel leaves NES and p-value empty for some gene sets when the
# ranking has many more positive than negative values (or the reverse), e.g.
# in kidney tumors. Those gene sets are rerun with the permutation method
# (fgseaSimple), and BH is recomputed over all gene sets. `method` records
# which method gave each row.
run_gsea <- function(ranking, gene_sets) {
  set.seed(gsea_seed)
  result <- suppressWarnings(fgsea::fgseaMultilevel(
    pathways = gene_sets,
    stats = ranking,
    minSize = min_gene_set_size,
    maxSize = max_gene_set_size,
    gseaParam = gsea_weight,
    BPPARAM = bpparam
  ))
  result <- as.data.frame(result)
  result$method <- "multilevel"
  missing <- is.na(result$pval)
  if (any(missing)) {
    set.seed(gsea_seed)
    rerun <- as.data.frame(suppressWarnings(fgsea::fgseaSimple(
      pathways = gene_sets[result$pathway[missing]],
      stats = ranking,
      nperm = gsea_permutations,
      minSize = min_gene_set_size,
      maxSize = max_gene_set_size,
      gseaParam = gsea_weight,
      BPPARAM = bpparam
    )))
    index <- match(rerun$pathway, result$pathway)
    for (column in c("pval", "ES", "NES")) {
      result[[column]][index] <- rerun[[column]]
    }
    result$leadingEdge[index] <- rerun$leadingEdge
    result$method[index] <- "permutation"
    result$padj <- stats::p.adjust(result$pval, method = "BH")
  }
  result$leading_edge <- ifelse(
    !is.na(result$padj) & result$padj < gsea_leading_edge_padj,
    vapply(result$leadingEdge, paste, character(1), collapse = ";"),
    ""
  )
  result$leadingEdge <- NULL
  names(result)[names(result) == "pathway"] <- "gene_set"
  names(result)[names(result) == "pval"] <- "p_value"
  result[order(result$p_value), , drop = FALSE]
}

statuses <- list()
for (collection in names(gene_set_collections)) {
  message("GSEA: ", collection)
  gene_set_table <- readRDS(file.path(gene_set_dir, paste0(collection, ".rds")))
  gene_sets <- lapply(split(gene_set_table$gene_symbol, gene_set_table$gene_set), unique)

  results <- lapply(rankings, function(r) {
    tryCatch(run_gsea(r$ranking, gene_sets), error = function(e) e)
  })

  info <- function(i) {
    data.frame(
      kind = rankings[[i]]$kind,
      cancer_project = rankings[[i]]$cancer_project,
      signature = rankings[[i]]$signature,
      stringsAsFactors = FALSE
    )
  }
  tables <- dplyr::bind_rows(lapply(seq_along(rankings), function(i) {
    if (inherits(results[[i]], "error") || !nrow(results[[i]])) return(NULL)
    cbind(info(i), collection = collection, results[[i]], row.names = NULL)
  }))
  write.csv(tables, file.path(results_dir, paste0("03_GSEA_", collection, ".csv")), row.names = FALSE)

  statuses[[collection]] <- dplyr::bind_rows(lapply(seq_along(rankings), function(i) {
    failed <- inherits(results[[i]], "error")
    cbind(info(i), data.frame(
      collection = collection,
      n_ranked_genes = length(rankings[[i]]$ranking),
      n_gene_sets_tested = if (failed) NA_integer_ else nrow(results[[i]]),
      n_significant = if (failed) NA_integer_ else sum(results[[i]]$padj < fdr_threshold, na.rm = TRUE),
      n_permutation = if (failed) NA_integer_ else sum(results[[i]]$method == "permutation"),
      status = if (failed) paste("failed:", conditionMessage(results[[i]])) else "completed",
      stringsAsFactors = FALSE
    ))
  }))
}

status <- dplyr::bind_rows(statuses)
write.csv(status, file.path(results_dir, "03_GSEA_status.csv"), row.names = FALSE)
message("GSEA results saved to: ", results_dir)
print(table(status$collection, status$status))

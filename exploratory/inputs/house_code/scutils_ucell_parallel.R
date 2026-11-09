#Fork-based parallel UCell scoring
#
#UCell parallelises through BiocParallel. Its default backends start socket workers, which need a local port, and
#sandboxed sessions block local ports. These helpers split the cells into ncores chunks, run UCell on each chunk inside
#parallel::mclapply (forked processes, no sockets) with BPPARAM = BiocParallel::SerialParam(), and put the results back
#in the original cell order. UCell ranks and scores each cell on its own, so chunking by cell gives the same numbers as
#one call on all cells. The caller's random-number stream is left exactly as it was: no draw is made, and .Random.seed
#is restored after the forks.

#Split column indices into at most n contiguous chunks
split_cells <- function(n_cells, n_chunks) {
  n_chunks <- max(1L, min(as.integer(n_chunks), n_cells))
  unname(split(seq_len(n_cells), rep(seq_len(n_chunks), each = ceiling(n_cells / n_chunks), length.out = n_cells)))
}

#Run fun on each chunk of column indices; forked when ncores > 1, serial otherwise
run_chunks <- function(chunks, fun, ncores) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  if (ncores > 1L && length(chunks) > 1L && .Platform$OS.type == "unix") {
    out <- parallel::mclapply(chunks, fun, mc.cores = ncores, mc.set.seed = FALSE, mc.preschedule = FALSE)
  } else {
    out <- lapply(chunks, fun)
  }
  failed <- vapply(out, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(failed)) {
    msg <- unique(vapply(out[failed], function(x) if (is.null(x)) "worker returned NULL (killed?)" else as.character(x), ""))
    stop("UCell failed in ", sum(failed), " of ", length(out), " chunks: ", paste(msg, collapse = "; "), call. = FALSE)
  }
  out
}

#Score UCell signatures on a gene x cell matrix, a stored-rankings matrix or a Seurat object, in forked chunks
#  x           gene x cell matrix (dense or sparse), or a Seurat object (scores are added to its meta.data)
#  signatures  named list of gene vectors (UCell 'features')
#  ncores      number of forked processes; 1 runs serially in this process
#  precalc.ranks  rankings from StoreRankings_UCell / store_rankings_parallel; if given, x is ignored (as in UCell)
#  ...         passed unchanged to UCell::ScoreSignatures_UCell (maxRank, w_neg, name, chunk.size, missing_genes,
#              ties.method, force.gc, assay for SingleCellExperiment input)
#  assay, slot Seurat only: as UCell::AddModuleScore_UCell (assay = DefaultAssay, slot = "counts")
#Returns a cell x signature matrix in the input cell order (Seurat: the object with the scores added).
score_ucell_parallel <- function(x = NULL, signatures, ncores = 1, precalc.ranks = NULL, ...,
                                 assay = NULL, slot = "counts") {
  dots <- list(...)
  if (any(c("BPPARAM", "ncores", "matrix", "features") %in% names(dots))) {
    stop("BPPARAM, ncores, matrix and features are set by score_ucell_parallel; do not pass them in ...", call. = FALSE)
  }
  score_one <- function(m, ranks) {
    args <- c(list(features = signatures, BPPARAM = BiocParallel::SerialParam(), ncores = 1), dots)
    if (is.null(ranks)) args$matrix <- m else args$precalc.ranks <- ranks
    do.call(UCell::ScoreSignatures_UCell, args)
  }
  score_matrix <- function(m, ranks) {
    src <- if (is.null(ranks)) m else ranks
    cells <- colnames(src)
    chunks <- split_cells(ncol(src), ncores)
    res <- run_chunks(chunks, function(j) {
      if (is.null(ranks)) score_one(m[, j, drop = FALSE], NULL) else score_one(NULL, ranks[, j, drop = FALSE])
    }, ncores)
    out <- do.call(rbind, res)
    if (!is.null(cells)) out <- out[cells, , drop = FALSE]
    out
  }
  if (!is.null(precalc.ranks)) return(score_matrix(NULL, precalc.ranks))
  if (inherits(x, "Seurat")) {
    if (is.null(assay)) assay <- Seurat::DefaultAssay(x)
    layers <- SeuratObject::Layers(x, assay = assay, search = slot)
    if (is.null(layers)) stop(sprintf("Cannot find layer %s in assay %s", slot, assay), call. = FALSE)
    scores <- do.call(rbind, lapply(layers, function(l) score_matrix(Seurat::GetAssayData(x, layer = l, assay = assay), NULL)))
    scores <- scores[colnames(x), , drop = FALSE]
    return(Seurat::AddMetaData(x, as.data.frame(scores)))
  }
  score_matrix(x, NULL)
}

#Store UCell rankings (as UCell::StoreRankings_UCell) in forked chunks; returns the genes x cells ranking matrix
#in the input cell order. ... is passed unchanged (maxRank, chunk.size, ties.method, force.gc, assay).
store_rankings_parallel <- function(x, ncores = 1, ...) {
  dots <- list(...)
  if (any(c("BPPARAM", "ncores", "matrix") %in% names(dots))) {
    stop("BPPARAM, ncores and matrix are set by store_rankings_parallel; do not pass them in ...", call. = FALSE)
  }
  cells <- colnames(x)
  chunks <- split_cells(ncol(x), ncores)
  res <- run_chunks(chunks, function(j) {
    do.call(UCell::StoreRankings_UCell, c(list(matrix = x[, j, drop = FALSE], BPPARAM = BiocParallel::SerialParam(), ncores = 1), dots))
  }, ncores)
  genes <- Reduce(union, lapply(res, rownames))
  res <- lapply(res, function(r) {
    miss <- setdiff(genes, rownames(r))
    if (length(miss)) r <- rbind(r, Matrix::Matrix(0, length(miss), ncol(r), sparse = TRUE, dimnames = list(miss, colnames(r))))
    r[genes, , drop = FALSE]
  })
  out <- do.call(cbind, res)
  if (!is.null(cells)) out <- out[, cells, drop = FALSE]
  out
}

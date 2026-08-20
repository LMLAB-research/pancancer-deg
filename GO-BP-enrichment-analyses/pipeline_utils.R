env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = if (default) "true" else "false")
  tolower(trimws(value)) %in% c("1", "true", "yes", "y")
}

env_number <- function(name, default) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (!is.finite(value)) stop("Environment variable ", name, " must be numeric.")
  value
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      ". Run: Rscript run_pipeline.R install"
    )
  }
}

read_gene_symbols <- function(file) {
  if (!file.exists(file)) stop("Missing input file: ", file)
  x <- utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"gene_symbol" %in% names(x)) {
    stop("Column 'gene_symbol' is missing in: ", file)
  }
  genes <- toupper(trimws(as.character(x$gene_symbol)))
  sort(unique(genes[!is.na(genes) & nzchar(genes)]))
}

parse_ratio <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  vapply(parts, function(value) {
    if (length(value) != 2L) return(NA_real_)
    numerator <- suppressWarnings(as.numeric(value[[1]]))
    denominator <- suppressWarnings(as.numeric(value[[2]]))
    if (!is.finite(numerator) || !is.finite(denominator) || denominator == 0) {
      return(NA_real_)
    }
    numerator / denominator
  }, numeric(1))
}

bind_nonempty <- function(x) {
  keep <- vapply(x, function(value) {
    is.data.frame(value) && nrow(value) > 0
  }, logical(1))
  if (!any(keep)) return(data.frame())
  dplyr::bind_rows(x[keep])
}

write_csv_utf8 <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}

collapse_sorted <- function(x, separator = "; ") {
  x <- sort(unique(trimws(as.character(x))))
  x <- x[!is.na(x) & nzchar(x)]
  paste(x, collapse = separator)
}

first_nonempty <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) x[[1]] else NA_character_
}

adjusted_rand_index <- function(reference, candidate) {
  if (length(reference) != length(candidate)) stop("Cluster vectors differ in length.")
  n <- length(reference)
  if (n < 2L) return(NA_real_)

  choose2 <- function(x) x * (x - 1) / 2
  tab <- table(reference, candidate)
  observed <- sum(choose2(tab))
  row_pairs <- sum(choose2(rowSums(tab)))
  col_pairs <- sum(choose2(colSums(tab)))
  total_pairs <- choose2(n)
  expected <- row_pairs * col_pairs / total_pairs
  maximum <- (row_pairs + col_pairs) / 2
  denominator <- maximum - expected
  if (denominator == 0) return(if (observed == expected) 1 else 0)
  (observed - expected) / denominator
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

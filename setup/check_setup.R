suppressPackageStartupMessages({
  ok_renv <- requireNamespace("renv", quietly = TRUE)
})

if (!ok_renv) {
  stop("Package 'renv' is not available. Run: bash setup/setup.sh")
}

project <- normalizePath(".", winslash = "/", mustWork = TRUE)
project_library <- renv::paths$library(project = project)
if (dir.exists(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}

required_packages <- c(
  "AnnotationDbi",
  "BiocParallel",
  "DESeq2",
  "EnhancedVolcano",
  "SummarizedExperiment",
  "TCGAbiolinks",
  "ashr",
  "dplyr",
  "ggplot2",
  "openxlsx",
  "org.Hs.eg.db",
  "pheatmap",
  "readr",
  "recount3"
)

available <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
missing <- names(available)[!available]

cat("R setup check\n")
cat("R version: ", as.character(getRversion()), "\n", sep = "")
cat("Platform: ", R.version$platform, "\n", sep = "")
cat("renv version: ", as.character(utils::packageVersion("renv")), "\n", sep = "")
cat("renv library: ", project_library, "\n", sep = "")
cat("renv.lock: ", if (file.exists("renv.lock")) "found" else "missing", "\n", sep = "")

if (length(missing) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing, collapse = ", "),
    "\nRun: bash setup/setup.sh"
  )
}

cat("Required packages: found\n")

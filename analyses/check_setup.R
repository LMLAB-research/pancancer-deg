suppressPackageStartupMessages({
  ok_renv <- requireNamespace("renv", quietly = TRUE)
})

if (!ok_renv) {
  stop("Package 'renv' is not available. Run: bash scripts/setup.sh")
}

project <- normalizePath(".", winslash = "/", mustWork = TRUE)
project_library <- renv::paths$library(project = project)
if (dir.exists(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}

required_packages <- c(
  "AnnotationDbi",
  "DESeq2",
  "SummarizedExperiment",
  "ggplot2",
  "org.Hs.eg.db",
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
    "\nRun: bash scripts/setup.sh"
  )
}

cat("Required packages: found\n")

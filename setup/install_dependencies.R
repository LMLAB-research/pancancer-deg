suppressPackageStartupMessages({
  if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv", repos = "https://cloud.r-project.org")
  }
})

project <- normalizePath(".", winslash = "/", mustWork = TRUE)
project_library <- renv::paths$library(project = project)

cran_packages <- c(
  "ashr",
  "dplyr",
  "ggplot2",
  "openxlsx"
)

bioc_packages <- c(
  "AnnotationDbi",
  "BiocParallel",
  "DESeq2",
  "org.Hs.eg.db",
  "SummarizedExperiment",
  "TCGAbiolinks"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  renv::install("BiocManager", prompt = FALSE, project = project)
}

if (dir.exists(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}

repos <- BiocManager::repositories()
repos["CRAN"] <- "https://cloud.r-project.org"
options(repos = repos)

packages <- c(cran_packages, paste0("bioc::", bioc_packages))
renv::install(packages, prompt = FALSE, project = project)

# Update renv.lock only on request (after adding or removing a package), so
# that running setup never overwrites the shared lockfile:
#   Rscript --vanilla setup/install_dependencies.R --snapshot
if ("--snapshot" %in% commandArgs(trailingOnly = TRUE)) {
  renv::snapshot(prompt = FALSE, project = project)
  cat("Installed project R dependencies and updated renv.lock.\n")
} else {
  cat("Installed project R dependencies (renv.lock unchanged).\n")
}

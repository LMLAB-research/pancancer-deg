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
  "openxlsx",
  "pheatmap",
  "readr"
)

bioc_packages <- c(
  "AnnotationDbi",
  "DESeq2",
  "EnhancedVolcano",
  "org.Hs.eg.db",
  "recount3",
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
renv::snapshot(prompt = FALSE, project = project)

cat("Installed and snapshotted project R dependencies.\n")

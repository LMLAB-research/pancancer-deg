file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
frame_files <- unlist(lapply(sys.frames(), function(x) {
  if (is.null(x$ofile)) character() else x$ofile
}))
script_file <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[[1]])
} else if (length(frame_files)) {
  tail(frame_files, 1)
} else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  rstudioapi::getActiveDocumentContext()$path
} else {
  file.path(getwd(), "install_packages.R")
}
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))

local_library <- file.path(script_dir, ".Rlib")
dir.create(local_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_library, .libPaths()))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = local_library)
}

cran_packages <- c("dplyr", "ggplot2", "scales")
bioconductor_packages <- c(
  "AnnotationDbi",
  "clusterProfiler",
  "GOSemSim",
  "org.Hs.eg.db",
  "GO.db",
  "rrvgo"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_cran)) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org", lib = local_library)
}

missing_bioc <- bioconductor_packages[
  !vapply(bioconductor_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, lib = local_library, ask = FALSE, update = FALSE)
}

message("Packages are ready in: ", local_library)

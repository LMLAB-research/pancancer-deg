#### RStudio launcher for the GO enrichment pipeline ####

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
frame_files <- unlist(lapply(sys.frames(), function(x) {
  if (is.null(x$ofile)) character() else x$ofile
}))
script_file <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[[1]])
} else if (length(frame_files)) {
  tail(frame_files, 1)
} else {
  file.path(getwd(), "run_pipeline.R")
}
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))

has_overlap_inputs <- function(path) {
  if (is.null(path) || !nzchar(path) || !dir.exists(path)) return(FALSE)
  cancer_dirs <- list.dirs(path, recursive = FALSE, full.names = TRUE)
  cancer_dirs <- cancer_dirs[grepl("^TCGA-", basename(cancer_dirs))]
  required <- c(
    "common_tested_genes.csv",
    "overlap_UP-UP.csv",
    "overlap_UP-DOWN.csv",
    "overlap_DOWN-UP.csv",
    "overlap_DOWN-DOWN.csv"
  )
  any(vapply(cancer_dirs, function(directory) {
    all(file.exists(file.path(directory, required)))
  }, logical(1)))
}

run_step <- function(number, file, label) {
  message("\n=== Step ", number, ": ", label, " ===")
  source(file.path(script_dir, file), chdir = TRUE)
}

run_go_pipeline <- function(mode = c("full", "plot", "install"),
                            input = NULL,
                            output = NULL) {
  mode <- match.arg(mode)

  if (identical(mode, "install")) {
    source(file.path(script_dir, "install_packages.R"), chdir = TRUE)
    return(invisible(TRUE))
  }

  output_root <- if (is.null(output)) {
    normalizePath(file.path(script_dir, "results"), mustWork = FALSE)
  } else {
    normalizePath(path.expand(output), mustWork = FALSE)
  }

  previous_output <- Sys.getenv("GO_PIPELINE_OUTPUT_ROOT", unset = NA_character_)
  previous_input <- Sys.getenv("GO_PIPELINE_INPUT_ROOT", unset = NA_character_)
  on.exit({
    if (is.na(previous_output)) {
      Sys.unsetenv("GO_PIPELINE_OUTPUT_ROOT")
    } else {
      Sys.setenv(GO_PIPELINE_OUTPUT_ROOT = previous_output)
    }
    if (is.na(previous_input)) {
      Sys.unsetenv("GO_PIPELINE_INPUT_ROOT")
    } else {
      Sys.setenv(GO_PIPELINE_INPUT_ROOT = previous_input)
    }
  }, add = TRUE)

  Sys.setenv(GO_PIPELINE_OUTPUT_ROOT = output_root)

  if (identical(mode, "full")) {
    candidates <- if (!is.null(input)) {
      normalizePath(path.expand(input), mustWork = FALSE)
    } else {
      c(
        Sys.getenv("GO_PIPELINE_INPUT_ROOT", unset = ""),
        file.path(script_dir, "data"),
        file.path(dirname(script_dir), "TCGA_OIS_overlapping_pipeline", "results")
      )
    }
    candidates <- unique(candidates[nzchar(candidates)])
    valid <- candidates[vapply(candidates, has_overlap_inputs, logical(1))]

    if (!length(valid)) {
      stop(
        "Cannot find the previous overlap results. Supply their folder with ",
        "run_go_pipeline('full', input = '/path/to/overlap/results')."
      )
    }

    Sys.setenv(
      GO_PIPELINE_INPUT_ROOT = normalizePath(valid[[1]], mustWork = TRUE)
    )
    run_step(1, "01_GO_enrichment.R", "GO enrichment")
    run_step(2, "02_global_semantic_themes.R", "global semantic themes")
  } else {
    required <- file.path(
      output_root,
      paste0("06_TOP", c(5L, 10L), "_cross_cancer_themes_per_class.csv")
    )
    missing <- required[!file.exists(required)]
    if (length(missing)) {
      stop("Run the full pipeline first. Missing: ", paste(missing, collapse = ", "))
    }
  }

  run_step(3, "03_GO_enrichment_map.R", "TOP5 and TOP10 figures")
  run_step(4, "04_validate_results.R", "validation")

  message(
    "\nPipeline complete.\n",
    "Review: ", file.path(output_root, "TOP5_theme_review.csv"), "\n",
    "Review: ", file.path(output_root, "TOP10_theme_review.csv"), "\n",
    "Figures: ", file.path(output_root, "figures")
  )
  invisible(TRUE)
}

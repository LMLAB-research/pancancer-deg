#### Step 4: reproducibility and integrity checks ####

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
  file.path(getwd(), "04_validate_results.R")
}
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))

local_library <- file.path(script_dir, ".Rlib")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
source(file.path(script_dir, "pipeline_utils.R"))
source(file.path(script_dir, "config.R"))
require_packages("dplyr")

common_files <- c(
  "00_GO_enrichment_run_log.csv",
  "01_GO_enrichment_all_results.csv",
  "02_GO_significant_terms.csv",
  "03_global_semantic_cluster_mapping.csv",
  "04_theme_per_cancer_class.csv",
  "05_cross_cancer_theme_summary.csv",
  "08_semantic_sensitivity_summary.csv",
  "09_semantic_sensitivity_assignments.csv"
)

top_files <- unlist(lapply(top_sizes, function(top_n) {
  c(
    top_table_filename(top_n),
    manual_top_table_filename(top_n),
    basename(theme_review_file(top_n)),
    file.path("figures", paste0("GO_TOP", top_n, "_theme_map_data.csv")),
    file.path("figures", paste0("GO_TOP", top_n, "_theme_map.png"))
  )
}))

required_files <- c(common_files, top_files)
missing_files <- required_files[!file.exists(file.path(output_root, required_files))]
if (length(missing_files)) {
  stop("Missing result files: ", paste(missing_files, collapse = ", "))
}

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop("VALIDATION FAILED: ", message)
}

run_log <- utils::read.csv(file.path(output_root, "00_GO_enrichment_run_log.csv"))
significant <- utils::read.csv(file.path(output_root, "02_GO_significant_terms.csv"))
mapping <- utils::read.csv(file.path(output_root, "03_global_semantic_cluster_mapping.csv"))
per_cancer <- utils::read.csv(file.path(output_root, "04_theme_per_cancer_class.csv"))
summary_table <- utils::read.csv(file.path(output_root, "05_cross_cancer_theme_summary.csv"))
sensitivity <- utils::read.csv(
  file.path(output_root, "08_semantic_sensitivity_summary.csv")
)

n_cancers <- dplyr::n_distinct(run_log$cancer)
expected_analyses <- n_cancers * length(classes)

assert(nrow(run_log) == expected_analyses, "run log is not cancers × four classes")
assert(
  !anyDuplicated(run_log[c("cancer", "overlap_class")]),
  "duplicate cancer/class rows in run log"
)
assert(!any(grepl("^ERROR:", run_log$status)), "at least one enrichment analysis failed")
assert(all(significant$p.adjust < bh_cutoff), "non-significant row in significant table")
assert(!anyDuplicated(mapping$GO_ID), "a GO ID has more than one global theme")
assert(
  setequal(unique(significant$ID), mapping$GO_ID),
  "significant GO union and semantic mapping differ"
)
assert(!anyNA(mapping$theme_ID), "missing theme assignment")
assert(
  !anyDuplicated(per_cancer[c("cancer", "overlap_class", "theme_ID")]),
  "duplicate cancer/class/theme summary rows"
)
assert(all(summary_table$n_cancers <= n_cancers), "theme occurs in too many cancers")
assert(
  any(
    sensitivity$analysis == "primary" &
      abs(sensitivity$adjusted_Rand_index_vs_primary - 1) < 1e-12
  ),
  "primary sensitivity reference is missing"
)

validate_top <- function(top_n) {
  top_table <- utils::read.csv(
    file.path(output_root, top_table_filename(top_n)),
    check.names = FALSE
  )
  named_table <- utils::read.csv(
    file.path(output_root, manual_top_table_filename(top_n)),
    check.names = FALSE
  )
  review <- utils::read.csv(theme_review_file(top_n), check.names = FALSE)
  figure_data <- utils::read.csv(
    file.path(
      output_root,
      "figures",
      paste0("GO_TOP", top_n, "_theme_map_data.csv")
    ),
    check.names = FALSE
  )

  top_counts <- table(factor(top_table$overlap_class, levels = classes))
  available_counts <- table(factor(summary_table$overlap_class, levels = classes))
  assert(
    all(as.integer(top_counts) == pmin(top_n, as.integer(available_counts))),
    paste0("unexpected number of TOP", top_n, " rows per class")
  )
  assert(
    !anyDuplicated(top_table[c("overlap_class", "rank_within_class")]),
    paste0("duplicate TOP", top_n, " ranks")
  )

  expected_top <- summary_table |>
    dplyr::group_by(.data$overlap_class) |>
    dplyr::arrange(
      dplyr::desc(.data$n_cancers),
      dplyr::desc(.data$median_log2_FE),
      dplyr::desc(.data$median_gene_ratio),
      .data$theme_ID,
      .by_group = TRUE
    ) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::ungroup()

  assert(
    setequal(
      paste(top_table$overlap_class, top_table$theme_ID),
      paste(expected_top$overlap_class, expected_top$theme_ID)
    ),
    paste0("TOP", top_n, " table does not match the ranking rule")
  )
  assert(
    setequal(
      paste(top_table$overlap_class, top_table$theme_ID),
      paste(named_table$overlap_class, named_table$theme_ID)
    ),
    paste0("manual names changed TOP", top_n, " membership")
  )
  assert(
    all(c("manual_theme_name", "manual_note") %in% names(review)),
    paste0("TOP", top_n, " review table is missing manual columns")
  )
  assert(
    !anyDuplicated(review[c("overlap_class", "theme_ID")]),
    paste0("TOP", top_n, " review table contains duplicate keys")
  )
  assert(
    !anyDuplicated(figure_data[c("row_ID", "occurrence_class")]),
    paste0("TOP", top_n, " figure repeats a term within one class")
  )

  display_rows <- unique(figure_data[c("display_key", "row_ID")])
  assert(
    !anyDuplicated(display_rows$display_key),
    paste0("TOP", top_n, " assigns one display name to multiple rows")
  )

  data.frame(
    top_n = top_n,
    rows = nrow(top_table),
    unique_figure_labels = dplyr::n_distinct(figure_data$row_ID)
  )
}

top_validation <- dplyr::bind_rows(lapply(top_sizes, validate_top))

report <- c(
  "GO enrichment pipeline validation: PASS",
  paste0("Cancers: ", n_cancers),
  paste0("Cancer/class analyses: ", nrow(run_log)),
  paste0("Significant cancer/class/GO rows: ", nrow(significant)),
  paste0("Unique significant GO IDs: ", dplyr::n_distinct(significant$ID)),
  paste0("Global semantic themes: ", dplyr::n_distinct(mapping$theme_ID)),
  unlist(lapply(seq_len(nrow(top_validation)), function(index) {
    c(
      paste0(
        "TOP", top_validation$top_n[[index]],
        " rows: ", top_validation$rows[[index]]
      ),
      paste0(
        "TOP", top_validation$top_n[[index]],
        " unique figure labels: ", top_validation$unique_figure_labels[[index]]
      )
    )
  })),
  paste0("Validated at: ", format(Sys.time(), tz = "UTC", usetz = TRUE))
)

writeLines(report, file.path(output_root, "10_validation_report.txt"), useBytes = TRUE)
message(paste(report, collapse = "\n"))

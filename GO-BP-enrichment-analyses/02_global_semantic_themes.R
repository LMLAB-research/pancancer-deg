#### Step 2: one global semantic reduction across all cancers and classes ####

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
  file.path(getwd(), "02_global_semantic_themes.R")
}
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))

local_library <- file.path(script_dir, ".Rlib")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))

source(file.path(script_dir, "pipeline_utils.R"))
source(file.path(script_dir, "config.R"))

require_packages(c(
  "dplyr", "GO.db", "GOSemSim", "org.Hs.eg.db", "rrvgo"
))

significant_file <- file.path(output_root, "02_GO_significant_terms.csv")
if (!file.exists(significant_file)) {
  stop("Run 01_GO_enrichment.R first. Missing: ", significant_file)
}

significant <- utils::read.csv(
  significant_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "ID", "Description", "p.adjust", "GeneRatio_numeric", "FoldEnrichment",
  "log2_FoldEnrichment", "cancer", "overlap_class"
)
missing_columns <- setdiff(required_columns, names(significant))
if (length(missing_columns)) {
  stop("Missing columns in significant-term table: ", paste(missing_columns, collapse = ", "))
}
if (!nrow(significant)) stop("No significant GO terms were found at BH < ", bh_cutoff)

significant <- significant |>
  dplyr::filter(
    !is.na(.data$ID), nzchar(.data$ID),
    !is.na(.data$p.adjust), .data$p.adjust < bh_cutoff,
    .data$overlap_class %in% classes
  ) |>
  dplyr::mutate(
    cancer = sub("^TCGA-", "", .data$cancer),
    analysis_ID = paste(.data$cancer, .data$overlap_class, sep = "__")
  )

go_ids <- sort(unique(significant$ID))
if (length(go_ids) < 2L) stop("At least two significant GO IDs are required.")

cache_dir <- file.path(output_root, "semantic_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

similarity_matrix <- function(ids, measure) {
  if (!measure %in% c("Wang", "Rel")) {
    stop("Supported semantic similarity measures are 'Wang' and 'Rel'.")
  }

  cache_measure <- if (identical(measure, "Rel")) "GOSemSim_Rel" else measure
  safe_measure <- gsub("[^A-Za-z0-9]+", "_", cache_measure)
  cache_file <- file.path(cache_dir, paste0("GO_BP_similarity_", safe_measure, ".rds"))

  if (file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    if (identical(rownames(cached), ids) && identical(colnames(cached), ids)) {
      message("Using cached ", measure, " similarity matrix.")
      return(cached)
    }
  }

  message("Calculating ", measure, " similarity for ", length(ids), " GO terms...")
  semantic_data <- GOSemSim::godata(
    annoDb = "org.Hs.eg.db",
    ont = "BP",
    computeIC = !identical(measure, "Wang")
  )
  matrix <- GOSemSim::mgoSim(
    ids,
    ids,
    semData = semantic_data,
    measure = measure,
    combine = NULL
  )
  matrix[!is.finite(matrix)] <- 0
  saveRDS(matrix, cache_file)
  matrix
}

cluster_terms <- function(ids, measure, cutoff) {
  matrix <- similarity_matrix(ids, measure)
  set.seed(20260820)
  reduction <- suppressMessages(rrvgo::reduceSimMatrix(
    matrix,
    scores = stats::setNames(rep(1, length(ids)), ids),
    threshold = cutoff,
    orgdb = "org.Hs.eg.db"
  ))
  assignment <- data.frame(
    id = reduction$go,
    cluster = reduction$cluster,
    stringsAsFactors = FALSE
  )
  assignment <- assignment[match(ids, assignment$id), , drop = FALSE]
  if (!identical(assignment$id, ids) || anyNA(assignment$cluster)) {
    stop("Semantic clustering did not return one cluster for every GO ID.")
  }
  assignment
}

settings <- sensitivity_settings
settings$analysis[[1]] <- "primary"
settings$measure[[1]] <- semantic_measure
settings$cutoff[[1]] <- semantic_cutoff
if (!run_sensitivity) settings <- settings[1, , drop = FALSE]

assignments <- vector("list", nrow(settings))
names(assignments) <- settings$analysis

for (index in seq_len(nrow(settings))) {
  message(
    "Semantic clustering: ", settings$analysis[[index]],
    " (", settings$measure[[index]], ", cutoff ", settings$cutoff[[index]], ")"
  )
  assignments[[index]] <- cluster_terms(
    go_ids,
    settings$measure[[index]],
    settings$cutoff[[index]]
  )
}

primary_assignment <- assignments[[1]]
names(primary_assignment) <- c("ID", "semantic_cluster")

term_statistics <- significant |>
  dplyr::group_by(.data$ID) |>
  dplyr::summarise(
    Description = first_nonempty(.data$Description),
    n_significant_lists = dplyr::n_distinct(.data$analysis_ID),
    n_cancers = dplyr::n_distinct(.data$cancer),
    median_BH = safe_median(.data$p.adjust),
    median_log2_FE = safe_median(.data$log2_FoldEnrichment),
    .groups = "drop"
  )

# The representative is chosen deterministically within each semantic cluster:
# widest list support, best median BH, largest median enrichment, then GO ID.
mapping_candidates <- primary_assignment |>
  dplyr::left_join(term_statistics, by = "ID") |>
  dplyr::arrange(
    .data$semantic_cluster,
    dplyr::desc(.data$n_significant_lists),
    .data$median_BH,
    dplyr::desc(.data$median_log2_FE),
    .data$ID
  )

representatives <- mapping_candidates |>
  dplyr::group_by(.data$semantic_cluster) |>
  dplyr::slice_head(n = 1L) |>
  dplyr::ungroup() |>
  dplyr::transmute(
    semantic_cluster = .data$semantic_cluster,
    theme_ID = .data$ID,
    representative_GO = .data$ID,
    auto_name = .data$Description
  )

cluster_sizes <- primary_assignment |>
  dplyr::count(.data$semantic_cluster, name = "cluster_size")

cluster_mapping <- mapping_candidates |>
  dplyr::left_join(representatives, by = "semantic_cluster") |>
  dplyr::left_join(cluster_sizes, by = "semantic_cluster") |>
  dplyr::transmute(
    theme_ID = .data$theme_ID,
    representative_GO = .data$representative_GO,
    auto_name = .data$auto_name,
    semantic_cluster = .data$semantic_cluster,
    cluster_size = .data$cluster_size,
    GO_ID = .data$ID,
    GO_description = .data$Description,
    GO_n_significant_lists = .data$n_significant_lists,
    GO_n_cancers = .data$n_cancers,
    GO_median_BH = .data$median_BH,
    GO_median_log2_FE = .data$median_log2_FE
  ) |>
  dplyr::arrange(.data$theme_ID, .data$GO_ID)

themed_terms <- significant |>
  dplyr::left_join(
    cluster_mapping |>
      dplyr::select(
        GO_ID, theme_ID, representative_GO, auto_name, semantic_cluster
      ),
    by = c("ID" = "GO_ID")
  )

if (anyNA(themed_terms$theme_ID)) stop("Some significant GO terms were not assigned to a theme.")

theme_per_cancer_class <- themed_terms |>
  dplyr::arrange(
    .data$cancer, .data$overlap_class, .data$theme_ID,
    .data$p.adjust, dplyr::desc(.data$log2_FoldEnrichment), .data$ID
  ) |>
  dplyr::group_by(
    .data$cancer, .data$overlap_class, .data$theme_ID,
    .data$representative_GO, .data$auto_name
  ) |>
  dplyr::summarise(
    leading_GO = dplyr::first(.data$ID),
    leading_description = dplyr::first(.data$Description),
    n_significant_GO_terms = dplyr::n_distinct(.data$ID),
    member_GO_IDs = collapse_sorted(.data$ID),
    member_GO_descriptions = collapse_sorted(.data$Description),
    min_BH = min(.data$p.adjust, na.rm = TRUE),
    median_BH = safe_median(.data$p.adjust),
    median_gene_ratio = safe_median(.data$GeneRatio_numeric),
    mean_gene_ratio = safe_mean(.data$GeneRatio_numeric),
    median_fold_enrichment = safe_median(.data$FoldEnrichment),
    median_log2_FE = safe_median(.data$log2_FoldEnrichment),
    .groups = "drop"
  )

class_members <- themed_terms |>
  dplyr::group_by(.data$overlap_class, .data$theme_ID) |>
  dplyr::summarise(
    n_member_GO_terms = dplyr::n_distinct(.data$ID),
    member_GO_IDs = collapse_sorted(.data$ID),
    member_GO_descriptions = collapse_sorted(.data$Description),
    .groups = "drop"
  )

run_log_file <- file.path(output_root, "00_GO_enrichment_run_log.csv")
if (file.exists(run_log_file)) {
  run_log <- utils::read.csv(run_log_file, stringsAsFactors = FALSE)
  total_cancers <- dplyr::n_distinct(run_log$cancer)
} else {
  total_cancers <- dplyr::n_distinct(significant$cancer)
}

cross_cancer_summary <- theme_per_cancer_class |>
  dplyr::group_by(
    .data$overlap_class, .data$theme_ID,
    .data$representative_GO, .data$auto_name
  ) |>
  dplyr::summarise(
    n_cancers = dplyr::n_distinct(.data$cancer),
    percent_cancers = 100 * .data$n_cancers / total_cancers,
    median_gene_ratio = safe_median(.data$median_gene_ratio),
    mean_gene_ratio = safe_mean(.data$median_gene_ratio),
    median_fold_enrichment = safe_median(.data$median_fold_enrichment),
    median_log2_FE = safe_median(.data$median_log2_FE),
    median_min_BH = safe_median(.data$min_BH),
    cancers = collapse_sorted(.data$cancer, separator = ", "),
    .groups = "drop"
  ) |>
  dplyr::left_join(class_members, by = c("overlap_class", "theme_ID")) |>
  dplyr::mutate(class_position = match(.data$overlap_class, classes)) |>
  dplyr::arrange(
    .data$class_position,
    dplyr::desc(.data$n_cancers),
    dplyr::desc(.data$median_log2_FE),
    dplyr::desc(.data$median_gene_ratio),
    .data$theme_ID
  ) |>
  dplyr::select(-class_position)

select_top_themes <- function(summary_table, top_n) {
  summary_table |>
    dplyr::group_by(.data$overlap_class) |>
    dplyr::arrange(
      dplyr::desc(.data$n_cancers),
      dplyr::desc(.data$median_log2_FE),
      dplyr::desc(.data$median_gene_ratio),
      .data$theme_ID,
      .by_group = TRUE
    ) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::mutate(rank_within_class = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::mutate(class_position = match(.data$overlap_class, classes)) |>
    dplyr::arrange(.data$class_position, .data$rank_within_class) |>
    dplyr::select(-class_position)
}

make_curation_template <- function(top_table) {
  top_table |>
    dplyr::transmute(
      overlap_class = .data$overlap_class,
      rank_within_class = .data$rank_within_class,
      theme_ID = .data$theme_ID,
      automatic_representative_GO = .data$representative_GO,
      automatic_GO_description = .data$auto_name,
      n_member_GO_terms = .data$n_member_GO_terms,
      member_GO_IDs = .data$member_GO_IDs,
      member_GO_descriptions = .data$member_GO_descriptions,
      n_cancers = .data$n_cancers,
      percent_cancers = .data$percent_cancers,
      median_gene_ratio = .data$median_gene_ratio,
      median_log2_FE = .data$median_log2_FE,
      manual_theme_name = "",
      manual_note = ""
    )
}

apply_curation <- function(table, curation_file) {
  if (!file.exists(curation_file)) {
    return(dplyr::mutate(
      table,
      manual_theme_name = NA_character_,
      display_name = .data$auto_name
    ))
  }

  curation <- utils::read.csv(
    curation_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required <- c("overlap_class", "theme_ID", "manual_theme_name")
  missing <- setdiff(required, names(curation))
  if (length(missing)) {
    stop("Missing columns in the theme review CSV: ", paste(missing, collapse = ", "))
  }
  if (anyDuplicated(curation[c("overlap_class", "theme_ID")])) {
    stop("The theme review CSV contains duplicate overlap_class + theme_ID keys.")
  }

  curation <- curation |>
    dplyr::mutate(
      overlap_class = as.character(.data$overlap_class),
      theme_ID = as.character(.data$theme_ID),
      manual_theme_name = trimws(as.character(.data$manual_theme_name))
    ) |>
    dplyr::filter(nzchar(.data$manual_theme_name)) |>
    dplyr::select(overlap_class, theme_ID, manual_theme_name)

  table |>
    dplyr::left_join(curation, by = c("overlap_class", "theme_ID")) |>
    dplyr::mutate(
      display_name = dplyr::if_else(
        !is.na(.data$manual_theme_name) & nzchar(.data$manual_theme_name),
        .data$manual_theme_name,
        .data$auto_name
      )
    )
}

refresh_curation_file <- function(template, curation_file) {
  if (!file.exists(curation_file)) return(template)

  previous <- utils::read.csv(
    curation_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required <- c("overlap_class", "theme_ID", "manual_theme_name", "manual_note")
  missing <- setdiff(required, names(previous))
  if (length(missing)) {
    stop("Missing columns in the theme review CSV: ", paste(missing, collapse = ", "))
  }
  if (anyDuplicated(previous[c("overlap_class", "theme_ID")])) {
    stop("The theme review CSV contains duplicate overlap_class + theme_ID keys.")
  }

  previous_manual <- previous |>
    dplyr::transmute(
      overlap_class = as.character(.data$overlap_class),
      theme_ID = as.character(.data$theme_ID),
      manual_theme_name_previous = trimws(as.character(.data$manual_theme_name)),
      manual_note_previous = trimws(as.character(.data$manual_note))
    )

  template |>
    dplyr::select(-manual_theme_name, -manual_note) |>
    dplyr::left_join(previous_manual, by = c("overlap_class", "theme_ID")) |>
    dplyr::mutate(
      manual_theme_name = dplyr::coalesce(.data$manual_theme_name_previous, ""),
      manual_note = dplyr::coalesce(.data$manual_note_previous, "")
    ) |>
    dplyr::select(-manual_theme_name_previous, -manual_note_previous)
}

write_csv_utf8(cluster_mapping, file.path(output_root, "03_global_semantic_cluster_mapping.csv"))
write_csv_utf8(theme_per_cancer_class, file.path(output_root, "04_theme_per_cancer_class.csv"))
write_csv_utf8(cross_cancer_summary, file.path(output_root, "05_cross_cancer_theme_summary.csv"))

write_top_outputs <- function(top_n) {
  top_table <- select_top_themes(cross_cancer_summary, top_n)
  review_file <- theme_review_file(top_n)
  review_table <- make_curation_template(top_table)
  review_table <- refresh_curation_file(review_table, review_file)

  write_csv_utf8(review_table, review_file)
  named_table <- apply_curation(top_table, review_file)
  write_csv_utf8(
    top_table,
    file.path(output_root, top_table_filename(top_n))
  )
  write_csv_utf8(
    named_table,
    file.path(output_root, manual_top_table_filename(top_n))
  )

  data.frame(top_n = top_n, rows = nrow(top_table))
}

top_output_log <- dplyr::bind_rows(lapply(top_sizes, write_top_outputs))

# Sensitivity validation: compare complete GO-term partitions with adjusted Rand index.
assignment_table <- data.frame(GO_ID = go_ids, stringsAsFactors = FALSE)
for (index in seq_along(assignments)) {
  assignment_table[[names(assignments)[[index]]]] <- assignments[[index]]$cluster
}

primary_labels <- assignment_table[["primary"]]
sensitivity_summary <- lapply(seq_len(nrow(settings)), function(index) {
  labels <- assignment_table[[settings$analysis[[index]]]]
  sizes <- as.integer(table(labels))
  data.frame(
    analysis = settings$analysis[[index]],
    measure = settings$measure[[index]],
    cutoff = settings$cutoff[[index]],
    n_GO_terms = length(labels),
    n_clusters = length(unique(labels)),
    median_cluster_size = stats::median(sizes),
    maximum_cluster_size = max(sizes),
    singleton_clusters = sum(sizes == 1L),
    percent_GO_terms_in_singletons = 100 * sum(sizes[sizes == 1L]) / length(labels),
    adjusted_Rand_index_vs_primary = adjusted_rand_index(primary_labels, labels),
    stringsAsFactors = FALSE
  )
}) |>
  dplyr::bind_rows()

write_csv_utf8(
  sensitivity_summary,
  file.path(output_root, "08_semantic_sensitivity_summary.csv")
)
write_csv_utf8(
  assignment_table,
  file.path(output_root, "09_semantic_sensitivity_assignments.csv")
)
capture.output(sessionInfo(), file = file.path(output_root, "00_sessionInfo_step2.txt"))

message(
  "Step 2 complete: ", length(go_ids), " unique GO terms -> ",
  dplyr::n_distinct(cluster_mapping$theme_ID), " global themes.\n",
  paste0("TOP", top_output_log$top_n, " rows: ", top_output_log$rows, collapse = "; ")
)

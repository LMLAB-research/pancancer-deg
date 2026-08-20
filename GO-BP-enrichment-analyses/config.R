# Portable defaults. Personal paths are supplied with environment variables,
# so they are never committed to GitHub.

input_root <- Sys.getenv(
  "GO_PIPELINE_INPUT_ROOT",
  unset = file.path(script_dir, "data")
)

output_root <- Sys.getenv(
  "GO_PIPELINE_OUTPUT_ROOT",
  unset = file.path(script_dir, "results")
)

classes <- c("UP-UP", "UP-DOWN", "DOWN-UP", "DOWN-DOWN")

class_labels <- c(
  "UP-UP" = "Up in cancer, up in OIS",
  "UP-DOWN" = "Up in cancer, down in OIS",
  "DOWN-UP" = "Down in cancer, up in OIS",
  "DOWN-DOWN" = "Down in cancer, down in OIS"
)

bh_cutoff <- env_number("GO_BH_CUTOFF", 0.05)
minimum_input_genes <- as.integer(env_number("GO_MINIMUM_INPUT_GENES", 5))
minimum_gene_set_size <- as.integer(env_number("GO_MINIMUM_GENE_SET_SIZE", 10))
maximum_gene_set_size <- as.integer(env_number("GO_MAXIMUM_GENE_SET_SIZE", 500))
top_sizes <- c(5L, 10L)

top_table_filename <- function(top_n) {
  paste0("06_TOP", top_n, "_cross_cancer_themes_per_class.csv")
}

manual_top_table_filename <- function(top_n) {
  paste0("07_TOP", top_n, "_cross_cancer_themes_with_manual_names.csv")
}

theme_review_file <- function(top_n) {
  file.path(output_root, paste0("TOP", top_n, "_theme_review.csv"))
}

# Primary global semantic-reduction settings.
semantic_measure <- Sys.getenv("GO_SEMANTIC_MEASURE", unset = "Wang")
semantic_cutoff <- env_number("GO_SEMANTIC_CUTOFF", 0.70)

# Sensitivity analyses compare the primary solution with lower/higher cutoffs
# and an alternative semantic similarity measure.
run_sensitivity <- env_flag("GO_RUN_SENSITIVITY", TRUE)
sensitivity_settings <- data.frame(
  analysis = c("primary", "wang_threshold_0.60", "wang_threshold_0.80", "rel_threshold_0.70"),
  measure = c(semantic_measure, "Wang", "Wang", "Rel"),
  cutoff = c(semantic_cutoff, 0.60, 0.80, 0.70),
  stringsAsFactors = FALSE
)

#### Compare Designs And Sex Models (step 06) ####
#
# 1. Design comparison: for each cancer and sex model, compares the raw DESeq2
# results of every pair of designs among naive, adjusted and paired, over their
# common genes:
# - Spearman correlation of log2 fold changes;
# - DEG counts (the deg column of step 05) and their agreement: DEGs in both
#   with the same direction, with opposite directions, and the Jaccard index
#   |D_A and D_B| / |D_A or D_B| (same direction only in the numerator).
#
# 2. Sex-model summary: for each cancer and version (naive, adjusted, paired),
# sample sizes and DEG counts of the all-samples, female, male and interaction
# models, and how the female and male tumor effects compare (Spearman
# correlation of log2 fold changes; DEGs in both sexes with the same or
# opposite direction). A gene that is a DEG in one sex but not the other is not
# evidence of a sex difference (power differs between sexes); only the
# interaction model tests it.

source(file.path("analyses", "TCGA_deg", "00_config.R"))

source(file.path(analysis_dir, "utils.R"))

adjustments <- c("naive", "adjusted", "paired")
design_pairs <- utils::combn(adjustments, 2, simplify = FALSE)

make_design_type <- function(adjustment, sex_model) {
  if (sex_model == "all") adjustment else paste(adjustment, sex_model, sep = "_")
}

# DEG direction per gene: 1 up, -1 down, 0 not a DEG, from the `deg` column
# written by step 05 (so that DEG calls match the test set in 00_config.R).
deg_direction <- function(results_df) {
  if (!"deg" %in% colnames(results_df)) {
    stop("Result table has no deg column; rerun step 05")
  }
  c(up = 1, down = -1, not_significant = 0)[results_df$deg]
}

compare_two_designs <- function(project_id, sex_model, design_a, design_b) {
  file_a <- tcga_deseq_raw_results_file(project_id, design_a)
  file_b <- tcga_deseq_raw_results_file(project_id, design_b)
  if (!file.exists(file_a) || !file.exists(file_b)) {
    return(NULL)
  }

  a <- read.csv(file_a, stringsAsFactors = FALSE)
  b <- read.csv(file_b, stringsAsFactors = FALSE)
  common <- intersect(a$gene_id, b$gene_id)
  a <- a[match(common, a$gene_id), ]
  b <- b[match(common, b$gene_id), ]

  dir_a <- deg_direction(a)
  dir_b <- deg_direction(b)
  n_same <- sum(dir_a != 0 & dir_a == dir_b)
  n_union <- sum(dir_a != 0 | dir_b != 0)

  data.frame(
    project_id = project_id,
    sex_model = sex_model,
    design_a = design_a,
    design_b = design_b,
    n_genes_common = length(common),
    spearman_lfc = stats::cor(
      a$log2FoldChange, b$log2FoldChange,
      method = "spearman", use = "complete.obs"
    ),
    n_deg_a = sum(dir_a != 0),
    n_deg_b = sum(dir_b != 0),
    n_deg_both_same_direction = n_same,
    n_deg_both_opposite_direction = sum(dir_a != 0 & dir_b != 0 & dir_a != dir_b),
    jaccard_deg = if (n_union > 0) n_same / n_union else NA_real_,
    stringsAsFactors = FALSE
  )
}

#### Compare All Projects ####

design_plan <- readRDS(tcga_design_plan_rds)
projects <- limit_projects(unique(design_plan$project_id))

comparison <- dplyr::bind_rows(lapply(projects, function(project_id) {
  dplyr::bind_rows(lapply(sex_models, function(sex_model) {
    dplyr::bind_rows(lapply(design_pairs, function(pair) {
      compare_two_designs(
        project_id,
        sex_model,
        make_design_type(pair[1], sex_model),
        make_design_type(pair[2], sex_model)
      )
    }))
  }))
}))

write.csv(comparison, tcga_design_comparison_file, row.names = FALSE)
message("Design comparison saved to: ", tcga_design_comparison_file)

#### Sex-Model Summary ####

read_result_table <- function(project_id, design_type) {
  file <- tcga_deseq_raw_results_file(project_id, design_type)
  if (!file.exists(file)) {
    return(NULL)
  }
  read.csv(file, stringsAsFactors = FALSE)
}

design_sample_size <- function(project_id, design_type, column) {
  row <- design_plan$project_id == project_id &
    design_plan$design_type == design_type &
    design_plan$eligible
  if (!any(row)) NA_integer_ else design_plan[[column]][which(row)[1]]
}

summarise_sex_models <- function(project_id, adjustment) {
  types <- vapply(sex_models, function(m) make_design_type(adjustment, m), character(1))
  tables <- lapply(types, function(type) read_result_table(project_id, type))
  names(tables) <- sex_models
  if (all(vapply(tables, is.null, logical(1)))) {
    return(NULL)
  }

  n_deg <- vapply(tables, function(x) {
    if (is.null(x)) NA_integer_ else sum(deg_direction(x) != 0)
  }, integer(1))

  # Female vs male tumor effects, over genes tested in both.
  female <- tables$female
  male <- tables$male
  sex_comparison <- if (is.null(female) || is.null(male)) {
    list(spearman = NA_real_, same = NA_integer_, opposite = NA_integer_)
  } else {
    common <- intersect(female$gene_id, male$gene_id)
    f <- female[match(common, female$gene_id), ]
    m <- male[match(common, male$gene_id), ]
    dir_f <- deg_direction(f)
    dir_m <- deg_direction(m)
    list(
      spearman = stats::cor(
        f$log2FoldChange, m$log2FoldChange,
        method = "spearman", use = "complete.obs"
      ),
      same = sum(dir_f != 0 & dir_f == dir_m),
      opposite = sum(dir_f != 0 & dir_m != 0 & dir_f != dir_m)
    )
  }

  data.frame(
    project_id = project_id,
    version = adjustment,
    n_normal_female = design_sample_size(project_id, types[["female"]], "n_normal"),
    n_tumor_female = design_sample_size(project_id, types[["female"]], "n_tumor"),
    n_normal_male = design_sample_size(project_id, types[["male"]], "n_normal"),
    n_tumor_male = design_sample_size(project_id, types[["male"]], "n_tumor"),
    n_deg_all = n_deg[["all"]],
    n_deg_female = n_deg[["female"]],
    n_deg_male = n_deg[["male"]],
    n_deg_sex_interaction = n_deg[["sex_interaction"]],
    spearman_lfc_female_male = sex_comparison$spearman,
    n_deg_female_male_same_direction = sex_comparison$same,
    n_deg_female_male_opposite_direction = sex_comparison$opposite,
    stringsAsFactors = FALSE
  )
}

sex_model_summary <- dplyr::bind_rows(lapply(projects, function(project_id) {
  dplyr::bind_rows(lapply(adjustments, function(adjustment) {
    summarise_sex_models(project_id, adjustment)
  }))
}))

write.csv(sex_model_summary, tcga_sex_model_summary_file, row.names = FALSE)
message("Sex-model summary saved to: ", tcga_sex_model_summary_file)

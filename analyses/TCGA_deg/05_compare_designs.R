#### Compare Naive, Adjusted And Paired Designs ####
#
# For each cancer and sex model, compares the raw DESeq2 results of every pair
# of designs among naive, adjusted and paired, over their common genes:
# - Spearman correlation of log2 fold changes;
# - DEG counts (padj < deseq_alpha and |log2FC| > deseq_lfc_threshold) and
#   their agreement: DEGs in both with the same direction, with opposite
#   directions, and the Jaccard index |D_A and D_B| / |D_A or D_B| (same
#   direction only in the numerator).

source(file.path("analyses", "TCGA_deg", "00_config.R"))

source(file.path(analysis_dir, "utils.R"))

adjustments <- c("naive", "adjusted", "paired")
design_pairs <- utils::combn(adjustments, 2, simplify = FALSE)

make_design_type <- function(adjustment, sex_model) {
  if (sex_model == "all") adjustment else paste(adjustment, sex_model, sep = "_")
}

# DEG direction per gene: 1 up, -1 down, 0 not a DEG.
deg_direction <- function(results_df) {
  significant <- !is.na(results_df$padj) &
    results_df$padj < deseq_alpha &
    abs(results_df$log2FoldChange) > deseq_lfc_threshold
  ifelse(significant, sign(results_df$log2FoldChange), 0)
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

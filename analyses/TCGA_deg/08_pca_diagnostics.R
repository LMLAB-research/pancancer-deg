#### PCA Diagnostics (step 08) ####
#
# What drives the main axes of variation in each cancer? For each cancer, PCA
# on variance-stabilized counts (naive fit, all samples; pca_n_top_genes most
# variable genes), first pca_n_pcs PCs. Each PC is related to sample variables
# by the adjusted R^2 of PC ~ variable (adjusted, so that variables with many
# levels such as plate are not favoured):
# - condition (tumor/normal);
# - clinical: sex, race, age;
# - technical: log10 library size, plate and tissue source site (from the
#   TCGA barcode);
# - paper: annotations from the TCGA marker papers (the `paper_*` columns
#   added by TCGAbiolinks: mostly molecular subtypes and histologies, but also
#   e.g. tissue source site), among tumors only.
# The same analysis is repeated on paired samples only (one normal and one
# tumor per patient), where tumor and normal contribute equally to the total
# variance.

source(file.path("analyses", "TCGA_deg", "00_config.R"))

source(file.path(analysis_dir, "utils.R"))

#### Helpers ####

# Adjusted R^2 of PC ~ variable; NA if too few samples or levels.
adjusted_r2 <- function(pc, x) {
  if (is.character(x) || is.factor(x)) {
    x <- as.character(x)
    x[x %in% metadata_missing_value_labels | x == ""] <- NA
  }
  ok <- !is.na(x) & is.finite(pc)
  if (sum(ok) < 20) {
    return(NA_real_)
  }
  x <- x[ok]
  if (!is.numeric(x)) {
    x <- factor(x)
    if (nlevels(x) < 2 || nlevels(x) >= sum(ok) / 2) {
      return(NA_real_)
    }
  } else if (length(unique(x)) < 3) {
    return(NA_real_)
  }
  summary(stats::lm(pc[ok] ~ x))$adj.r.squared
}

run_pca <- function(vst_matrix) {
  top <- head(
    order(matrixStats::rowVars(vst_matrix), decreasing = TRUE),
    pca_n_top_genes
  )
  pca <- stats::prcomp(t(vst_matrix[top, , drop = FALSE]))
  n_pcs <- min(pca_n_pcs, ncol(pca$x))
  list(
    scores = pca$x[, seq_len(n_pcs), drop = FALSE],
    var_explained = (pca$sdev^2 / sum(pca$sdev^2))[seq_len(n_pcs)]
  )
}

# Sample variables, grouped. Paper annotations are kept if they have 2-15
# levels among tumors.
sample_variables <- function(cd, counts) {
  barcode <- as.character(cd$barcode)
  variables <- list(
    condition = list(group = "condition", values = as.character(cd$condition)),
    sex = list(group = "clinical", values = as.character(cd$sex)),
    race = list(group = "clinical", values = as.character(cd$race)),
    age = list(group = "clinical", values = suppressWarnings(as.numeric(cd$age_at_diagnosis_numeric))),
    library_size = list(group = "technical", values = log10(colSums(counts))),
    plate = list(group = "technical", values = substr(barcode, 22, 25)),
    tissue_source_site = list(group = "technical", values = substr(barcode, 6, 7))
  )
  tumor <- cd$condition == "tumor"
  for (column in grep("^paper_", colnames(cd), value = TRUE)) {
    x <- as.character(cd[[column]])
    x[!tumor] <- NA
    n_levels <- length(unique(x[!is.na(x) & x != ""]))
    if (n_levels >= 2 && n_levels <= 15) {
      variables[[column]] <- list(group = "paper", values = x)
    }
  }
  variables
}

associate <- function(project_id, analysis, pca, variables) {
  dplyr::bind_rows(lapply(seq_along(pca$var_explained), function(i) {
    dplyr::bind_rows(lapply(names(variables), function(name) {
      data.frame(
        project_id = project_id,
        analysis = analysis,
        pc = i,
        var_explained = pca$var_explained[i],
        variable = name,
        variable_group = variables[[name]]$group,
        adj_r2 = adjusted_r2(pca$scores[, i], variables[[name]]$values),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

plot_condition_pcas <- function(project_id, all_pca, all_cd, paired_pca, paired_cd) {
  panel <- function(pca, cd, a, b, label) {
    if (is.null(pca) || ncol(pca$scores) < b) {
      return(NULL)
    }
    data.frame(
      x = pca$scores[, a], y = pca$scores[, b],
      condition = as.character(cd$condition),
      panel = sprintf(
        "%s: PC%d (%.0f%%) vs PC%d (%.0f%%)",
        label, a, 100 * pca$var_explained[a], b, 100 * pca$var_explained[b]
      )
    )
  }
  d <- dplyr::bind_rows(
    panel(all_pca, all_cd, 1, 2, "All samples"),
    panel(all_pca, all_cd, 3, 4, "All samples"),
    panel(paired_pca, paired_cd, 1, 2, "Paired samples")
  )
  d$panel <- factor(d$panel, levels = unique(d$panel))
  figure <- ggplot(d, aes(x, y, color = condition)) +
    geom_point(alpha = 0.6, size = 1) +
    facet_wrap(~panel, scales = "free", nrow = 1) +
    scale_color_manual(values = c(normal = "#3266A7", tumor = "#A73232")) +
    theme_minimal(base_size = 11) +
    labs(title = paste(project_id, "PCA by condition"), x = NULL, y = NULL, color = NULL)
  ggsave(
    tcga_pca_condition_figure_file(project_id), figure,
    width = 4 * nlevels(d$panel), height = 4, dpi = plot_dpi, bg = "white"
  )
}

plot_paper_pca <- function(project_id, all_pca, variables, best_paper) {
  if (is.na(best_paper)) {
    return(invisible(NULL))
  }
  annotation <- variables[[best_paper]]$values
  keep <- !is.na(annotation) & annotation != ""
  d <- data.frame(
    x = all_pca$scores[keep, 1], y = all_pca$scores[keep, 2],
    annotation = annotation[keep]
  )
  figure <- ggplot(d, aes(x, y, color = annotation)) +
    geom_point(alpha = 0.7, size = 1.2) +
    theme_minimal(base_size = 11) +
    labs(
      title = paste(project_id, "tumors: PCA by", sub("^paper_", "", best_paper)),
      x = sprintf("PC1 (%.0f%%)", 100 * all_pca$var_explained[1]),
      y = sprintf("PC2 (%.0f%%)", 100 * all_pca$var_explained[2]),
      color = NULL
    )
  ggsave(
    tcga_pca_paper_figure_file(project_id), figure,
    width = 7, height = 5, dpi = plot_dpi, bg = "white"
  )
}

#### Per-Project Analysis ####

pca_project <- function(project_id) {
  message("PCA diagnostics for ", project_id)
  dds <- readRDS(tcga_project_fit_file(project_id, "naive"))
  vst_matrix <- assay(vst(dds, blind = FALSE))
  cd <- as.data.frame(colData(dds))

  all_pca <- run_pca(vst_matrix)
  all_variables <- sample_variables(cd, counts(dds))
  associations <- associate(project_id, "all_samples", all_pca, all_variables)

  # Paired samples only: patients with both a normal and a tumor sample.
  patient <- as.character(cd$patient_barcode)
  pairs <- intersect(patient[cd$condition == "normal"], patient[cd$condition == "tumor"])
  paired_pca <- NULL
  paired_cd <- NULL
  if (length(pairs) >= min_normal_samples) {
    keep <- patient %in% pairs
    paired_cd <- cd[keep, , drop = FALSE]
    paired_pca <- run_pca(vst_matrix[, keep, drop = FALSE])
    paired_variables <- sample_variables(paired_cd, counts(dds)[, keep, drop = FALSE])
    associations <- dplyr::bind_rows(
      associations,
      associate(project_id, "paired_samples", paired_pca, paired_variables)
    )
  }

  paper_pc1 <- associations[
    associations$analysis == "all_samples" & associations$pc == 1 &
      associations$variable_group == "paper" & !is.na(associations$adj_r2),
  ]
  best_paper <- if (nrow(paper_pc1)) {
    paper_pc1$variable[which.max(paper_pc1$adj_r2)]
  } else {
    NA_character_
  }

  plot_condition_pcas(project_id, all_pca, cd, paired_pca, paired_cd)
  plot_paper_pca(project_id, all_pca, all_variables, best_paper)

  list(
    variance = data.frame(
      project_id = project_id,
      pc = seq_along(all_pca$var_explained),
      var_explained = all_pca$var_explained
    ),
    associations = associations,
    n_samples = ncol(dds),
    n_normal = sum(cd$condition == "normal"),
    n_pairs = length(pairs)
  )
}

# One row per cancer: where tumor/normal sits, and what drives PC1.
summarise_project <- function(result) {
  a <- result$associations
  all_samples <- a[a$analysis == "all_samples" & !is.na(a$adj_r2), ]
  condition <- all_samples[all_samples$variable == "condition", ]
  pc1 <- all_samples[all_samples$pc == 1 & all_samples$variable != "condition", ]
  pc1_paper <- pc1[pc1$variable_group == "paper", ]
  paired_condition <- a[
    a$analysis == "paired_samples" & a$pc == 1 & a$variable == "condition",
  ]
  top_of <- function(d) if (nrow(d)) d[which.max(d$adj_r2), ] else NULL
  pc1_top <- top_of(pc1)
  pc1_top_paper <- top_of(pc1_paper)
  best_condition <- top_of(condition)

  data.frame(
    project_id = result$variance$project_id[1],
    n_samples = result$n_samples,
    n_normal = result$n_normal,
    n_pairs = result$n_pairs,
    pc1_var_explained = result$variance$var_explained[1],
    pc1_condition_adj_r2 = if (any(condition$pc == 1)) condition$adj_r2[condition$pc == 1][1] else NA_real_,
    best_condition_pc = if (is.null(best_condition)) NA_integer_ else best_condition$pc,
    best_condition_adj_r2 = if (is.null(best_condition)) NA_real_ else best_condition$adj_r2,
    pc1_top_other_variable = if (is.null(pc1_top)) NA_character_ else pc1_top$variable,
    pc1_top_other_adj_r2 = if (is.null(pc1_top)) NA_real_ else pc1_top$adj_r2,
    pc1_top_paper_annotation = if (is.null(pc1_top_paper)) NA_character_ else pc1_top_paper$variable,
    pc1_top_paper_annotation_adj_r2 = if (is.null(pc1_top_paper)) NA_real_ else pc1_top_paper$adj_r2,
    paired_pc1_condition_adj_r2 = if (nrow(paired_condition)) paired_condition$adj_r2 else NA_real_,
    stringsAsFactors = FALSE
  )
}

#### Run ####

design_plan <- readRDS(tcga_design_plan_rds)
projects <- limit_projects(unique(
  design_plan$project_id[design_plan$design_type == "naive" & design_plan$eligible]
))
projects <- projects[file.exists(vapply(
  projects, function(p) tcga_project_fit_file(p, "naive"), character(1)
))]

results <- lapply(projects, function(project_id) {
  tryCatch(pca_project(project_id), error = function(e) {
    warning(project_id, ": ", conditionMessage(e))
    NULL
  })
})
results <- results[!vapply(results, is.null, logical(1))]

write.csv(
  dplyr::bind_rows(lapply(results, `[[`, "variance")),
  tcga_pca_variance_file, row.names = FALSE
)
write.csv(
  dplyr::bind_rows(lapply(results, `[[`, "associations")),
  tcga_pca_association_file, row.names = FALSE
)
write.csv(
  dplyr::bind_rows(lapply(results, summarise_project)),
  tcga_pca_summary_file, row.names = FALSE
)
message("PCA diagnostics saved to: ", tcga_pca_summary_file)

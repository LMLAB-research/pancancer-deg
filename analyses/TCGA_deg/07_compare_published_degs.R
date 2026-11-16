#### Compare Our DEGs With Published Cancer-DEGs (step 07) ####
#
# Compares our tumor-vs-normal results with the cancer-DEGs of Chatsirisupachai
# et al. 2019 (Aging Cell, PMC6826163, Data S2): 10 TCGA projects, RSEM counts
# from FireBrowse (GRCh37), limma-voom, primary tumor vs solid normal without
# covariates, DEG = BH-adjusted p < 0.01 and |log2FC| > 1. Data S2 lists only
# their DEGs, so all measures are over their DEGs. Genes are matched by Ensembl
# ID (without version).
#
# Per cancer and design (naive = their design; paired for reference):
# - direction_concordance: share of their DEGs whose log2FC has the same sign
#   in our results (threshold-free);
# - spearman_lfc: Spearman correlation of their and our log2FC over their DEGs;
# - recovery_*: share of their up (down) DEGs that are up (down) DEGs in our
#   results, with our default calls (`deg` column of step 05) and with their
#   rule applied to our fits (beta = 0 test, padj < 0.01, |log2FC| > 1);
# - precision_their_rule: share of our protein-coding DEGs (their rule) that
#   are their DEGs in the same direction (genes they filtered out count as
#   misses, so this is a lower bound);
# - jaccard_their_rule: same-direction DEGs in both / DEGs in either
#   (protein-coding, their rule on both sides).

source(file.path("analyses", "TCGA_deg", "00_config.R"))

source(file.path(analysis_dir, "utils.R"))

#### Published DEGs ####

# Download Data S2 from the Europe PMC supplementary files (a zip containing
# the article's supporting-tables zip) and check its checksum.
fetch_published_degs <- function() {
  if (file.exists(published_deg_file)) {
    return(invisible(published_deg_file))
  }
  dir.create(published_deg_dir, recursive = TRUE, showWarnings = FALSE)
  tmp_dir <- tempfile("published_degs_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  message("Downloading published DEGs from ", published_deg_url)
  # Europe PMC builds the zip on request; R's default 60 s timeout is too short.
  old_options <- options(timeout = max(600, getOption("timeout")))
  on.exit(options(old_options), add = TRUE)
  outer_zip <- file.path(tmp_dir, "supplementary.zip")
  utils::download.file(published_deg_url, outer_zip, mode = "wb", quiet = TRUE)
  utils::unzip(outer_zip, exdir = tmp_dir)
  inner_zip <- list.files(tmp_dir, "s002\\.zip$", full.names = TRUE)
  if (length(inner_zip) != 1) stop("Supporting-tables zip not found in download")
  utils::unzip(inner_zip, exdir = tmp_dir)
  data_s2 <- list.files(tmp_dir, "^Data S2\\.xlsx$", recursive = TRUE, full.names = TRUE)
  data_s2 <- data_s2[!grepl("__MACOSX", data_s2)]
  if (length(data_s2) != 1) stop("Data S2.xlsx not found in download")

  if (unname(tools::md5sum(data_s2)) != published_deg_md5) {
    stop("Unexpected checksum for downloaded Data S2; check the source")
  }
  file.copy(data_s2, published_deg_file)
  invisible(published_deg_file)
}

# One table per cancer: Ensembl ID, log2FC and direction, one row per gene.
read_published_degs <- function(file) {
  sheets <- setdiff(openxlsx::getSheetNames(file), "Data_S2")
  tables <- lapply(sheets, function(sheet) {
    x <- openxlsx::read.xlsx(file, sheet = sheet)
    x <- data.frame(
      ensembl_id = as.character(x$ensembl_gene_id),
      published_log2FoldChange = as.numeric(x$logFC),
      stringsAsFactors = FALSE
    )
    # biomaRt mapping repeats some genes; keep one row per Ensembl ID.
    x <- x[!is.na(x$ensembl_id) & !duplicated(x$ensembl_id), , drop = FALSE]
    x$published_direction <- ifelse(x$published_log2FoldChange > 0, "up", "down")
    x
  })
  stats::setNames(tables, paste0("TCGA-", sheets))
}

#### Comparison ####

strip_ensembl_version <- function(gene_id) sub("\\.[0-9]+$", "", gene_id)

share <- function(numerator, denominator) {
  if (denominator > 0) numerator / denominator else NA_real_
}

compare_with_published <- function(project_id, design_type, published) {
  fit_file <- tcga_project_fit_file(project_id, design_type)
  raw_file <- tcga_deseq_raw_results_file(project_id, design_type)
  design_file <- tcga_project_design_file(project_id, design_type)
  if (!all(file.exists(c(fit_file, raw_file, design_file)))) {
    return(NULL)
  }

  dds <- readRDS(fit_file)
  coefficient <- readRDS(design_file)$coefficient

  # Their rule on our fit: test beta = 0, padj < 0.01, |log2FC| > 1.
  res <- as.data.frame(results(dds, name = coefficient, alpha = published_deg_alpha))
  res$ensembl_id <- strip_ensembl_version(rownames(res))
  res$gene_type <- as.character(rowData(dds)$gene_type)
  their_rule_deg <- !is.na(res$padj) & res$padj < published_deg_alpha &
    abs(res$log2FoldChange) > published_deg_lfc_threshold
  res$deg_their_rule <- ifelse(
    their_rule_deg,
    ifelse(res$log2FoldChange > 0, "up", "down"),
    "not_significant"
  )

  # Our default DEG calls (step 05).
  raw <- read.csv(raw_file, stringsAsFactors = FALSE)
  res$deg_default <- raw$deg[match(res$ensembl_id, strip_ensembl_version(raw$gene_id))]
  res$deg_default[is.na(res$deg_default)] <- "not_significant"

  m <- merge(published, res, by = "ensembl_id")
  m <- m[is.finite(m$log2FoldChange), , drop = FALSE]

  recovery <- function(direction, call_column) {
    in_direction <- m$published_direction == direction
    share(sum(in_direction & m[[call_column]] == direction), sum(in_direction))
  }

  # Precision and Jaccard over protein-coding genes (their gene universe).
  ours_pc <- res[res$gene_type %in% "protein_coding" & res$deg_their_rule != "not_significant", ]
  ours_key <- paste(ours_pc$ensembl_id, ours_pc$deg_their_rule)
  theirs_key <- paste(published$ensembl_id, published$published_direction)
  n_both <- length(intersect(ours_key, theirs_key))

  data.frame(
    project_id = project_id,
    design_type = design_type,
    n_published_deg = nrow(published),
    n_published_up = sum(published$published_direction == "up"),
    n_published_down = sum(published$published_direction == "down"),
    n_matched = nrow(m),
    direction_concordance = mean(sign(m$log2FoldChange) == sign(m$published_log2FoldChange)),
    spearman_lfc = stats::cor(m$log2FoldChange, m$published_log2FoldChange, method = "spearman"),
    recovery_up_default = recovery("up", "deg_default"),
    recovery_down_default = recovery("down", "deg_default"),
    recovery_up_their_rule = recovery("up", "deg_their_rule"),
    recovery_down_their_rule = recovery("down", "deg_their_rule"),
    n_our_deg_their_rule_protein_coding = nrow(ours_pc),
    precision_their_rule = share(n_both, nrow(ours_pc)),
    jaccard_their_rule = share(n_both, length(union(ours_key, theirs_key))),
    stringsAsFactors = FALSE
  )
}

#### Run ####

published_degs <- read_published_degs(fetch_published_degs())
projects <- limit_projects(names(published_degs))

comparison <- dplyr::bind_rows(lapply(projects, function(project_id) {
  dplyr::bind_rows(lapply(published_deg_designs, function(design_type) {
    message("Comparing ", design_type, " results for ", project_id)
    compare_with_published(project_id, design_type, published_degs[[project_id]])
  }))
}))

write.csv(comparison, tcga_published_deg_comparison_file, row.names = FALSE)
message("Published-DEG comparison saved to: ", tcga_published_deg_comparison_file)

#### Figure ####

metric_labels <- c(
  direction_concordance = "Direction concordance",
  spearman_lfc = "Spearman of log2FC",
  recovery_up_their_rule = "Recovery of up DEGs\n(their rule)",
  recovery_down_their_rule = "Recovery of down DEGs\n(their rule)",
  precision_their_rule = "Precision\n(their rule)"
)
plot_df <- do.call(rbind, lapply(names(metric_labels), function(metric) {
  data.frame(
    cancer = sub("^TCGA-", "", comparison$project_id),
    design = comparison$design_type,
    metric = metric_labels[[metric]],
    value = comparison[[metric]],
    stringsAsFactors = FALSE
  )
}))
plot_df$metric <- factor(plot_df$metric, levels = unname(metric_labels))

figure <- ggplot(plot_df, aes(x = value, y = cancer, color = design)) +
  geom_point(size = 2.5) +
  facet_wrap(~metric, nrow = 1) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_y_discrete(limits = rev) +
  theme_minimal(base_size = 11) +
  labs(
    title = "Our tumor-vs-normal results vs Chatsirisupachai et al. 2019 (Data S2)",
    x = NULL,
    y = NULL,
    color = "Our design"
  ) +
  theme(legend.position = "bottom", panel.spacing = grid::unit(1, "lines"))

ggsave(
  tcga_published_deg_figure_file, figure,
  width = 12, height = 4.5, dpi = plot_dpi, bg = "white"
)
message("Figure saved to: ", tcga_published_deg_figure_file)

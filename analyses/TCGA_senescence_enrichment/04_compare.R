#### Compare Cancer And Signature Enrichment (step 04) ####
#
# For each collection, cancer x signature and gene set, compares the GSEA
# results of the tumor and the signature (step 03):
# - direction: "same" (significant on both sides, same NES sign), "opposite"
#   (significant on both, opposite signs), "cancer_only", "signature_only";
#   rows significant on neither side are counted but not written;
# - leading-edge Jaccard index (gene sets significant on both sides): do the
#   same genes drive the enrichment on both sides?
# - within-gene-set correlation rho_set (gene sets significant on both sides):
#   Spearman correlation of tumor and signature log2 fold changes over the
#   genes of the gene set, next to the genome-wide correlation rho_genome of
#   the pair. rho_set well above rho_genome: the genes of the set agree more
#   than the genome overall (shared regulation at gene level).
# Per pair and collection, a summary gives the counts and the Spearman
# correlation of NES across all gene sets tested on both sides.

source(file.path("analyses", "TCGA_senescence_enrichment", "00_config.R"))

source(file.path("analyses", "TCGA_senescence_enrichment", "utils.R"))

library(ggplot2)

#### Helpers ####

split_genes <- function(x) {
  if (is.na(x) || !nzchar(x)) character() else strsplit(x, ";", fixed = TRUE)[[1]]
}

jaccard <- function(a, b) {
  a <- split_genes(a)
  b <- split_genes(b)
  if (!length(a) && !length(b)) NA_real_ else length(intersect(a, b)) / length(union(a, b))
}

spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) NA_real_ else stats::cor(x[ok], y[ok], method = "spearman")
}

direction_category <- function(nes_cancer, padj_cancer, nes_signature, padj_signature) {
  sig_cancer <- !is.na(padj_cancer) & padj_cancer < fdr_threshold
  sig_signature <- !is.na(padj_signature) & padj_signature < fdr_threshold
  ifelse(
    sig_cancer & sig_signature,
    ifelse(sign(nes_cancer) == sign(nes_signature), "same", "opposite"),
    ifelse(sig_cancer, "cancer_only", ifelse(sig_signature, "signature_only", "neither"))
  )
}

#### Rankings (for within-gene-set correlations) ####

gsea_status <- read.csv(file.path(results_dir, "03_GSEA_status.csv"), stringsAsFactors = FALSE)
cancer_projects <- sort(unique(stats::na.omit(gsea_status$cancer_project)))
cancer_rankings <- stats::setNames(lapply(cancer_projects, cancer_ranking), cancer_projects)
signature_rankings <- stats::setNames(lapply(signatures, signature_ranking), signatures)

genome_rho <- do.call(rbind, lapply(cancer_projects, function(cp) {
  do.call(rbind, lapply(signatures, function(sg) {
    genes <- intersect(names(cancer_rankings[[cp]]), names(signature_rankings[[sg]]))
    data.frame(
      cancer_project = cp, signature = sg,
      rho_genome = spearman(cancer_rankings[[cp]][genes], signature_rankings[[sg]][genes]),
      stringsAsFactors = FALSE
    )
  }))
}))

#### Compare ####

summaries <- list()
for (collection in names(gene_set_collections)) {
  message("Comparing: ", collection)
  gsea <- read.csv(file.path(results_dir, paste0("03_GSEA_", collection, ".csv")), stringsAsFactors = FALSE)
  gene_sets <- read_gene_sets(collection)

  cancer <- gsea[gsea$kind == "cancer", c("cancer_project", "gene_set", "NES", "padj", "leading_edge")]
  signature <- gsea[gsea$kind == "signature", c("signature", "gene_set", "NES", "padj", "leading_edge")]
  names(cancer)[3:5] <- paste0(names(cancer)[3:5], "_cancer")
  names(signature)[3:5] <- paste0(names(signature)[3:5], "_signature")

  m <- merge(cancer, signature, by = "gene_set")
  m$direction <- direction_category(m$NES_cancer, m$padj_cancer, m$NES_signature, m$padj_signature)

  summary <- do.call(rbind, lapply(split(m, list(m$cancer_project, m$signature), drop = TRUE), function(d) {
    data.frame(
      cancer_project = d$cancer_project[1],
      signature = d$signature[1],
      collection = collection,
      n_gene_sets = nrow(d),
      n_same = sum(d$direction == "same"),
      n_opposite = sum(d$direction == "opposite"),
      n_cancer_only = sum(d$direction == "cancer_only"),
      n_signature_only = sum(d$direction == "signature_only"),
      spearman_nes = spearman(d$NES_cancer, d$NES_signature),
      stringsAsFactors = FALSE
    )
  }))
  summaries[[collection]] <- merge(summary, genome_rho, by = c("cancer_project", "signature"))

  m <- m[m$direction != "neither", , drop = FALSE]
  both <- m$direction %in% c("same", "opposite")
  m$leading_edge_jaccard <- NA_real_
  m$rho_set <- NA_real_
  m$n_genes_set <- NA_integer_
  for (i in which(both)) {
    m$leading_edge_jaccard[i] <- jaccard(m$leading_edge_cancer[i], m$leading_edge_signature[i])
    cr <- cancer_rankings[[m$cancer_project[i]]]
    sr <- signature_rankings[[m$signature[i]]]
    genes <- intersect(intersect(gene_sets[[m$gene_set[i]]], names(cr)), names(sr))
    m$n_genes_set[i] <- length(genes)
    m$rho_set[i] <- spearman(cr[genes], sr[genes])
  }
  m <- merge(m, genome_rho, by = c("cancer_project", "signature"), all.x = TRUE)
  m$collection <- collection
  m <- m[order(m$cancer_project, m$signature, m$direction, -abs(m$NES_cancer)), c(
    "collection", "cancer_project", "signature", "gene_set", "direction",
    "NES_cancer", "padj_cancer", "NES_signature", "padj_signature",
    "leading_edge_jaccard", "n_genes_set", "rho_set", "rho_genome",
    "leading_edge_cancer", "leading_edge_signature"
  )]
  write.csv(m, file.path(results_dir, paste0("04_compare_", collection, ".csv")), row.names = FALSE)
}

summary <- do.call(rbind, summaries)
write.csv(summary, file.path(results_dir, "04_compare_summary.csv"), row.names = FALSE)

#### Figures ####

cancer_labels <- sub("^TCGA-", "", cancer_projects)

# 1. HALLMARK NES across tumors and signatures (red = up, blue = down;
#    dot = BH-significant). Gene sets ordered by clustering.
hallmark <- read.csv(file.path(results_dir, "03_GSEA_HALLMARK.csv"), stringsAsFactors = FALSE)
hallmark$column <- ifelse(hallmark$kind == "cancer", sub("^TCGA-", "", hallmark$cancer_project), hallmark$signature)
hallmark$gene_set <- sub("^HALLMARK_", "", hallmark$gene_set)
nes <- reshape(hallmark[, c("gene_set", "column", "NES")], idvar = "gene_set", timevar = "column", direction = "wide")
rownames(nes) <- nes$gene_set
nes_matrix <- as.matrix(nes[, -1])
nes_matrix[is.na(nes_matrix)] <- 0
row_order <- rownames(nes_matrix)[stats::hclust(stats::dist(nes_matrix))$order]
hallmark$gene_set <- factor(hallmark$gene_set, levels = row_order)
hallmark$column <- factor(hallmark$column, levels = c(cancer_labels, signatures))
nes_limit <- max(abs(hallmark$NES), na.rm = TRUE)

nes_figure <- ggplot(hallmark, aes(x = column, y = gene_set, fill = NES)) +
  geom_tile(color = "white") +
  geom_point(data = hallmark[which(hallmark$padj < fdr_threshold), ], size = 0.6, color = "black") +
  geom_vline(xintercept = length(cancer_labels) + 0.5, linewidth = 0.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-nes_limit, nes_limit)) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(
    title = "HALLMARK GSEA: tumors (left) and signatures (right)",
    subtitle = "NES > 0 (red): up; NES < 0 (blue): down; dot: BH < 0.05",
    x = NULL, y = NULL
  )
ggsave(file.path(figures_dir, "04_HALLMARK_NES_heatmap.png"), nes_figure,
       width = 9, height = 10, dpi = 300, bg = "white")

# 2. Pathway-level concordance: Spearman correlation of NES per cancer x
#    signature, per collection (blue = same direction, red = opposite).
summary$cancer <- factor(sub("^TCGA-", "", summary$cancer_project), levels = rev(cancer_labels))
summary$signature <- factor(summary$signature, levels = signatures)
summary$collection <- factor(summary$collection, levels = names(gene_set_collections))
correlation_figure <- ggplot(summary, aes(x = signature, y = cancer, fill = spearman_nes)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", spearman_nes)), size = 2.4) +
  facet_wrap(~collection, nrow = 1) +
  scale_fill_gradient2(low = "#A73232", mid = "white", high = "#3266A7", limits = c(-1, 1)) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(
    title = "Correlation of tumor and signature NES across gene sets",
    subtitle = "Spearman rho; > 0 (blue): pathways move in the same direction; < 0 (red): opposite",
    x = NULL, y = NULL, fill = "rho"
  )
ggsave(file.path(figures_dir, "04_NES_correlation_heatmap.png"), correlation_figure,
       width = 14, height = 5.5, dpi = 300, bg = "white")

message("Comparison saved to: ", results_dir)

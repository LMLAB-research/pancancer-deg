library(ggplot2)

results_root <- "/Users/danabiruk/Documents/all_TCGA_OIS_overlaps"
class_order <- c("UP-UP", "UP-DOWN", "DOWN-UP", "DOWN-DOWN")

# Find and combine statistics from all cancers
files <- list.files(
  results_root,
  pattern = "overlap_statistics\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

data <- do.call(rbind, lapply(files, function(f) {
  x <- read.csv(f)
  x$cancer_type <- sub("^TCGA-", "", basename(dirname(f)))
  x
}))

# BH correction across ALL cancer & class tests
data$FDR <- p.adjust(data$p_value, method = "BH")

# Prepare heatmap values
data$overlap_class <- factor(data$overlap_class, levels = class_order)
data$log2_FE <- log2(data$fold_enrichment)

data$label <- paste0(
  data$observed_overlap,
  "\np=", signif(data$FDR, 2)
)

cancer_order <- sort(unique(data$cancer_type))
data$cancer_type <- factor(
  data$cancer_type,
  levels = rev(cancer_order)
)

# Plot
plot_can_ois <- ggplot(
  data,
  aes(overlap_class, cancer_type, fill = log2_FE)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 2.5) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "log2 fold\nenrichment"
  ) +
  labs(
    title = "Directional gene overlap between TCGA cancers and OIS",
    x = "Cancer–OIS overlap",
    y = "Cancer type"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

print(plot_can_ois)

ggsave(
  file.path(results_root, "TCGA_OIS_four_class_heatmap.png"),
  plot_can_ois,
  width = 8,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(results_root, "TCGA_OIS_four_class_heatmap.pdf"),
  plot_can_ois,
  width = 8,
  height = 8
)



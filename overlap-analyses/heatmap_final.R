script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1]])
} else {
  sys.frame(1)$ofile
}
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "config.R"))

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required. Install it with: install.packages('ggplot2')")
}
if (!requireNamespace("scales", quietly = TRUE)) {
  stop("Package 'scales' is required. Install it with: install.packages('scales')")
}

library(ggplot2)

files <- list.files(
  results_dir,
  pattern = "overlap_statistics\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
files <- files[grepl("/TCGA-[^/]+/overlap_statistics\\.csv$", files)]
if (!length(files)) stop("No overlap_statistics.csv files found in: ", results_dir)

read_statistics <- function(file) {
  x <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("overlap_class", "odds_ratio", "p_value")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("Missing columns in ", file, ": ", paste(missing, collapse = ", "))
  }
  x$cancer_type <- sub("^TCGA-", "", basename(dirname(file)))
  x
}

data <- do.call(rbind, lapply(files, read_statistics))
rownames(data) <- NULL
data <- data[
  data$cancer_type %in% cancer_order & data$overlap_class %in% class_order,
  ,
  drop = FALSE
]

expected_rows <- length(cancer_order) * length(class_order)
if (nrow(data) != expected_rows) {
  stop("Expected ", expected_rows, " heatmap rows but found ", nrow(data), ".")
}
if (anyDuplicated(data[c("cancer_type", "overlap_class")])) {
  stop("Duplicate cancer/overlap-class rows found.")
}

# Global BH correction across all cancer-by-class tests in the heatmap.
data$BH <- p.adjust(data$p_value, method = "BH")
data$significant <- !is.na(data$BH) & data$BH < 0.05
data$log2_OR <- ifelse(data$odds_ratio > 0, log2(data$odds_ratio), NA_real_)
data$label <- cut(
  data$BH,
  breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf),
  labels = c("****", "***", "**", "*", "N.S."),
  right = FALSE
)

finite_positive <- data$log2_OR[
  data$significant & is.finite(data$log2_OR) & data$log2_OR > 0
]
fill_limit <- if (length(finite_positive)) max(finite_positive) else 1
data$plot_log2_OR <- pmax(data$log2_OR, 0)
data$plot_log2_OR[data$log2_OR == Inf] <- fill_limit
data$text_color <- ifelse(
  data$significant & data$plot_log2_OR >= 0.55 * fill_limit,
  "white",
  "black"
)

data$cancer_type <- factor(data$cancer_type, levels = rev(cancer_order))
data$overlap_class <- factor(data$overlap_class, levels = class_order)

heatmap <- ggplot(data, aes(overlap_class, cancer_type)) +
  geom_tile(fill = "white", color = "grey78", linewidth = 0.7) +
  geom_tile(
    data = data[data$significant & !is.na(data$plot_log2_OR), , drop = FALSE],
    aes(fill = plot_log2_OR),
    color = "grey78",
    linewidth = 0.7
  ) +
  geom_text(
    aes(label = label, color = text_color),
    size = 5.5,
    fontface = "bold"
  ) +
  scale_color_identity() +
  scale_x_discrete(
    labels = class_labels[class_order],
    position = "top",
    expand = expansion(add = 0)
  ) +
  scale_y_discrete(
    labels = cancer_labels[rev(cancer_order)],
    expand = expansion(add = 0)
  ) +
  scale_fill_gradientn(
    colors = c("white", "#B4C2DE", "#8DA2C9", "#3266A7", "#275392", "#1B427D"),
    limits = c(0, fill_limit),
    oob = scales::squish,
    name = "log2 odds\nratio"
  ) +
  labs(
    title = "Overlap analysis between cancer-DEGs and OIS-DEGs",
    x = NULL,
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 17) +
  theme(
    plot.title = element_text(size = 28, face = "bold", margin = margin(b = 20)),
    plot.title.position = "plot",
    axis.text.x = element_text(face = "bold", size = 13, lineheight = 0.9),
    axis.text.y = element_text(face = "bold", size = 12.5, margin = margin(r = 6)),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 13),
    plot.margin = margin(18, 24, 14, 18)
  )

print(heatmap)

output_base <- file.path(results_dir, "TCGA_OIS_heatmap")
png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"

ggsave(
  paste0(output_base, ".png"), heatmap,
  device = png_device, width = 11, height = 13, units = "in", dpi = 600,
  bg = "white"
)
write.csv(data, paste0(output_base, "_data.csv"), row.names = FALSE)

message("Saved heatmap and plot data to: ", results_dir)

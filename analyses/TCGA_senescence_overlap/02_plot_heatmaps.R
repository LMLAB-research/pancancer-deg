source(file.path("analyses", "TCGA_senescence_overlap", "00_config.R"))

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 first.")
library(ggplot2)

status_file <- file.path(results_dir, "batch_run_status.csv")
if (!file.exists(status_file)) stop("Run 01_compute_overlaps.R first.")
status <- read.csv(status_file, stringsAsFactors = FALSE)
if (!all(c("signature", "cancer_project", "status") %in% names(status)) || !nrow(status)) {
  stop("Invalid batch_run_status.csv.")
}
if (anyNA(status$status) || any(status$status != "COMPLETED")) {
  stop("The batch contains failed comparisons. Fix them before plotting.")
}
files <- file.path(results_dir, status$signature, status$cancer_project, "overlap_statistics.csv")
if (any(!file.exists(files))) stop("Missing statistics files. Rerun the analysis.")

read_statistics <- function(file) {

  x <- read.csv( file, stringsAsFactors = FALSE )

  required <- c("overlap_class", "p_value", "odds_ratio")
  if (!all(required %in% names(x))) stop("Missing statistics columns: ", file)
  if (!is.numeric(x$p_value) || anyNA(x$p_value) ||
      any(!is.finite(x$p_value)) || any(x$p_value < 0 | x$p_value > 1)) {
    stop("Invalid p-values: ", file)
  }
  if (!is.numeric(x$odds_ratio) || anyNA(x$odds_ratio) || any(x$odds_ratio < 0)) {
    stop("Invalid odds ratios: ", file)
  }
  if (nrow(x) != length(class_order) || !setequal(x$overlap_class, class_order)) {
    stop("Expected one row for each overlap class: ", file)
  }
  x$cancer_type <- sub( "^TCGA-", "", basename(dirname(file)) )

  x$signature <- basename( dirname(dirname(file)) )

  x
}

# Combine result tables
data <- do.call( rbind, lapply(files, read_statistics) )

unknown_cancers <- setdiff(unique(data$cancer_type), cancer_order)
if (length(unknown_cancers)) {
  stop("Add these cancers to cancer_order and cancer_labels in 00_config.R: ",
       paste(unknown_cancers, collapse = ", "))
}
if (!setequal(unique(data$signature), expected_signatures)) {
  stop("Batch signatures differ from expected_signatures. Check 00_config.R.")
}
if (anyNA(cancer_labels[unique(data$cancer_type)])) stop("Missing cancer labels in 00_config.R.")

plots_dir <- figures_dir

dir.create( plots_dir, recursive = TRUE, showWarnings = FALSE )

#### Heatmap helper ####

# Heatmap of log2 odds ratios with cancers as rows. `d` has one row per cell,
# with columns cancer_type, `x_column`, p_value and odds_ratio. BH is computed
# within each level of `bh_group`; only significant cells are colored
# (blue = log2 OR > 0, red = log2 OR < 0). Saves the PNG to plots_dir and the
# plotted values to results_dir, both named `<file_stem>_heatmap*`.
plot_or_heatmap <- function(
    d,
    x_column,
    x_levels,
    bh_group,
    title,
    legend_name,
    file_stem,
    x_labels = waiver(),
    x_text = element_text(size = 9.2, face = "plain", lineheight = 0.9),
    legend_title_size = 14
) {
  if (!nrow(d)) stop("No results for heatmap: ", file_stem)
  if (anyDuplicated(d[c("cancer_type", x_column)])) stop("Duplicate heatmap rows.")

  d$BH <- ave(d$p_value, bh_group, FUN = function(p) p.adjust(p, method = "BH"))
  d$significant <- !is.na(d$BH) & d$BH < 0.05

  # OR = 0 (no shared genes) gives -Inf; it is drawn at the lower end of the scale.
  d$log2_OR <- log2(d$odds_ratio)

  d$label <- as.character(cut(
    d$BH,
    breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf),
    labels = c("****", "***", "**", "*", "N.S."),
    right = FALSE
  ))
  d$label[is.na(d$odds_ratio)] <- "NA"

  # Symmetric diverging scale around 0.
  finite_OR <- d$log2_OR[d$significant & is.finite(d$log2_OR)]
  legend_max <- ceiling((if (length(finite_OR)) max(abs(finite_OR)) else 1) * 2) / 2
  legend_breaks <- pretty(c(-legend_max, legend_max))
  legend_breaks <- legend_breaks[abs(legend_breaks) <= legend_max]
  d$plot_log2_OR <- pmin(pmax(d$log2_OR, -legend_max), legend_max)
  d$text_color <- ifelse(d$significant & abs(d$plot_log2_OR) >= 0.55 * legend_max, "white", "black")

  plot_cancers <- cancer_order[cancer_order %in% d$cancer_type]
  d$cancer_type <- factor(d$cancer_type, levels = rev(plot_cancers))
  d[[x_column]] <- factor(d[[x_column]], levels = x_levels)

  heatmap <- ggplot(d, aes(x = .data[[x_column]], y = cancer_type)) +
    geom_tile(fill = "white", color = "grey78", linewidth = 0.6) +
    geom_tile(
      data = d[d$significant & !is.na(d$plot_log2_OR), ],
      aes(fill = plot_log2_OR),
      color = "grey78",
      linewidth = 0.6
    ) +
    geom_text(aes(label = label, color = text_color), size = 4) +
    scale_color_identity() +
    scale_x_discrete(labels = x_labels, position = "top", expand = expansion(add = 0)) +
    scale_y_discrete(labels = cancer_labels[rev(plot_cancers)], expand = expansion(add = 0)) +
    scale_fill_gradientn(
      colors = c("#7D1B1B", "#A73232", "#C98D8D", "white", "#8DA2C9", "#3266A7", "#1B427D"),
      values = c(0, 1/6, 2/6, 0.5, 4/6, 5/6, 1),
      limits = c(-legend_max, legend_max),
      breaks = legend_breaks,
      name = legend_name,
      guide = guide_colorbar(
        barheight = grid::unit(5, "cm"),
        barwidth = grid::unit(0.65, "cm"),
        ticks = TRUE,
        draw.llim = FALSE,
        draw.ulim = FALSE
      )
    ) +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(size = 16.5, face = "bold", margin = margin(b = 14)),
      plot.title.position = "plot",
      axis.text.x = x_text,
      axis.text.y = element_text(size = 9.2, face = "plain", margin = margin(r = 5)),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = legend_title_size, face = "bold"),
      legend.text = element_text(size = 11.5),
      plot.margin = margin(10, 12, 10, 10)
    )

  ggsave(
    file.path(plots_dir, paste0(file_stem, "_heatmap.png")),
    heatmap,
    width = 7.2, height = max(3.5, 1.5 + 0.41 * length(plot_cancers)),
    units = "in", dpi = 600, bg = "white"
  )
  write.csv(d, file.path(results_dir, paste0(file_stem, "_heatmap_data.csv")), row.names = FALSE)
}

#### Class-overlap heatmaps (one per signature, cancers x classes) ####

# BH covers all cancer-by-class tests within each signature.
for (signature in expected_signatures) {
  d <- data[data$signature == signature, ]
  plot_or_heatmap(
    d,
    x_column = "overlap_class",
    x_levels = class_order,
    bh_group = d$signature,
    title = paste0("Overlap analysis between cancer-DEGs and ", signature, "-DEGs"),
    legend_name = "log2 odds\nratio",
    file_stem = paste0("TCGA_", signature),
    x_labels = make_class_labels(signature)[class_order]
  )
}

#### Share and direction heatmaps (cancers x signatures) ####

share_direction_files <- file.path(
  results_dir, status$signature, status$cancer_project, "share_direction_statistics.csv"
)
if (any(!file.exists(share_direction_files))) {
  stop("Missing share_direction_statistics.csv files. Rerun 01_compute_overlaps.R.")
}
share_direction <- do.call(rbind, lapply(share_direction_files, function(file) {
  x <- read.csv(file, stringsAsFactors = FALSE)
  if (!setequal(x$test, c("share", "direction"))) stop("Expected share and direction rows: ", file)
  x$cancer_type <- sub("^TCGA-", "", basename(dirname(file)))
  x
}))

# BH covers all cancers within each signature.
for (test_name in c("share", "direction")) {
  d <- share_direction[share_direction$test == test_name, ]
  plot_or_heatmap(
    d,
    x_column = "signature",
    x_levels = expected_signatures,
    bh_group = d$signature,
    title = if (test_name == "share") {
      "Sharing: DE in cancer vs DE in signature"
    } else {
      "Direction among genes DE in both"
    },
    legend_name = if (test_name == "share") {
      "log2 odds ratio\n(> 0: DE genes\nshared)"
    } else {
      "log2 odds ratio\n(> 0 same,\n< 0 opposite)"
    },
    file_stem = paste0("TCGA_", test_name),
    x_text = element_text(size = 11, face = "bold"),
    legend_title_size = 12
  )
}

message("Heatmaps saved to: ", plots_dir)

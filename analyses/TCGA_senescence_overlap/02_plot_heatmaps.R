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

# Build a separate heatmap 
for (signature in expected_signatures) {

  d <- data[
    data$signature == signature,
  ]

  # BH correction covers all cancer-by-class tests within this signature
  if (!nrow(d)) stop("No results for signature: ", signature)
  if (anyDuplicated(d[c("cancer_type", "overlap_class")])) stop("Duplicate heatmap rows.")
  d$BH <- p.adjust( d$p_value, method = "BH" )

  d$significant <- !is.na(d$BH) & d$BH < 0.05

  # OR = 0 (no shared genes) gives -Inf; it is drawn at the lower end of the scale.
  d$log2_OR <- log2(d$odds_ratio)

  d$label <- cut(
    d$BH,
    breaks = c( -Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf ),
    labels = c( "****", "***", "**", "*", "N.S." ),
    right = FALSE
  )

  # Symmetric diverging scale: blue = enrichment (log2 OR > 0), red = depletion.
  finite_OR <- d$log2_OR[
    d$significant &
      is.finite(d$log2_OR)
  ]

  fill_limit <- if (length(finite_OR)) {
    max(abs(finite_OR))
  } else {
    1
  }

  legend_max <- ceiling(fill_limit * 2) / 2

  legend_breaks <- pretty(c(-legend_max, legend_max)); legend_breaks <- legend_breaks[abs(legend_breaks) <= legend_max]

  d$plot_log2_OR <- pmin( pmax( d$log2_OR, -legend_max ), legend_max )

  d$text_color <- ifelse( d$significant & abs(d$plot_log2_OR) >= 0.55 * legend_max, "white", "black" )

  plot_cancers <- cancer_order[cancer_order %in% d$cancer_type]
  d$cancer_type <- factor(d$cancer_type, levels = rev(plot_cancers))

  d$overlap_class <- factor( d$overlap_class, levels = class_order )

  labels <- make_class_labels(signature)

  # Color significant cells by log2 odds ratio and overlay significance labels.
  heatmap <- ggplot( d, aes( x = overlap_class, y = cancer_type ) ) +

    geom_tile( fill = "white", color = "grey78", linewidth = 0.6 ) +

    geom_tile(
      data = d[
        d$significant &
          !is.na(d$plot_log2_OR),
      ],
      aes(fill = plot_log2_OR),
      color = "grey78",
      linewidth = 0.6
    ) +

    geom_text( aes( label = label, color = text_color ), size = 4 ) +

    scale_color_identity() +

    scale_x_discrete( labels = labels[class_order], position = "top", expand = expansion(add = 0) ) +

    scale_y_discrete( labels = cancer_labels[rev(plot_cancers)], expand = expansion(add = 0) ) +

    scale_fill_gradientn(
      colors = c(
        "#7D1B1B",
        "#A73232",
        "#C98D8D",
        "white",
        "#8DA2C9",
        "#3266A7",
        "#1B427D"
      ),
      values = c( 0, 1/6, 2/6, 0.5, 4/6, 5/6, 1 ),
      limits = c( -legend_max, legend_max ),
      breaks = legend_breaks,
      name = "log2 odds\nratio",
      guide = guide_colorbar(
        barheight = grid::unit(5, "cm"),
        barwidth = grid::unit(0.65, "cm"),
        ticks = TRUE,
        draw.llim = FALSE,
        draw.ulim = FALSE
      )
    ) +

    labs(
      title = paste0( "Overlap analysis between cancer-DEGs and ", signature, "-DEGs" ),
      x = NULL,
      y = NULL
    ) +

    theme_classic( base_size = 13 ) +

    theme(
      plot.title = element_text( size = 16.5, face = "bold", margin = margin(b = 14) ),

      plot.title.position = "plot",

      axis.text.x = element_text( size = 9.2, face = "plain", lineheight = 0.9 ),

      axis.text.y = element_text( size = 9.2, face = "plain", margin = margin(r = 5) ),

      axis.ticks = element_blank(),
      axis.line = element_blank(),

      legend.position = "right",

      legend.title = element_text( size = 14, face = "bold" ),

      legend.text = element_text( size = 11.5 ),

      plot.margin = margin( 10, 12, 10, 10 )
    )

  # Export the heatmap and the data 
  output_file <- file.path( plots_dir, paste0( "TCGA_", signature, "_heatmap.png" ) )

  ggsave( output_file, heatmap, width = 7.2, height = max(3.5, 1.5 + 0.41 * length(plot_cancers)), units = "in", dpi = 600, bg = "white" )

  write.csv(
    d,
    file.path( results_dir, paste0( "TCGA_", signature, "_heatmap_data.csv" ) ),
    row.names = FALSE
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

# One heatmap of log2 odds ratios for one test. Stars use BH across cancers
# within each signature; only significant cells are colored.
plot_share_direction_heatmap <- function(test_name, title, legend_name) {
  d <- share_direction[share_direction$test == test_name, ]
  if (anyDuplicated(d[c("cancer_type", "signature")])) stop("Duplicate heatmap rows.")

  d$BH <- ave(d$p_value, d$signature, FUN = function(p) p.adjust(p, method = "BH"))
  d$significant <- !is.na(d$BH) & d$BH < 0.05
  d$log2_OR <- log2(d$odds_ratio)
  d$label <- as.character(cut(
    d$BH,
    breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf),
    labels = c("****", "***", "**", "*", "N.S."),
    right = FALSE
  ))
  d$label[is.na(d$odds_ratio)] <- "NA"

  finite_OR <- d$log2_OR[d$significant & is.finite(d$log2_OR)]
  legend_max <- ceiling((if (length(finite_OR)) max(abs(finite_OR)) else 1) * 2) / 2
  legend_breaks <- pretty(c(-legend_max, legend_max)); legend_breaks <- legend_breaks[abs(legend_breaks) <= legend_max]
  d$plot_log2_OR <- pmin(pmax(d$log2_OR, -legend_max), legend_max)
  d$text_color <- ifelse(d$significant & abs(d$plot_log2_OR) >= 0.55 * legend_max, "white", "black")

  plot_cancers <- cancer_order[cancer_order %in% d$cancer_type]
  d$cancer_type <- factor(d$cancer_type, levels = rev(plot_cancers))
  d$signature <- factor(d$signature, levels = expected_signatures)

  heatmap <- ggplot(d, aes(x = signature, y = cancer_type)) +
    geom_tile(fill = "white", color = "grey78", linewidth = 0.6) +
    geom_tile(
      data = d[d$significant & !is.na(d$plot_log2_OR), ],
      aes(fill = plot_log2_OR),
      color = "grey78",
      linewidth = 0.6
    ) +
    geom_text(aes(label = label, color = text_color), size = 4) +
    scale_color_identity() +
    scale_x_discrete(position = "top", expand = expansion(add = 0)) +
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
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 9.2, margin = margin(r = 5)),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11.5),
      plot.margin = margin(10, 12, 10, 10)
    )

  ggsave(
    file.path(plots_dir, paste0("TCGA_", test_name, "_heatmap.png")),
    heatmap,
    width = 7.2, height = max(3.5, 1.5 + 0.41 * length(plot_cancers)),
    units = "in", dpi = 600, bg = "white"
  )
  write.csv(
    d,
    file.path(results_dir, paste0("TCGA_", test_name, "_heatmap_data.csv")),
    row.names = FALSE
  )
}

plot_share_direction_heatmap(
  "share",
  "Sharing: DE in cancer vs DE in signature",
  "log2 odds ratio\n(> 0: DE genes\nshared)"
)
plot_share_direction_heatmap(
  "direction",
  "Direction among genes DE in both",
  "log2 odds ratio\n(> 0 same,\n< 0 opposite)"
)

message( "Heatmaps saved to: ", plots_dir )

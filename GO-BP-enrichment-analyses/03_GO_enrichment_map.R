#### Step 3: cross-cancer GO theme maps for TOP5 and TOP10 ####

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
frame_files <- unlist(lapply(sys.frames(), function(x) {
  if (is.null(x$ofile)) character() else x$ofile
}))
script_file <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[[1]])
} else if (length(frame_files)) {
  tail(frame_files, 1)
} else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  rstudioapi::getActiveDocumentContext()$path
} else {
  file.path(getwd(), "03_GO_enrichment_map.R")
}
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))

local_library <- file.path(script_dir, ".Rlib")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))

source(file.path(script_dir, "pipeline_utils.R"))
source(file.path(script_dir, "config.R"))
require_packages(c("dplyr", "ggplot2", "scales"))

read_manual_names <- function(file) {
  if (!file.exists(file)) {
    return(data.frame(
      overlap_class = character(),
      theme_ID = character(),
      manual_theme_name = character()
    ))
  }

  review <- utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("overlap_class", "theme_ID", "manual_theme_name")
  missing <- setdiff(required, names(review))
  if (length(missing)) {
    stop("Missing columns in ", basename(file), ": ", paste(missing, collapse = ", "))
  }
  if (anyDuplicated(review[c("overlap_class", "theme_ID")])) {
    stop(basename(file), " contains duplicate overlap_class + theme_ID keys.")
  }

  review |>
    dplyr::transmute(
      overlap_class = as.character(.data$overlap_class),
      theme_ID = as.character(.data$theme_ID),
      manual_theme_name = trimws(as.character(.data$manual_theme_name))
    ) |>
    dplyr::filter(nzchar(.data$manual_theme_name))
}

create_theme_map <- function(top_n) {
  top_file <- file.path(output_root, top_table_filename(top_n))
  if (!file.exists(top_file)) {
    stop("Run 02_global_semantic_themes.R first. Missing: ", top_file)
  }

  top_table <- utils::read.csv(
    top_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  manual_names <- read_manual_names(theme_review_file(top_n))

  selected_rows <- top_table |>
    dplyr::left_join(manual_names, by = c("overlap_class", "theme_ID")) |>
    dplyr::mutate(
      manual_theme_name = dplyr::coalesce(.data$manual_theme_name, ""),
      display_name = dplyr::if_else(
        nzchar(.data$manual_theme_name),
        .data$manual_theme_name,
        .data$auto_name
      ),
      display_key = tolower(gsub(
        "[[:space:]]+", " ", trimws(.data$display_name)
      )),
      selected_class = .data$overlap_class,
      selected_class_position = match(.data$selected_class, classes)
    ) |>
    dplyr::arrange(
      .data$selected_class_position,
      .data$rank_within_class,
      .data$theme_ID
    )

  row_definitions <- selected_rows |>
    dplyr::distinct(.data$display_key, .keep_all = TRUE) |>
    dplyr::transmute(
      display_key = .data$display_key,
      row_ID = sprintf("TERM_%02d", dplyr::row_number()),
      row_label = .data$display_name
    )

  selected_rows <- selected_rows |>
    dplyr::left_join(row_definitions, by = "display_key")

  plot_data <- selected_rows |>
    dplyr::mutate(
      occurrence_class = .data$selected_class,
      occurrence_label = unname(class_labels[.data$selected_class])
    ) |>
    dplyr::distinct(.data$row_ID, .data$occurrence_class, .keep_all = TRUE)

  row_labels <- stats::setNames(
    row_definitions$row_label,
    row_definitions$row_ID
  )
  plot_data$row_ID <- factor(
    plot_data$row_ID,
    levels = rev(row_definitions$row_ID)
  )
  plot_data$occurrence_label <- factor(
    plot_data$occurrence_label,
    levels = unname(class_labels[classes])
  )

  colour_max <- max(plot_data$median_log2_FE, na.rm = TRUE)
  if (!is.finite(colour_max) || colour_max <= 0) colour_max <- 1
  y_text_size <- if (top_n <= 5L) 13 else 11.5

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$occurrence_label, y = .data$row_ID)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = .data$median_gene_ratio,
        fill = .data$median_log2_FE
      ),
      shape = 21,
      colour = "white",
      stroke = 0.4
    ) +
    ggplot2::scale_y_discrete(labels = row_labels, drop = FALSE) +
    ggplot2::scale_size_continuous(
      name = "Median GeneRatio",
      range = c(4, 11),
      labels = scales::label_number(accuracy = 0.001)
    ) +
    ggplot2::scale_fill_gradientn(
      name = "Median log2 FE",
      colours = c("#2C7BB6", "#8E4A9D", "#D7191C"),
      limits = c(0, colour_max),
      oob = scales::squish
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(
        order = 1,
        override.aes = list(fill = "black", colour = "black")
      ),
      fill = ggplot2::guide_colourbar(order = 2)
    ) +
    ggplot2::labs(
      title = paste0(
        "GO enrichment analyses of overlapping genes\n",
        "between cancer-DEGs and OIS-DEGs"
      ),
      x = NULL,
      y = NULL
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 17) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = 28, lineheight = 1.1,
        hjust = 0, margin = ggplot2::margin(b = 20)
      ),
      plot.title.position = "plot",
      axis.text.x = ggplot2::element_text(
        face = "bold", size = 15, angle = 45, hjust = 1, vjust = 1,
        colour = "#202020"
      ),
      axis.text.y = ggplot2::element_text(
        size = y_text_size, colour = "#202020",
        margin = ggplot2::margin(r = 8)
      ),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "#333333", linewidth = 0.7),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(
        colour = "#E3E3E3", linewidth = 0.5
      ),
      panel.grid.major.y = ggplot2::element_line(
        colour = "#D5D5D5", linewidth = 0.55
      ),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold", size = 15),
      legend.text = ggplot2::element_text(size = 13),
      legend.key.height = grid::unit(0.65, "cm"),
      plot.margin = ggplot2::margin(18, 24, 18, 18)
    )

  figure_dir <- file.path(output_root, "figures")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  write_csv_utf8(
    selected_rows |>
      dplyr::select(
        -selected_class_position,
        -display_key,
        -selected_class,
        -row_ID,
        -row_label
      ),
    file.path(output_root, manual_top_table_filename(top_n))
  )
  write_csv_utf8(
    plot_data |>
      dplyr::mutate(
        row_ID = as.character(.data$row_ID),
        occurrence_label = as.character(.data$occurrence_label)
      ),
    file.path(figure_dir, paste0("GO_TOP", top_n, "_theme_map_data.csv"))
  )

  figure_file <- file.path(figure_dir, paste0("GO_TOP", top_n, "_theme_map.png"))
  figure_height <- max(13, 0.50 * nrow(row_definitions) + 5)
  png_device <- if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png
  } else {
    "png"
  }
  ggplot2::ggsave(
    figure_file,
    plot,
    device = png_device,
    width = 11,
    height = figure_height,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )

  data.frame(
    top_n = top_n,
    rows = nrow(selected_rows),
    unique_labels = nrow(row_definitions),
    figure = figure_file,
    stringsAsFactors = FALSE
  )
}

plot_log <- dplyr::bind_rows(lapply(top_sizes, create_theme_map))
capture.output(sessionInfo(), file = file.path(output_root, "00_sessionInfo_step3.txt"))

message(
  "Step 3 complete:\n",
  paste0("TOP", plot_log$top_n, ": ", plot_log$figure, collapse = "\n")
)

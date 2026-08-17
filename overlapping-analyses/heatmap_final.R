library(ggplot2)


#### 1. Settings ####

results_root <- "/Users/danabiruk/Documents/all_TCGA_OIS_overlaps"

cancer_order <- c(
  "HNSC", "ESCA", "LUAD", "LUSC",
  "STAD", "COAD", "READ", "LIHC", "CHOL",
  "KICH", "KIRC", "KIRP", "BLCA",
  "THCA", "BRCA", "PRAD", "UCEC"
)

cancer_labels <- c(
  HNSC = "HNSC (head and neck)",
  ESCA = "ESCA (esophagus)",
  LUAD = "LUAD (lung adenocarcinoma)",
  LUSC = "LUSC (lung squamous)",
  STAD = "STAD (stomach)",
  COAD = "COAD (colon)",
  READ = "READ (rectum)",
  LIHC = "LIHC (liver)",
  CHOL = "CHOL (bile duct)",
  KICH = "KICH (kidney chromophobe)",
  KIRC = "KIRC (kidney clear cell)",
  KIRP = "KIRP (kidney papillary)",
  BLCA = "BLCA (bladder)",
  THCA = "THCA (thyroid)",
  BRCA = "BRCA (breast)",
  PRAD = "PRAD (prostate)",
  UCEC = "UCEC (uterus)"
)

class_order <- c(
  "UP-UP",
  "UP-DOWN",
  "DOWN-UP",
  "DOWN-DOWN"
)

class_labels <- c(
  "UP-UP" = "Up_cancer,\nUp_OIS",
  "UP-DOWN" = "Up_cancer,\nDown_OIS",
  "DOWN-UP" = "Down_cancer,\nUp_OIS",
  "DOWN-DOWN" = "Down_cancer,\nDown_OIS"
)

fill_palette <- c(
  "white",
  "#B4C2DE",
  "#8DA2C9",
  "#3266A7",
  "#275392",
  "#1B427D"
)


#### 2. Read statistics ####

files <- list.files(
  results_root,
  pattern = "overlap_statistics\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

files <- files[
  grepl(
    "/TCGA-.+/overlap_statistics\\.csv$",
    files
  )
]

if (length(files) == 0) {
  stop("No overlap_statistics.csv files found.")
}

data <- do.call(
  rbind,
  lapply(files, function(file) {
    x <- read.csv(
      file,
      stringsAsFactors = FALSE
    )
    
    x$cancer_type <- sub(
      "^TCGA-",
      "",
      basename(dirname(file))
    )
    
    x
  })
)

rownames(data) <- NULL

data <- data[
  data$cancer_type %in% cancer_order &
    data$overlap_class %in% class_order,
]


#### 3. Statistics ####

# Global BH correction across:
# 17 cancers × 4 classes = 68 tests.
data$BH <- p.adjust(
  data$p_value,
  method = "BH"
)

data$significant <- data$BH < 0.05

data$log2_OR <- log2(
  data$odds_ratio
)


#### 4. Standard significance symbols ####

# *       BH < 0.05
# **      BH < 0.01
# ***     BH < 0.001
# ****    BH < 0.0001
# N.S.    BH >= 0.05

data$label <- ifelse(
  data$BH < 0.0001, "****",
  ifelse(
    data$BH < 0.001, "***",
    ifelse(
      data$BH < 0.01, "**",
      ifelse(
        data$BH < 0.05, "*",
        "N.S."
      )
    )
  )
)


#### 5. Calculate actual cell colors ####

max_log2_OR <- max(
  data$log2_OR[data$significant],
  na.rm = TRUE
)

color_mapper <- scales::col_numeric(
  palette = fill_palette,
  domain = c(
    0,
    max_log2_OR
  )
)

data$cell_fill <- color_mapper(
  data$log2_OR
)


#### 6. Calculate color luminance ####

get_luminance <- function(color) {
  rgb <- grDevices::col2rgb(color) / 255
  
  rgb <- ifelse(
    rgb <= 0.04045,
    rgb / 12.92,
    ((rgb + 0.055) / 1.055)^2.4
  )
  
  as.numeric(
    0.2126 * rgb[1, ] +
      0.7152 * rgb[2, ] +
      0.0722 * rgb[3, ]
  )
}

data$luminance <- get_luminance(
  data$cell_fill
)

data$contrast_white <- 1.05 / (
  data$luminance + 0.05
)


#### 7. Select black or white text ####

# White text is used only when its contrast
# with the cell is at least 3:1.
# Otherwise black text is used.
# N.S. is always black.
data$text_color <- ifelse(
  !data$significant,
  "black",
  ifelse(
    data$contrast_white >= 3,
    "white",
    "black"
  )
)


#### 8. Set plot order ####

data$cancer_type <- factor(
  data$cancer_type,
  levels = cancer_order,
  labels = cancer_labels[cancer_order]
)

data$overlap_class <- factor(
  data$overlap_class,
  levels = rev(class_order),
  labels = rev(class_labels[class_order])
)


#### 9. Create heatmap ####

heatmap <- ggplot(
  data,
  aes(
    x = cancer_type,
    y = overlap_class
  )
) +
  
  # White background for every cell.
  geom_tile(
    fill = "white",
    color = "grey78",
    linewidth = 0.7
  ) +
  
  # Blue fill only for significant cells.
  geom_tile(
    data = data[data$significant, ],
    aes(fill = log2_OR),
    color = "grey78",
    linewidth = 0.7
  ) +
  
  # Large bold significance stars.
  geom_text(
    data = data[data$significant, ],
    aes(
      label = label,
      color = text_color
    ),
    size = 8.5,
    fontface = "bold"
  ) +
  
  # Separate N.S. layer.
  geom_text(
    data = data[!data$significant, ],
    aes(label = label),
    size = 7,
    fontface = "bold",
    color = "black"
  ) +
  
  scale_color_identity() +
  
  scale_x_discrete(
    position = "top"
  ) +
  
  scale_fill_gradientn(
    colors = fill_palette,
    limits = c(
      0,
      max_log2_OR
    ),
    oob = scales::squish,
    name = "log2 odds\nratio"
  ) +
  
  guides(
    fill = guide_colorbar(
      title.position = "top",
      
      # Left-align the legend title.
      title.hjust = 0,
      
      barheight = grid::unit(
        4,
        "cm"
      ),
      
      barwidth = grid::unit(
        0.65,
        "cm"
      ),
      
      ticks = TRUE,
      
      frame.colour = "grey60"
    )
  ) +
  
  labs(
    title = "Overlap analysis between cancer-DEGs and OIS-DEGs",
    x = NULL,
    y = NULL
  ) +
  
  coord_fixed(
    ratio = 0.95,
    clip = "off"
  ) +
  
  theme_classic(
    base_size = 17
  ) +
  
  theme(
    plot.title = element_text(
      size = 32,
      face = "bold",
      hjust = 0,
      margin = margin(
        b = 20
      )
    ),
    
    plot.title.position = "plot",
    
    axis.text.x = element_text(
      angle = 50,
      hjust = 0,
      vjust = 0,
      face = "bold",
      size = 17
    ),
    
    axis.text.y = element_text(
      face = "bold",
      size = 20,
      lineheight = 0.9
    ),
    
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold",
      size = 17,
      hjust = 0,
      lineheight = 0.9
    ),
    
    legend.text = element_text(
      size = 15
    ),
    
    legend.margin = margin(
      0, 0, 0, 12
    ),
    
    plot.margin = margin(
      18, 24, 14, 18
    )
  )


#### 10. Show heatmap ####

print(heatmap)


#### 11. Output files ####

png_file <- file.path(
  results_root,
  "TCGA_OIS_heatmap_Figma_high_quality.png"
)

csv_file <- file.path(
  results_root,
  "TCGA_OIS_heatmap_Figma_data.csv"
)


#### 12. Save high-quality PNG ####

# Final size:
# 20 inches × 8.5 inches at 600 dpi
# approximately 12,000 × 5,100 pixels.

if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(
    filename = png_file,
    plot = heatmap,
    device = ragg::agg_png,
    width = 20,
    height = 8.5,
    units = "in",
    dpi = 600,
    bg = "white"
  )
} else {
  ggsave(
    filename = png_file,
    plot = heatmap,
    device = "png",
    width = 20,
    height = 8.5,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}


#### 13. Save plotted data ####

write.csv(
  data,
  csv_file,
  row.names = FALSE
)

cat(
  "\nSaved:\n",
  "High-quality PNG: ", png_file, "\n",
  "Data: ", csv_file, "\n",
  sep = ""
)
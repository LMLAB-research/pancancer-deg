#### BLCA–OIS directional gene-overlap analysis ####

if (!requireNamespace("GeneOverlap", quietly = TRUE)) {
  stop(
    "Package 'GeneOverlap' is required. Install it with: ",
    "BiocManager::install('GeneOverlap')"
  )
}

library(GeneOverlap)


#### 1. Configuration ####

blca_file <- paste0(
  "/Users/danabiruk/Downloads/",
  "TCGA-BLCA_DGE_Raw_Results_Tumor_vs_Normal_Adjusted.csv"
)

ois_file <- paste0(
  "/Users/danabiruk/Downloads/",
  "senescence_OIS_genes.csv"
)

output_dir <- file.path(
  "outputs",
  "gene_overlap_results"
)

padj_threshold <- 0.05
lfc_threshold <- log2(1.5)


#### 2. Helper functions ####

classify_regulation <- function(
    padj,
    lfc,
    padj_cutoff = 0.05,
    lfc_cutoff = log2(1.5)
) {
  ifelse(
    !is.na(padj) & !is.na(lfc) &
      padj < padj_cutoff & lfc > lfc_cutoff,
    "UP",
    ifelse(
      !is.na(padj) & !is.na(lfc) &
        padj < padj_cutoff & lfc < -lfc_cutoff,
      "DOWN",
      "NS"
    )
  )
}


resolve_symbol_duplicates <- function(
    x,
    source_id_column,
    symbol_column = "gene_symbol",
    padj_column = "padj",
    lfc_column = "log2FoldChange",
    regulation_column = "regulation"
) {
  # Remove rows without a usable gene symbol.
  x <- x[
    !is.na(x[[symbol_column]]) & x[[symbol_column]] != "",
    ,
    drop = FALSE
  ]
  
  groups <- split(x, x[[symbol_column]])
  selected <- vector("list", length(groups))
  audit <- vector("list", length(groups))
  names(selected) <- names(groups)
  names(audit) <- names(groups)
  
  for (symbol in names(groups)) {
    group <- groups[[symbol]]
    
    significant_directions <- unique(
      group[[regulation_column]][
        group[[regulation_column]] %in% c("UP", "DOWN")
      ]
    )
    
    direction_conflict <- length(significant_directions) > 1
    
    # Representative selection order:
    # 1. non-missing padj;
    # 2. lowest padj;
    # 3. largest absolute LFC;
    # 4. lexical source ID.
    row_order <- order(
      is.na(group[[padj_column]]),
      group[[padj_column]],
      -abs(group[[lfc_column]]),
      group[[source_id_column]],
      na.last = TRUE
    )
    
    selected_row <- row_order[[1]]
    representative <- group[selected_row, , drop = FALSE]
    
    duplicate_status <- if (nrow(group) == 1) {
      "UNIQUE"
    } else if (direction_conflict) {
      "CONFLICT_OPPOSITE_DIRECTIONS"
    } else {
      "REVIEW_DUPLICATE_SAME_DIRECTION"
    }
    
    representative$duplicate_n <- nrow(group)
    representative$duplicate_status <- duplicate_status
    representative$needs_review <- nrow(group) > 1
    representative$exclude_from_overlap <- direction_conflict
    
    group$duplicate_n <- nrow(group)
    group$duplicate_status <- duplicate_status
    group$needs_review <- nrow(group) > 1
    group$direction_conflict <- direction_conflict
    group$exclude_from_overlap <- direction_conflict
    group$selected_representative <- seq_len(nrow(group)) == selected_row
    
    selected[[symbol]] <- representative
    audit[[symbol]] <- group
  }
  
  list(
    data = do.call(rbind, selected),
    audit = do.call(rbind, audit)
  )
}


run_overlap_test <- function(
    blca_genes,
    ois_genes,
    universe,
    overlap_class
) {
  # Every tested list must be a subset of the stated universe.
  blca_genes <- intersect(unique(blca_genes), universe)
  ois_genes <- intersect(unique(ois_genes), universe)
  
  test <- testGeneOverlap(
    newGeneOverlap(
      blca_genes,
      ois_genes,
      genome.size = length(universe)
    )
  )
  
  observed_overlap <- length(
    intersect(blca_genes, ois_genes)
  )
  
  expected_overlap <-
    length(blca_genes) *
    length(ois_genes) /
    length(universe)
  
  data.frame(
    overlap_class = overlap_class,
    universe_n = length(universe),
    blca_set_n = length(blca_genes),
    ois_set_n = length(ois_genes),
    observed_overlap = observed_overlap,
    expected_overlap = expected_overlap,
    fold_enrichment = observed_overlap / expected_overlap,
    odds_ratio = getOddsRatio(test),
    p_value = getPval(test),
    stringsAsFactors = FALSE
  )
}


#### 3. Read and validate input files ####

blca_deg <- read.csv(
  blca_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ois_deg <- read.csv(
  ois_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_blca_columns <- c(
  "gene_id",
  "gene_name",
  "gene_type",
  "log2FoldChange",
  "padj"
)

required_ois_columns <- c(
  "gene",
  "log2FoldChange",
  "padj"
)

missing_blca_columns <- setdiff(
  required_blca_columns,
  names(blca_deg)
)

missing_ois_columns <- setdiff(
  required_ois_columns,
  names(ois_deg)
)

if (length(missing_blca_columns) > 0) {
  stop(
    "Missing BLCA columns: ",
    paste(missing_blca_columns, collapse = ", ")
  )
}

if (length(missing_ois_columns) > 0) {
  stop(
    "Missing OIS columns: ",
    paste(missing_ois_columns, collapse = ", ")
  )
}


#### 4. Harmonize comparable genes ####

# OIS contains protein-coding genes only, so use the same gene type in BLCA.
blca_deg <- blca_deg[
  blca_deg$gene_type == "protein_coding",
  ,
  drop = FALSE
]

# Preserve the original Ensembl ID and add a version-free ID for auditing.
blca_deg$ensembl_id <- sub(
  "\\.[0-9]+$",
  "",
  blca_deg$gene_id
)

# Create matching symbol columns without overwriting the original columns.
blca_deg$gene_symbol <- toupper(
  trimws(as.character(blca_deg$gene_name))
)

ois_deg$gene_symbol <- toupper(
  trimws(as.character(ois_deg$gene))
)

blca_deg$gene_symbol[blca_deg$gene_symbol == ""] <- NA_character_
ois_deg$gene_symbol[ois_deg$gene_symbol == ""] <- NA_character_


#### 5. Apply the same DEG definition ####

blca_deg$regulation <- classify_regulation(
  blca_deg$padj,
  blca_deg$log2FoldChange,
  padj_cutoff = padj_threshold,
  lfc_cutoff = lfc_threshold
)

ois_deg$regulation <- classify_regulation(
  ois_deg$padj,
  ois_deg$log2FoldChange,
  padj_cutoff = padj_threshold,
  lfc_cutoff = lfc_threshold
)


#### 6. Resolve duplicated gene symbols ####

blca_resolved <- resolve_symbol_duplicates(
  blca_deg,
  source_id_column = "ensembl_id"
)

ois_resolved <- resolve_symbol_duplicates(
  ois_deg,
  source_id_column = "gene_symbol"
)

blca_deg <- blca_resolved$data[
  !blca_resolved$data$exclude_from_overlap,
  ,
  drop = FALSE
]

ois_deg <- ois_resolved$data[
  !ois_resolved$data$exclude_from_overlap,
  ,
  drop = FALSE
]


#### 7. Define the common statistical universe ####

# Only genes present in both processed tables can contribute to an overlap.
common_tested_genes <- intersect(
  blca_deg$gene_symbol,
  ois_deg$gene_symbol
)

universe_size <- length(common_tested_genes)

if (universe_size == 0) {
  stop("The two datasets have no common tested genes.")
}


#### 8. Create directional DEG sets ####

# Restrict every statistical list to the same common universe.
blca_up <- intersect(
  blca_deg$gene_symbol[blca_deg$regulation == "UP"],
  common_tested_genes
)

blca_down <- intersect(
  blca_deg$gene_symbol[blca_deg$regulation == "DOWN"],
  common_tested_genes
)

ois_up <- intersect(
  ois_deg$gene_symbol[ois_deg$regulation == "UP"],
  common_tested_genes
)

ois_down <- intersect(
  ois_deg$gene_symbol[ois_deg$regulation == "DOWN"],
  common_tested_genes
)

up_up <- intersect(blca_up, ois_up)
up_down <- intersect(blca_up, ois_down)
down_up <- intersect(blca_down, ois_up)
down_down <- intersect(blca_down, ois_down)


#### 9. Test overlap enrichment ####

overlap_statistics <- rbind(
  run_overlap_test(
    blca_up,
    ois_up,
    common_tested_genes,
    "UP-UP"
  ),
  run_overlap_test(
    blca_up,
    ois_down,
    common_tested_genes,
    "UP-DOWN"
  ),
  run_overlap_test(
    blca_down,
    ois_up,
    common_tested_genes,
    "DOWN-UP"
  ),
  run_overlap_test(
    blca_down,
    ois_down,
    common_tested_genes,
    "DOWN-DOWN"
  )
)

overlap_statistics$padj_BH <- p.adjust(
  overlap_statistics$p_value,
  method = "BH"
)

overlap_statistics$significant_BH <-
  overlap_statistics$padj_BH < 0.05


#### 10. Build the detailed overlap table ####

common_genes <- data.frame(
  gene_symbol = c(
    up_up,
    up_down,
    down_up,
    down_down
  ),
  overlap_class = c(
    rep("UP-UP", length(up_up)),
    rep("UP-DOWN", length(up_down)),
    rep("DOWN-UP", length(down_up)),
    rep("DOWN-DOWN", length(down_down))
  ),
  stringsAsFactors = FALSE
)

if (anyDuplicated(common_genes$gene_symbol)) {
  stop("At least one gene was assigned to more than one overlap class.")
}

blca_for_merge <- blca_deg[
  blca_deg$gene_symbol %in% common_genes$gene_symbol,
  c(
    "gene_symbol",
    "ensembl_id",
    "gene_name",
    "log2FoldChange",
    "padj",
    "regulation",
    "duplicate_n",
    "duplicate_status",
    "needs_review"
  ),
  drop = FALSE
]

ois_for_merge <- ois_deg[
  ois_deg$gene_symbol %in% common_genes$gene_symbol,
  c(
    "gene_symbol",
    "gene",
    "log2FoldChange",
    "padj",
    "regulation",
    "duplicate_n",
    "duplicate_status",
    "needs_review"
  ),
  drop = FALSE
]

overlap_table <- merge(
  blca_for_merge,
  ois_for_merge,
  by = "gene_symbol",
  suffixes = c("_BLCA", "_OIS"),
  all = FALSE,
  sort = FALSE
)

overlap_table <- merge(
  overlap_table,
  common_genes,
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)

overlap_table$needs_manual_review <-
  overlap_table$needs_review_BLCA |
  overlap_table$needs_review_OIS

class_order <- c(
  "UP-UP",
  "UP-DOWN",
  "DOWN-UP",
  "DOWN-DOWN"
)

overlap_table$overlap_class <- factor(
  overlap_table$overlap_class,
  levels = class_order
)

overlap_table <- overlap_table[
  order(
    overlap_table$overlap_class,
    overlap_table$gene_symbol
  ),
  ,
  drop = FALSE
]

overlap_table$overlap_class <- as.character(
  overlap_table$overlap_class
)

overlap_counts <- data.frame(
  overlap_class = class_order,
  n_genes = c(
    length(up_up),
    length(up_down),
    length(down_up),
    length(down_down)
  ),
  stringsAsFactors = FALSE
)


#### 11. Save results ####

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  overlap_table,
  file.path(output_dir, "overlap_all_classes.csv"),
  row.names = FALSE
)

write.csv(
  overlap_counts,
  file.path(output_dir, "overlap_counts.csv"),
  row.names = FALSE
)

write.csv(
  overlap_statistics,
  file.path(output_dir, "overlap_statistics.csv"),
  row.names = FALSE
)

for (class_name in class_order) {
  class_table <- overlap_table[
    overlap_table$overlap_class == class_name,
    ,
    drop = FALSE
  ]
  
  write.csv(
    class_table,
    file.path(
      output_dir,
      paste0("overlap_", class_name, ".csv")
    ),
    row.names = FALSE
  )
}

write.csv(
  blca_resolved$audit[blca_resolved$audit$needs_review, ],
  file.path(output_dir, "BLCA_duplicates_for_review.csv"),
  row.names = FALSE
)

write.csv(
  ois_resolved$audit[ois_resolved$audit$needs_review, ],
  file.path(output_dir, "OIS_duplicates_for_review.csv"),
  row.names = FALSE
)


#### 12. Print concise summary ####

cat("Common tested-gene universe:", universe_size, "genes\n\n")
print(overlap_counts, row.names = FALSE)
cat("\n")
print(overlap_statistics, row.names = FALSE)
cat("\nResults saved to:", normalizePath(output_dir), "\n")




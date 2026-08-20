script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1]])
} else {
  sys.frame(1)$ofile
}
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "config.R"))

if (!requireNamespace("GeneOverlap", quietly = TRUE)) {
  stop(
    "Package 'GeneOverlap' is required. Install it with: ",
    "BiocManager::install('GeneOverlap')"
  )
}

library(GeneOverlap)

if (!dir.exists(cancer_input_dir)) {
  stop("Input directory does not exist: ", cancer_input_dir)
}


classify_regulation <- function(padj, lfc) {
  ifelse(
    !is.na(padj) & !is.na(lfc) &
      padj < padj_threshold & lfc > lfc_threshold,
    "UP",
    ifelse(
      !is.na(padj) & !is.na(lfc) &
        padj < padj_threshold & lfc < -lfc_threshold,
      "DOWN",
      "NS"
    )
  )
}


validate_numeric <- function(x, column, file) {
  converted <- suppressWarnings(as.numeric(x))
  invalid <- !is.na(x) & trimws(as.character(x)) != "" & is.na(converted)
  if (any(invalid)) {
    stop("Non-numeric values in column '", column, "' of: ", file)
  }
  converted
}


resolve_symbol_duplicates <- function(
    x,
    source_id_column,
    symbol_column = "gene_symbol",
    padj_column = "padj",
    lfc_column = "log2FoldChange",
    regulation_column = "regulation"
) {
  x <- x[
    !is.na(x[[symbol_column]]) & x[[symbol_column]] != "",
    ,
    drop = FALSE
  ]

  if (!nrow(x)) {
    stop("No valid gene symbols remain after filtering.")
  }
  
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
    cancer_genes,
    ois_genes,
    universe,
    overlap_class
) {
  cancer_genes <- intersect(unique(cancer_genes), universe)
  ois_genes <- intersect(unique(ois_genes), universe)
  
  test <- testGeneOverlap(
    newGeneOverlap(
      cancer_genes,
      ois_genes,
      genome.size = length(universe)
    )
  )
  
  observed_overlap <- length(
    intersect(cancer_genes, ois_genes)
  )
  
  expected_overlap <-
    length(cancer_genes) * length(ois_genes) / length(universe)
  
  data.frame(
    overlap_class = overlap_class,
    universe_n = length(universe),
    cancer_set_n = length(cancer_genes),
    ois_set_n = length(ois_genes),
    observed_overlap = observed_overlap,
    expected_overlap = expected_overlap,
    fold_enrichment = if (expected_overlap > 0) {
      observed_overlap / expected_overlap
    } else {
      NA_real_
    },
    odds_ratio = getOddsRatio(test),
    p_value = getPval(test),
    stringsAsFactors = FALSE
  )
}


if (!file.exists(ois_file)) {
  stop("OIS file does not exist: ", ois_file)
}

ois_original <- read.csv(
  ois_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_ois_columns <- c(
  "gene",
  "log2FoldChange",
  "padj"
)

missing_ois_columns <- setdiff(
  required_ois_columns,
  names(ois_original)
)

if (length(missing_ois_columns) > 0) {
  stop(
    "Missing OIS columns: ",
    paste(missing_ois_columns, collapse = ", ")
  )
}

ois_original$gene_symbol <- toupper(
  trimws(as.character(ois_original$gene))
)

ois_original$gene_symbol[ois_original$gene_symbol == ""] <- NA_character_
ois_original$source_row <- seq_len(nrow(ois_original))
ois_original$padj <- validate_numeric(ois_original$padj, "padj", ois_file)
ois_original$log2FoldChange <- validate_numeric(
  ois_original$log2FoldChange,
  "log2FoldChange",
  ois_file
)

ois_original$regulation <- classify_regulation(
  ois_original$padj,
  ois_original$log2FoldChange
)

ois_resolved <- resolve_symbol_duplicates(
  ois_original,
  source_id_column = "source_row"
)

ois_prepared <- ois_resolved$data[
  !ois_resolved$data$exclude_from_overlap,
  ,
  drop = FALSE
]


cancer_files <- list.files(
  cancer_input_dir,
  pattern = cancer_file_pattern,
  full.names = TRUE
)

cancer_files <- sort(cancer_files)

if (length(cancer_files) == 0) {
  stop("No matching TCGA cancer files found in: ", cancer_input_dir)
}

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  ois_resolved$audit[ois_resolved$audit$needs_review, , drop = FALSE],
  file.path(results_dir, "OIS_duplicates_for_review.csv"),
  row.names = FALSE
)


analyze_one_cancer <- function(cancer_file) {
  cancer_project <- sub(
    "_DGE_Raw_Results_Tumor_vs_Normal_Adjusted\\.csv$",
    "",
    basename(cancer_file)
  )
  
  cancer_output_dir <- file.path(
    results_dir,
    cancer_project
  )
  
  cancer_deg <- read.csv(
    cancer_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  required_cancer_columns <- c(
    "gene_id",
    "gene_name",
    "gene_type",
    "log2FoldChange",
    "padj"
  )
  
  missing_columns <- setdiff(
    required_cancer_columns,
    names(cancer_deg)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "Missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  # Match the protein-coding OIS gene universe.
  cancer_deg <- cancer_deg[
    cancer_deg$gene_type == "protein_coding",
    ,
    drop = FALSE
  ]
  
  cancer_deg$ensembl_id <- sub(
    "\\.[0-9]+$",
    "",
    cancer_deg$gene_id
  )
  
  cancer_deg$gene_symbol <- toupper(
    trimws(as.character(cancer_deg$gene_name))
  )
  
  cancer_deg$gene_symbol[cancer_deg$gene_symbol == ""] <- NA_character_
  cancer_deg$padj <- validate_numeric(cancer_deg$padj, "padj", cancer_file)
  cancer_deg$log2FoldChange <- validate_numeric(
    cancer_deg$log2FoldChange,
    "log2FoldChange",
    cancer_file
  )
  
  cancer_deg$regulation <- classify_regulation(
    cancer_deg$padj,
    cancer_deg$log2FoldChange
  )
  
  cancer_resolved <- resolve_symbol_duplicates(
    cancer_deg,
    source_id_column = "ensembl_id"
  )
  
  cancer_prepared <- cancer_resolved$data[
    !cancer_resolved$data$exclude_from_overlap,
    ,
    drop = FALSE
  ]
  
  common_tested_genes <- intersect(
    cancer_prepared$gene_symbol,
    ois_prepared$gene_symbol
  )
  
  if (length(common_tested_genes) == 0) {
    stop("No common tested genes with OIS.")
  }
  
  cancer_up <- intersect(
    cancer_prepared$gene_symbol[cancer_prepared$regulation == "UP"],
    common_tested_genes
  )
  
  cancer_down <- intersect(
    cancer_prepared$gene_symbol[cancer_prepared$regulation == "DOWN"],
    common_tested_genes
  )
  
  ois_up <- intersect(
    ois_prepared$gene_symbol[ois_prepared$regulation == "UP"],
    common_tested_genes
  )
  
  ois_down <- intersect(
    ois_prepared$gene_symbol[ois_prepared$regulation == "DOWN"],
    common_tested_genes
  )
  
  up_up <- sort(intersect(cancer_up, ois_up))
  up_down <- sort(intersect(cancer_up, ois_down))
  down_up <- sort(intersect(cancer_down, ois_up))
  down_down <- sort(intersect(cancer_down, ois_down))
  
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
  
  overlap_statistics <- rbind(
    run_overlap_test(cancer_up, ois_up, common_tested_genes, "UP-UP"),
    run_overlap_test(cancer_up, ois_down, common_tested_genes, "UP-DOWN"),
    run_overlap_test(cancer_down, ois_up, common_tested_genes, "DOWN-UP"),
    run_overlap_test(cancer_down, ois_down, common_tested_genes, "DOWN-DOWN")
  )
  
  overlap_statistics$BH_within_cancer <- p.adjust(
    overlap_statistics$p_value,
    method = "BH"
  )
  overlap_statistics$significant_BH_within_cancer <-
    overlap_statistics$BH_within_cancer < 0.05
  
  common_genes <- data.frame(
    gene_symbol = c(up_up, up_down, down_up, down_down),
    overlap_class = c(
      rep("UP-UP", length(up_up)),
      rep("UP-DOWN", length(up_down)),
      rep("DOWN-UP", length(down_up)),
      rep("DOWN-DOWN", length(down_down))
    ),
    stringsAsFactors = FALSE
  )
  
  if (anyDuplicated(common_genes$gene_symbol)) {
    stop("A gene was assigned to multiple overlap classes.")
  }
  
  cancer_for_merge <- cancer_prepared[
    cancer_prepared$gene_symbol %in% common_genes$gene_symbol,
    c(
      "gene_symbol", "ensembl_id", "gene_name", "log2FoldChange",
      "padj", "regulation", "duplicate_n", "duplicate_status",
      "needs_review"
    ),
    drop = FALSE
  ]
  
  ois_for_merge <- ois_prepared[
    ois_prepared$gene_symbol %in% common_genes$gene_symbol,
    c(
      "gene_symbol", "gene", "log2FoldChange", "padj", "regulation",
      "duplicate_n", "duplicate_status", "needs_review"
    ),
    drop = FALSE
  ]
  
  overlap_table <- merge(
    cancer_for_merge,
    ois_for_merge,
    by = "gene_symbol",
    suffixes = c("_CANCER", "_OIS"),
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
    overlap_table$needs_review_CANCER |
    overlap_table$needs_review_OIS
  
  overlap_table$overlap_class <- factor(
    overlap_table$overlap_class,
    levels = class_order
  )
  
  overlap_table <- overlap_table[
    order(overlap_table$overlap_class, overlap_table$gene_symbol),
    ,
    drop = FALSE
  ]
  
  overlap_table$overlap_class <- as.character(
    overlap_table$overlap_class
  )
  
  dir.create(
    cancer_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  write.csv(
    overlap_table,
    file.path(cancer_output_dir, "overlap_all_classes.csv"),
    row.names = FALSE
  )
  
  write.csv(
    overlap_counts,
    file.path(cancer_output_dir, "overlap_counts.csv"),
    row.names = FALSE
  )
  
  write.csv(
    overlap_statistics,
    file.path(cancer_output_dir, "overlap_statistics.csv"),
    row.names = FALSE
  )
  
  # Exact background used in every overlap test for this cancer.
  write.csv(
    data.frame(
      gene_symbol = sort(unique(common_tested_genes)),
      stringsAsFactors = FALSE
    ),
    file.path(cancer_output_dir, "common_tested_genes.csv"),
    row.names = FALSE
  )
  
  for (class_name in class_order) {
    write.csv(
      overlap_table[
        overlap_table$overlap_class == class_name,
        ,
        drop = FALSE
      ],
      file.path(
        cancer_output_dir,
        paste0("overlap_", class_name, ".csv")
      ),
      row.names = FALSE
    )
  }
  
  write.csv(
    cancer_resolved$audit[cancer_resolved$audit$needs_review, ],
    file.path(cancer_output_dir, "CANCER_duplicates_for_review.csv"),
    row.names = FALSE
  )
  
  data.frame(
    cancer_project = cancer_project,
    input_file = cancer_file,
    output_dir = cancer_output_dir,
    status = "COMPLETED",
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}


run_status <- vector("list", length(cancer_files))

for (i in seq_along(cancer_files)) {
  cancer_file <- cancer_files[[i]]
  cancer_name <- sub(
    "_DGE_Raw_Results_Tumor_vs_Normal_Adjusted\\.csv$",
    "",
    basename(cancer_file)
  )
  
  message(
    "[", i, "/", length(cancer_files), "] Running ",
    cancer_name
  )
  
  run_status[[i]] <- tryCatch(
    analyze_one_cancer(cancer_file),
    error = function(error) {
      warning(
        "Analysis failed for ", cancer_name, ": ",
        conditionMessage(error)
      )
      
      data.frame(
        cancer_project = cancer_name,
        input_file = cancer_file,
        output_dir = file.path(results_dir, cancer_name),
        status = "FAILED",
        error_message = conditionMessage(error),
        stringsAsFactors = FALSE
      )
    }
  )
}

run_status <- do.call(rbind, run_status)

write.csv(
  run_status,
  file.path(results_dir, "batch_run_status.csv"),
  row.names = FALSE
)

cat("\nProcessing complete.\n")
print(
  run_status[, c("cancer_project", "status")],
  row.names = FALSE
)
cat("\nResult folders saved in:\n", results_dir, "\n")

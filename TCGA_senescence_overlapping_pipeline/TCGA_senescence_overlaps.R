script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else sys.frame(1)$ofile
project_dir <- dirname(normalizePath(script_file))
source(file.path(project_dir, "config.R"))

status_file <- file.path(results_dir, "batch_run_status.csv")
if (file.exists(status_file) && !file.remove(status_file)) stop("Cannot replace batch status.")

if (!requireNamespace("GeneOverlap", quietly = TRUE)) {
  stop("Install GeneOverlap with BiocManager::install('GeneOverlap')")
}

library(GeneOverlap)

# Assign UP, DOWN, or NS using adjusted p-value and fold-change thresholds
classify_regulation <- function(padj, lfc) {
  ifelse(
    !is.na(padj) & !is.na(lfc) &
      padj < padj_threshold & lfc > lfc_threshold,
    "UP",
    ifelse( !is.na(padj) & !is.na(lfc) & padj < padj_threshold & lfc < -lfc_threshold, "DOWN", "NS" )
  )
}

validate_numeric <- function(x, column, file) {
  y <- suppressWarnings(as.numeric(x))

  bad <- !is.na(x) &
    trimws(as.character(x)) != "" &
    is.na(y)

  if (any(bad)) {
    stop("Non-numeric values in ", column, ": ", file)
  }

  if (any(!is.na(y) & !is.finite(y))) stop("Non-finite values in ", column, ": ", file)
  if (column == "padj" && any(y < 0 | y > 1, na.rm = TRUE)) {
    stop("Adjusted p-values must be between 0 and 1: ", file)
  }
  y
}

# One row per gene symbol and retain all duplicate rows for review
resolve_symbol_duplicates <- function(
    x,
    source_id_column,
    symbol_column = "gene_symbol"
) {
  x <- x[
    !is.na(x[[symbol_column]]) &
      x[[symbol_column]] != "",
    ,
    drop = FALSE
  ]

  if (!nrow(x)) {
    stop("No valid gene symbols remain.")
  }

  groups <- split(x, x[[symbol_column]])
  selected <- list()
  audit <- list()

  for (symbol in names(groups)) {
    group <- groups[[symbol]]

    directions <- unique( group$regulation[ group$regulation %in% c("UP", "DOWN") ] )

    # Opposite significant directions exclude this symbol from the analysis.
    conflict <- length(directions) > 1

    # Prefer lowest padj, then largest absolute fold change, then source ID.
    row_order <- order(
      is.na(group$padj),
      group$padj,
      -abs(group$log2FoldChange),
      group[[source_id_column]],
      na.last = TRUE
    )

    selected_row <- row_order[[1]]
    representative <- group[selected_row, , drop = FALSE]

    status <- if (nrow(group) == 1) {
      "UNIQUE"
    } else if (conflict) {
      "CONFLICT_OPPOSITE_DIRECTIONS"
    } else {
      "REVIEW_DUPLICATE_SAME_DIRECTION"
    }

    representative$duplicate_n <- nrow(group)
    representative$duplicate_status <- status
    representative$needs_review <- nrow(group) > 1
    representative$exclude_from_overlap <- conflict

    group$duplicate_n <- nrow(group)
    group$duplicate_status <- status
    group$needs_review <- nrow(group) > 1
    group$direction_conflict <- conflict
    group$exclude_from_overlap <- conflict
    group$selected_representative <-
      seq_len(nrow(group)) == selected_row

    selected[[symbol]] <- representative
    audit[[symbol]] <- group
  }

  list( data = do.call(rbind, selected), audit = do.call(rbind, audit) )
}

# Test whether two directional gene sets overlap more than expected by chance.
run_overlap_test <- function(
    cancer_genes,
    signature_genes,
    universe,
    overlap_class
) {
  cancer_genes <- intersect(unique(cancer_genes), universe)
  signature_genes <- intersect(unique(signature_genes), universe)

  test <- testGeneOverlap(
    newGeneOverlap( cancer_genes, signature_genes, genome.size = length(universe) )
  )

  observed <- length( intersect(cancer_genes, signature_genes) )

  expected <-
    length(cancer_genes) *
    length(signature_genes) /
    length(universe)

  data.frame(
    overlap_class = overlap_class,
    universe_n = length(universe),
    cancer_set_n = length(cancer_genes),
    signature_set_n = length(signature_genes),
    observed_overlap = observed,
    expected_overlap = expected,
    fold_enrichment = if (expected > 0) observed / expected else NA_real_,
    odds_ratio = getOddsRatio(test),
    p_value = getPval(test),
    stringsAsFactors = FALSE
  )
}

read_background <- function(file) {
  x <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  columns <- intersect(c("gene", "gene_symbol", "gene_name"), names(x))
  if (!length(columns) && ncol(x) == 1) columns <- names(x)
  if (!length(columns)) {
    stop("Background needs a gene, gene_symbol, or gene_name column: ", file)
  }
  genes <- unique(toupper(trimws(as.character(x[[columns[1]]]))))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes)) stop("Empty background: ", file)
  genes
}

prepare_genes <- function(file, cancer = FALSE) {
  x <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("log2FoldChange", "padj", if (cancer) {
    c("gene_id", "gene_name", "gene_type")
  } else "gene")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("Missing columns in ", file, ": ", paste(missing, collapse = ", "))
  }
  if (cancer) {
    x <- x[!is.na(x$gene_type) & x$gene_type == "protein_coding", , drop = FALSE]
    x$ensembl_id <- sub("\\.[0-9]+$", "", x$gene_id)
  } else {
    x$source_row <- seq_len(nrow(x))
  }
  x$gene_symbol <- toupper(trimws(as.character(x[[if (cancer) "gene_name" else "gene"]])))
  x$gene_symbol[x$gene_symbol == ""] <- NA_character_
  x$padj <- validate_numeric(x$padj, "padj", file)
  x$log2FoldChange <- validate_numeric(x$log2FoldChange, "log2FoldChange", file)
  if (cancer) {
    # Only genes with usable DEG statistics could enter a directional set.
    x <- x[!is.na(x$padj) & !is.na(x$log2FoldChange), , drop = FALSE]
  }
  x$regulation <- classify_regulation(x$padj, x$log2FoldChange)
  resolved <- resolve_symbol_duplicates(x, if (cancer) "ensembl_id" else "source_row")
  list(
    data = resolved$data[!resolved$data$exclude_from_overlap, , drop = FALSE],
    audit = resolved$audit
  )
}

prepare_signature <- function(file, background_file) {
  signature <- prepare_genes(file)
  background <- read_background(background_file)
  outside <- setdiff(signature$audit$gene_symbol, background)
  if (length(outside)) {
    stop("Signature genes absent from background: ", paste(head(outside, 10), collapse = ", "),
         " in ", background_file)
  }
  conflicts <- signature$audit$gene_symbol[signature$audit$exclude_from_overlap]
  signature$background <- setdiff(background, conflicts)
  signature
}

prepare_cancer <- function(file) {
  prepare_genes(file, cancer = TRUE)
}

# Compare one cancer with one signature and export its overlap results
analyze_one_cancer <- function(
    cancer_file,
    signature,
    signature_name,
    output_root
) {
  cancer_project <- sub( "_DGE_Raw_Results_Tumor_vs_Normal_Adjusted\\.csv$", "", basename(cancer_file) )

  output_dir <- file.path( output_root, cancer_project )

  cancer <- prepare_cancer(cancer_file)

  universe <- intersect( cancer$data$gene_symbol, signature$background )

  if (!length(universe)) {
    stop("No common tested genes.")
  }

  cancer_sets <- lapply(c("UP", "DOWN"), function(direction) {
    intersect(cancer$data$gene_symbol[cancer$data$regulation == direction], universe)
  })
  signature_sets <- lapply(c("UP", "DOWN"), function(direction) {
    intersect(signature$data$gene_symbol[signature$data$regulation == direction], universe)
  })
  names(cancer_sets) <- names(signature_sets) <- c("UP", "DOWN")

  directions <- strsplit(class_order, "-", fixed = TRUE)
  overlaps <- lapply(directions, function(d) {
    sort(intersect(cancer_sets[[d[1]]], signature_sets[[d[2]]]))
  })
  names(overlaps) <- class_order
  counts <- data.frame(overlap_class = class_order, n_genes = lengths(overlaps))
  statistics <- do.call(rbind, lapply(seq_along(class_order), function(i) {
    d <- directions[[i]]
    run_overlap_test(cancer_sets[[d[1]]], signature_sets[[d[2]]], universe, class_order[i])
  }))

  statistics$signature <- signature_name

  # Correct the four directional tests within this cancer/signature pair
  statistics$BH_within_cancer <- p.adjust( statistics$p_value, method = "BH" )

  genes <- data.frame(
    gene_symbol = unlist(overlaps, use.names = FALSE),
    overlap_class = rep(class_order, lengths(overlaps))
  )

  if (anyDuplicated(genes$gene_symbol)) {
    stop("Gene assigned to multiple overlap classes.")
  }

  # Attach fold changes, p-values, and duplicate flags to the shared genes
  cancer_table <- cancer$data[
    cancer$data$gene_symbol %in% genes$gene_symbol,
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

  signature_table <- signature$data[
    signature$data$gene_symbol %in% genes$gene_symbol,
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

  names(signature_table)[
    names(signature_table) == "gene"
  ] <- "signature_gene"

  overlap_table <- merge(
    cancer_table,
    signature_table,
    by = "gene_symbol",
    suffixes = c("_CANCER", "_SIGNATURE"),
    all = FALSE,
    sort = FALSE
  )

  overlap_table <- merge( overlap_table, genes, by = "gene_symbol", all.x = TRUE, sort = FALSE )

  overlap_table$signature <- rep(signature_name, nrow(overlap_table))

  overlap_table$needs_manual_review <-
    overlap_table$needs_review_CANCER |
    overlap_table$needs_review_SIGNATURE

  overlap_table$overlap_class <- factor( overlap_table$overlap_class, levels = class_order )

  overlap_table <- overlap_table[
    order( overlap_table$overlap_class, overlap_table$gene_symbol ),
    ,
    drop = FALSE
  ]

  overlap_table$overlap_class <-
    as.character(overlap_table$overlap_class)

  # Save combined and per-class gene lists
  dir.create( output_dir, recursive = TRUE, showWarnings = FALSE )

  tables <- list(
    overlap_all_classes = overlap_table,
    overlap_counts = counts,
    overlap_statistics = statistics,
    common_tested_genes = data.frame(gene_symbol = sort(unique(universe))),
    CANCER_duplicates_for_review = cancer$audit[cancer$audit$needs_review, , drop = FALSE]
  )
  for (class_name in class_order) {
    tables[[paste0("overlap_", class_name)]] <-
      overlap_table[overlap_table$overlap_class == class_name, , drop = FALSE]
  }
  for (name in names(tables)) {
    write.csv(tables[[name]], file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
  }

  data.frame(
    signature = signature_name,
    cancer_project = cancer_project,
    status = "COMPLETED",
    error_message = NA_character_
  )
}

for (folder in c(senescence_input_dir, cancer_input_dir)) {
  if (!dir.exists(folder)) stop("Missing input directory: ", folder)
}
signature_files <- list.files(
  senescence_input_dir, "^(senescence|quiescence)_[A-Za-z0-9_-]+_genes\\.csv$",
  full.names = TRUE, ignore.case = TRUE
)
signature_names <- toupper(sub(
  "_genes\\.csv$", "", sub("^(senescence|quiescence)_", "", basename(signature_files),
                          ignore.case = TRUE), ignore.case = TRUE
))
keep <- signature_names %in% expected_signatures
signature_files <- signature_files[keep]
signature_names <- signature_names[keep]
if (anyDuplicated(signature_names)) stop("Duplicate signature filenames found.")
names(signature_files) <- signature_names
missing <- setdiff(expected_signatures, signature_names)
if (length(missing)) stop("Missing signatures: ", paste(missing, collapse = ", "))
signature_files <- signature_files[expected_signatures]
background_files <- sub("_genes\\.csv$", "_background.csv", signature_files, ignore.case = TRUE)
names(background_files) <- expected_signatures
if (any(!file.exists(background_files))) {
  stop("Missing backgrounds: ", paste(background_files[!file.exists(background_files)], collapse = ", "))
}
cancer_files <- sort(list.files(cancer_input_dir, cancer_file_pattern, full.names = TRUE))
if (!length(cancer_files)) stop("No matching TCGA files in: ", cancer_input_dir)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Process each signature once, then record the outcome for every cancer.
run_status <- list()
for (signature_name in expected_signatures) {
  message("Running ", signature_name)
  signature <- prepare_signature(signature_files[[signature_name]], background_files[[signature_name]])
  output_root <- file.path(results_dir, signature_name)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  write.csv(
    signature$audit[signature$audit$needs_review, , drop = FALSE],
    file.path(output_root, paste0(signature_name, "_duplicates_for_review.csv")), row.names = FALSE
  )
  for (cancer_file in cancer_files) {
    cancer_project <- sub("_DGE_Raw_Results_Tumor_vs_Normal_Adjusted\\.csv$", "", basename(cancer_file))
    message("  ", cancer_project)
    run_status[[length(run_status) + 1]] <- tryCatch(
      analyze_one_cancer(cancer_file, signature, signature_name, output_root),
      error = function(e) {
        warning(signature_name, " / ", cancer_project, ": ", conditionMessage(e))
        data.frame(signature = signature_name, cancer_project = cancer_project,
                   status = "FAILED", error_message = conditionMessage(e))
      }
    )
  }
}
run_status <- do.call(rbind, run_status)
write.csv(run_status, file.path(results_dir, "batch_run_status.csv"), row.names = FALSE)
print(run_status[, c("signature", "cancer_project", "status")], row.names = FALSE)
if (any(run_status$status == "FAILED")) stop("Some comparisons failed. See batch_run_status.csv.")

#### TCGA cancer DEG versus OIS DEG directional overlap analysis ####

# Run from the repository root with:
# Rscript --vanilla overlapping-analyses/TCGA_OIS_overlaps_fisher.R
#
# The test universe is cancer-specific: protein-coding genes with usable DEG
# statistics that also occur in the supplied OIS tested-gene background.
# Fisher's exact tests are one-sided tests for enrichment (odds ratio > 1).


#### 1. Configuration ####

cancer_input_dir <- file.path("results", "tables", "TCGA_deg")
ois_results_file <- file.path(
  "data", "senescence_meta_analysis_deg", "senescence_OIS_genes.csv"
)
ois_background_file <- file.path(
  "data", "senescence_meta_analysis_deg", "senescence_OIS_background.csv"
)
output_root <- file.path(
  "results", "tables", "overlapping-analyses", "TCGA_OIS_overlap_fisher"
)

cancer_file_pattern <- paste0(
  "^TCGA-[A-Za-z0-9-]+_DGE_Raw_Results_",
  "Tumor_vs_Normal_Adjusted\\.csv$"
)

padj_threshold <- 0.05
lfc_threshold <- log2(1.5)
class_order <- c("UP-UP", "UP-DOWN", "DOWN-UP", "DOWN-DOWN")


#### 2. Helpers ####

read_csv_checked <- function(path, required_columns) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path)
  }

  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing_columns <- setdiff(required_columns, names(x))

  if (length(missing_columns) > 0) {
    stop(
      "Missing columns in ", path, ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  x
}


clean_symbol <- function(x) {
  symbol <- toupper(trimws(as.character(x)))
  symbol[is.na(symbol) | symbol == ""] <- NA_character_
  symbol
}


classify_regulation <- function(padj, lfc) {
  regulation <- rep("NS", length(padj))
  regulation[padj < padj_threshold & lfc > lfc_threshold] <- "UP"
  regulation[padj < padj_threshold & lfc < -lfc_threshold] <- "DOWN"
  regulation
}


resolve_symbol_duplicates <- function(
    x,
    source_id_column,
    symbol_column = "gene_symbol"
) {
  x <- x[!is.na(x[[symbol_column]]), , drop = FALSE]
  groups <- split(x, x[[symbol_column]])

  selected <- lapply(groups, function(group) {
    # Use one deterministic representative per symbol. Selection does not
    # remove a gene from the universe based on its observed DEG direction.
    row_order <- order(
      group$padj,
      -abs(group$log2FoldChange),
      as.character(group[[source_id_column]]),
      na.last = TRUE
    )
    representative <- group[row_order[[1]], , drop = FALSE]
    significant_directions <- unique(
      group$regulation[group$regulation %in% c("UP", "DOWN")]
    )

    representative$duplicate_n <- nrow(group)
    representative$duplicate_direction_conflict <-
      length(significant_directions) > 1
    representative
  })

  audit <- do.call(rbind, lapply(groups, function(group) {
    if (nrow(group) == 1) {
      return(NULL)
    }

    row_order <- order(
      group$padj,
      -abs(group$log2FoldChange),
      as.character(group[[source_id_column]]),
      na.last = TRUE
    )
    significant_directions <- unique(
      group$regulation[group$regulation %in% c("UP", "DOWN")]
    )
    group$duplicate_n <- nrow(group)
    group$duplicate_direction_conflict <- length(significant_directions) > 1
    group$selected_representative <- seq_len(nrow(group)) == row_order[[1]]
    group
  }))

  list(
    data = do.call(rbind, selected),
    audit = audit
  )
}


run_fisher_overlap <- function(
    cancer_genes,
    ois_genes,
    universe,
    overlap_class
) {
  universe <- unique(universe)
  cancer_genes <- intersect(unique(cancer_genes), universe)
  ois_genes <- intersect(unique(ois_genes), universe)

  overlap_n <- length(intersect(cancer_genes, ois_genes))
  cancer_only_n <- length(cancer_genes) - overlap_n
  ois_only_n <- length(ois_genes) - overlap_n
  neither_n <- length(universe) - overlap_n - cancer_only_n - ois_only_n

  contingency_table <- matrix(
    c(overlap_n, cancer_only_n, ois_only_n, neither_n),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      cancer_set = c("in", "out"),
      ois_set = c("in", "out")
    )
  )

  fisher_result <- fisher.test(
    contingency_table,
    alternative = "greater",
    conf.level = 0.95
  )

  expected_overlap <- length(cancer_genes) * length(ois_genes) /
    length(universe)

  data.frame(
    overlap_class = overlap_class,
    universe_n = length(universe),
    cancer_set_n = length(cancer_genes),
    ois_set_n = length(ois_genes),
    overlap_n = overlap_n,
    cancer_only_n = cancer_only_n,
    ois_only_n = ois_only_n,
    neither_n = neither_n,
    expected_overlap = expected_overlap,
    fold_enrichment = if (expected_overlap > 0) {
      overlap_n / expected_overlap
    } else {
      NA_real_
    },
    odds_ratio = unname(fisher_result$estimate),
    odds_ratio_ci_lower = unname(fisher_result$conf.int[[1]]),
    odds_ratio_ci_upper = unname(fisher_result$conf.int[[2]]),
    p_value = fisher_result$p.value,
    alternative = "greater",
    stringsAsFactors = FALSE
  )
}


write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "")
}


#### 3. Prepare OIS results and tested-gene background ####

ois_original <- read_csv_checked(
  ois_results_file,
  c("gene", "log2FoldChange", "padj")
)
ois_background <- read_csv_checked(ois_background_file, "gene")

ois_original$gene_symbol <- clean_symbol(ois_original$gene)
ois_background$gene_symbol <- clean_symbol(ois_background$gene)
ois_background_symbols <- sort(unique(na.omit(ois_background$gene_symbol)))

# A gene must have usable statistics to have had an opportunity to satisfy the
# DEG definition. The explicit OIS background provides an additional guard
# against accidentally using a DEG-only results table as the universe.
ois_eligible <- ois_original[
  !is.na(ois_original$gene_symbol) &
    !is.na(ois_original$padj) &
    !is.na(ois_original$log2FoldChange) &
    ois_original$gene_symbol %in% ois_background_symbols,
  ,
  drop = FALSE
]

if (nrow(ois_eligible) == 0) {
  stop("No eligible OIS genes remain after applying the tested-gene background.")
}

ois_eligible$regulation <- classify_regulation(
  ois_eligible$padj,
  ois_eligible$log2FoldChange
)
ois_resolved <- resolve_symbol_duplicates(ois_eligible, "gene")
ois_prepared <- ois_resolved$data

ois_background_missing_results <- setdiff(
  ois_background_symbols,
  unique(ois_eligible$gene_symbol)
)
ois_results_missing_background <- setdiff(
  unique(na.omit(ois_original$gene_symbol)),
  ois_background_symbols
)


#### 4. Find adjusted TCGA DEG result files ####

cancer_files <- sort(list.files(
  cancer_input_dir,
  pattern = cancer_file_pattern,
  full.names = TRUE
))

if (length(cancer_files) == 0) {
  stop("No adjusted TCGA DEG result files found in: ", cancer_input_dir)
}

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)


#### 5. Analyze one cancer ####

analyze_one_cancer <- function(cancer_file) {
  cancer_project <- sub(
    "_DGE_Raw_Results_Tumor_vs_Normal_Adjusted\\.csv$",
    "",
    basename(cancer_file)
  )
  cancer_output_dir <- file.path(output_root, cancer_project)

  cancer_original <- read_csv_checked(
    cancer_file,
    c("gene_id", "gene_name", "gene_type", "log2FoldChange", "padj")
  )
  cancer_original$gene_symbol <- clean_symbol(cancer_original$gene_name)
  cancer_original$ensembl_id <- sub(
    "\\.[0-9]+$", "", as.character(cancer_original$gene_id)
  )

  cancer_eligible <- cancer_original[
    cancer_original$gene_type == "protein_coding" &
      !is.na(cancer_original$gene_symbol) &
      !is.na(cancer_original$padj) &
      !is.na(cancer_original$log2FoldChange),
    ,
    drop = FALSE
  ]

  if (nrow(cancer_eligible) == 0) {
    stop("No eligible protein-coding cancer genes with usable DEG statistics.")
  }

  cancer_eligible$regulation <- classify_regulation(
    cancer_eligible$padj,
    cancer_eligible$log2FoldChange
  )
  cancer_resolved <- resolve_symbol_duplicates(cancer_eligible, "ensembl_id")
  cancer_prepared <- cancer_resolved$data

  common_tested_genes <- sort(intersect(
    cancer_prepared$gene_symbol,
    ois_prepared$gene_symbol
  ))

  if (length(common_tested_genes) == 0) {
    stop("No common eligible genes with OIS.")
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

  overlap_sets <- list(
    "UP-UP" = intersect(cancer_up, ois_up),
    "UP-DOWN" = intersect(cancer_up, ois_down),
    "DOWN-UP" = intersect(cancer_down, ois_up),
    "DOWN-DOWN" = intersect(cancer_down, ois_down)
  )

  statistics <- rbind(
    run_fisher_overlap(cancer_up, ois_up, common_tested_genes, "UP-UP"),
    run_fisher_overlap(cancer_up, ois_down, common_tested_genes, "UP-DOWN"),
    run_fisher_overlap(cancer_down, ois_up, common_tested_genes, "DOWN-UP"),
    run_fisher_overlap(cancer_down, ois_down, common_tested_genes, "DOWN-DOWN")
  )
  statistics$cancer_project <- cancer_project
  statistics$padj_within_cancer <- p.adjust(statistics$p_value, method = "BH")
  statistics <- statistics[, c(
    "cancer_project", "overlap_class", "universe_n", "cancer_set_n",
    "ois_set_n", "overlap_n", "cancer_only_n", "ois_only_n", "neither_n",
    "expected_overlap", "fold_enrichment", "odds_ratio",
    "odds_ratio_ci_lower", "odds_ratio_ci_upper", "p_value",
    "padj_within_cancer", "alternative"
  )]

  overlap_table <- do.call(rbind, lapply(class_order, function(class_name) {
    symbols <- sort(overlap_sets[[class_name]])
    if (length(symbols) == 0) {
      return(NULL)
    }

    cancer_rows <- cancer_prepared[
      match(symbols, cancer_prepared$gene_symbol),
      c(
        "gene_symbol", "ensembl_id", "log2FoldChange", "padj", "regulation",
        "duplicate_n", "duplicate_direction_conflict"
      )
    ]
    ois_rows <- ois_prepared[
      match(symbols, ois_prepared$gene_symbol),
      c(
        "gene_symbol", "log2FoldChange", "padj", "regulation",
        "duplicate_n", "duplicate_direction_conflict"
      )
    ]

    data.frame(
      gene_symbol = symbols,
      overlap_class = class_name,
      ensembl_id_CANCER = cancer_rows$ensembl_id,
      log2FoldChange_CANCER = cancer_rows$log2FoldChange,
      padj_CANCER = cancer_rows$padj,
      regulation_CANCER = cancer_rows$regulation,
      duplicate_n_CANCER = cancer_rows$duplicate_n,
      duplicate_direction_conflict_CANCER =
        cancer_rows$duplicate_direction_conflict,
      log2FoldChange_OIS = ois_rows$log2FoldChange,
      padj_OIS = ois_rows$padj,
      regulation_OIS = ois_rows$regulation,
      duplicate_n_OIS = ois_rows$duplicate_n,
      duplicate_direction_conflict_OIS = ois_rows$duplicate_direction_conflict,
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(overlap_table)) {
    overlap_table <- data.frame(
      gene_symbol = character(),
      overlap_class = character(),
      stringsAsFactors = FALSE
    )
  }

  dir.create(cancer_output_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(
    data.frame(gene_symbol = common_tested_genes),
    file.path(cancer_output_dir, "common_tested_genes.csv")
  )
  write_csv(
    overlap_table,
    file.path(cancer_output_dir, "overlap_all_classes.csv")
  )

  for (class_name in class_order) {
    write_csv(
      overlap_table[overlap_table$overlap_class == class_name, , drop = FALSE],
      file.path(cancer_output_dir, paste0("overlap_", class_name, ".csv"))
    )
  }

  if (!is.null(cancer_resolved$audit)) {
    write_csv(
      cancer_resolved$audit,
      file.path(cancer_output_dir, "CANCER_symbol_duplicates_audit.csv")
    )
  }

  list(
    cancer_project = cancer_project,
    output_dir = cancer_output_dir,
    statistics = statistics,
    universe_n = length(common_tested_genes),
    cancer_eligible_n = nrow(cancer_prepared),
    cancer_duplicate_symbols_n = if (is.null(cancer_resolved$audit)) {
      0L
    } else {
      length(unique(cancer_resolved$audit$gene_symbol))
    }
  )
}


#### 6. Run all cancers ####

analysis_results <- vector("list", length(cancer_files))
run_status <- vector("list", length(cancer_files))

for (i in seq_along(cancer_files)) {
  cancer_file <- cancer_files[[i]]
  cancer_project <- sub(
    "_DGE_Raw_Results_Tumor_vs_Normal_Adjusted\\.csv$",
    "",
    basename(cancer_file)
  )
  message("[", i, "/", length(cancer_files), "] Running ", cancer_project)

  result <- tryCatch(
    analyze_one_cancer(cancer_file),
    error = function(error) error
  )

  if (inherits(result, "error")) {
    warning("Analysis failed for ", cancer_project, ": ", conditionMessage(result))
    run_status[[i]] <- data.frame(
      cancer_project = cancer_project,
      input_file = cancer_file,
      status = "FAILED",
      universe_n = NA_integer_,
      cancer_eligible_n = NA_integer_,
      cancer_duplicate_symbols_n = NA_integer_,
      error_message = conditionMessage(result),
      stringsAsFactors = FALSE
    )
  } else {
    analysis_results[[i]] <- result
    run_status[[i]] <- data.frame(
      cancer_project = cancer_project,
      input_file = cancer_file,
      status = "COMPLETED",
      universe_n = result$universe_n,
      cancer_eligible_n = result$cancer_eligible_n,
      cancer_duplicate_symbols_n = result$cancer_duplicate_symbols_n,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
  }
}

completed_results <- Filter(Negate(is.null), analysis_results)

if (length(completed_results) == 0) {
  stop("All cancer overlap analyses failed.")
}


#### 7. Global multiple-testing correction and outputs ####

all_statistics <- do.call(rbind, lapply(completed_results, `[[`, "statistics"))
all_statistics$padj_global <- p.adjust(all_statistics$p_value, method = "BH")
all_statistics$significant_global_0.05 <- all_statistics$padj_global < 0.05

# Write the final statistics, including global BH values, both as one combined
# table and back into each cancer-specific directory.
write_csv(
  all_statistics,
  file.path(output_root, "all_cancers_overlap_statistics.csv")
)

for (result in completed_results) {
  cancer_statistics <- all_statistics[
    all_statistics$cancer_project == result$cancer_project,
    ,
    drop = FALSE
  ]
  write_csv(
    cancer_statistics,
    file.path(result$output_dir, "overlap_statistics.csv")
  )
}

run_status <- do.call(rbind, run_status)
write_csv(run_status, file.path(output_root, "batch_run_status.csv"))

ois_audit <- data.frame(
  metric = c(
    "OIS results rows",
    "OIS unique result symbols",
    "OIS background symbols",
    "OIS eligible symbols",
    "OIS background symbols without eligible results",
    "OIS result symbols absent from background",
    "OIS duplicate symbols"
  ),
  value = c(
    nrow(ois_original),
    length(unique(na.omit(ois_original$gene_symbol))),
    length(ois_background_symbols),
    nrow(ois_prepared),
    length(ois_background_missing_results),
    length(ois_results_missing_background),
    if (is.null(ois_resolved$audit)) {
      0L
    } else {
      length(unique(ois_resolved$audit$gene_symbol))
    }
  ),
  stringsAsFactors = FALSE
)
write_csv(ois_audit, file.path(output_root, "OIS_input_audit.csv"))

if (!is.null(ois_resolved$audit)) {
  write_csv(
    ois_resolved$audit,
    file.path(output_root, "OIS_symbol_duplicates_audit.csv")
  )
}

writeLines(
  c(
    paste("Cancer DEG adjusted-p threshold:", padj_threshold),
    paste("OIS DEG adjusted-p threshold:", padj_threshold),
    paste("Absolute log2-fold-change threshold:", lfc_threshold),
    "Cancer input: adjusted-design, raw (unshrunken) DESeq2 result tables",
    "Cancer genes: protein-coding genes with non-missing padj and log2FoldChange",
    "OIS genes: supplied background genes with non-missing padj and log2FoldChange",
    "Identifier: normalized uppercase gene symbol",
    "Test: one-sided Fisher exact test for enrichment (alternative = greater)",
    "Multiplicity: BH within each cancer and globally across all completed tests"
  ),
  file.path(output_root, "analysis_parameters.txt")
)

cat("\nOverlap analysis complete.\n")
print(run_status[, c("cancer_project", "status", "universe_n")], row.names = FALSE)
cat("\nResults saved in:", output_root, "\n")

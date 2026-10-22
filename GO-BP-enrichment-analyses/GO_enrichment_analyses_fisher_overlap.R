#### GO BP enrichment of directional TCGA-OIS overlap gene sets ####

# Run from the repository root after TCGA_OIS_overlaps_fisher.R:
# Rscript GO-BP-enrichment-analyses/GO_enrichment_analyses_fisher_overlap.R


#### 1. Configuration ####

input_root <- file.path(
  "results", "tables", "overlapping-analyses", "TCGA_OIS_overlap_fisher"
)
table_output_root <- file.path(
  "results", "tables", "GO-BP-enrichment-analyses",
  "TCGA_OIS_overlap_fisher"
)
figure_output_root <- file.path(
  "results", "figures", "GO-BP-enrichment-analyses",
  "TCGA_OIS_overlap_fisher"
)

overlap_classes <- c("UP-UP", "UP-DOWN", "DOWN-UP", "DOWN-DOWN")
bh_cutoff <- 0.05
minimum_input_genes <- 5L
minimum_go_set_size <- 10L
maximum_go_set_size <- 500L
top_terms_to_plot <- 14L


#### 2. Dependency and input checks ####

required_packages <- c(
  "AnnotationDbi", "clusterProfiler", "ggplot2", "org.Hs.eg.db"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install and snapshot them in the project environment before running."
  )
}

if (!dir.exists(input_root)) {
  stop(
    "Overlap result directory does not exist: ", input_root,
    ". Run overlapping-analyses/TCGA_OIS_overlaps_fisher.R first."
  )
}

dir.create(table_output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(
  file.path(table_output_root, "audit"),
  recursive = TRUE,
  showWarnings = FALSE
)


#### 3. Helpers ####

read_gene_symbols <- function(path) {
  if (!file.exists(path)) {
    stop("Required input file does not exist: ", path)
  }

  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"gene_symbol" %in% names(x)) {
    stop("Missing gene_symbol column in: ", path)
  }

  symbols <- toupper(trimws(as.character(x$gene_symbol)))
  sort(unique(symbols[!is.na(symbols) & symbols != ""]))
}


write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "")
}


ratio_to_numeric <- function(x) {
  pieces <- strsplit(as.character(x), "/", fixed = TRUE)
  vapply(
    pieces,
    function(piece) {
      if (length(piece) != 2L) {
        return(NA_real_)
      }
      as.numeric(piece[[1]]) / as.numeric(piece[[2]])
    },
    numeric(1)
  )
}


ratio_denominator <- function(x) {
  pieces <- strsplit(as.character(x), "/", fixed = TRUE)
  vapply(
    pieces,
    function(piece) {
      if (length(piece) != 2L) {
        return(NA_real_)
      }
      as.numeric(piece[[2]])
    },
    numeric(1)
  )
}


empty_go_results <- function() {
  data.frame(
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p_adjust_within_list = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = integer(),
    GeneRatio_numeric = numeric(),
    BgRatio_numeric = numeric(),
    expected_gene_count = numeric(),
    FoldEnrichment = numeric(),
    log2_FoldEnrichment = numeric(),
    cancer_project = character(),
    overlap_class = character(),
    input_genes_n = integer(),
    mapped_input_genes_n = integer(),
    background_genes_n = integer(),
    mapped_background_genes_n = integer(),
    stringsAsFactors = FALSE
  )
}


combine_nonempty <- function(x, empty_value = data.frame()) {
  keep <- vapply(x, function(value) !is.null(value) && nrow(value) > 0, logical(1))
  if (!any(keep)) {
    return(empty_value)
  }
  do.call(rbind, x[keep])
}


#### 4. Discover cancer inputs ####

cancer_dirs <- sort(list.dirs(input_root, recursive = FALSE, full.names = TRUE))
cancer_dirs <- cancer_dirs[grepl("^TCGA-[A-Za-z0-9-]+$", basename(cancer_dirs))]

if (length(cancer_dirs) == 0) {
  stop("No TCGA cancer directories found under: ", input_root)
}

all_results <- list()
mapping_qc <- list()
run_log <- list()
result_index <- 1L
log_index <- 1L


#### 5. Run every cancer and directional overlap class ####

for (cancer_dir in cancer_dirs) {
  cancer_project <- basename(cancer_dir)
  message("Preparing ", cancer_project)

  cancer_table_dir <- file.path(table_output_root, "tables", cancer_project)
  cancer_figure_dir <- file.path(figure_output_root, cancer_project)
  cancer_audit_dir <- file.path(table_output_root, "audit", cancer_project)
  dir.create(cancer_table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cancer_figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cancer_audit_dir, recursive = TRUE, showWarnings = FALSE)

  background <- tryCatch(
    read_gene_symbols(file.path(cancer_dir, "common_tested_genes.csv")),
    error = function(error) error
  )

  if (inherits(background, "error")) {
    for (overlap_class in overlap_classes) {
      run_log[[log_index]] <- data.frame(
        cancer_project = cancer_project,
        overlap_class = overlap_class,
        input_genes_n = NA_integer_,
        mapped_input_genes_n = NA_integer_,
        background_genes_n = NA_integer_,
        mapped_background_genes_n = NA_integer_,
        tested_terms_n = 0L,
        significant_within_list_n = 0L,
        status = "FAILED",
        error_message = conditionMessage(background),
        stringsAsFactors = FALSE
      )
      log_index <- log_index + 1L
    }
    next
  }

  # Determine which background symbols have at least one GO BP annotation.
  background_annotations <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = background,
    columns = c("GO", "ONTOLOGY"),
    keytype = "SYMBOL"
  ))
  background_annotations <- background_annotations[
    !is.na(background_annotations$GO) &
      background_annotations$ONTOLOGY == "BP",
    ,
    drop = FALSE
  ]
  mapped_background <- sort(unique(background_annotations$SYMBOL))

  background_audit <- data.frame(
    gene_symbol = background,
    has_GO_BP_annotation = background %in% mapped_background,
    stringsAsFactors = FALSE
  )
  write_csv(
    background_audit,
    file.path(cancer_audit_dir, "background_GO_BP_mapping_audit.csv")
  )

  for (overlap_class in overlap_classes) {
    message("  Running ", overlap_class)
    overlap_file <- file.path(
      cancer_dir,
      paste0("overlap_", overlap_class, ".csv")
    )

    genes <- tryCatch(
      read_gene_symbols(overlap_file),
      error = function(error) error
    )

    if (inherits(genes, "error")) {
      run_log[[log_index]] <- data.frame(
        cancer_project = cancer_project,
        overlap_class = overlap_class,
        input_genes_n = NA_integer_,
        mapped_input_genes_n = NA_integer_,
        background_genes_n = length(background),
        mapped_background_genes_n = length(mapped_background),
        tested_terms_n = 0L,
        significant_within_list_n = 0L,
        status = "FAILED",
        error_message = conditionMessage(genes),
        stringsAsFactors = FALSE
      )
      log_index <- log_index + 1L
      next
    }

    genes_outside_background <- setdiff(genes, background)
    genes <- intersect(genes, background)
    mapped_genes <- intersect(genes, mapped_background)

    foreground_audit <- data.frame(
      gene_symbol = sort(unique(c(genes, genes_outside_background))),
      in_cancer_background = sort(unique(c(genes, genes_outside_background))) %in%
        background,
      has_GO_BP_annotation = sort(unique(c(genes, genes_outside_background))) %in%
        mapped_background,
      stringsAsFactors = FALSE
    )
    write_csv(
      foreground_audit,
      file.path(
        cancer_audit_dir,
        paste0(overlap_class, "_input_gene_audit.csv")
      )
    )

    mapping_qc[[log_index]] <- data.frame(
      cancer_project = cancer_project,
      overlap_class = overlap_class,
      input_file_genes_n = length(genes) + length(genes_outside_background),
      genes_outside_background_n = length(genes_outside_background),
      input_genes_n = length(genes),
      mapped_input_genes_n = length(mapped_genes),
      input_GO_BP_mapping_rate = if (length(genes) > 0) {
        length(mapped_genes) / length(genes)
      } else {
        NA_real_
      },
      background_genes_n = length(background),
      mapped_background_genes_n = length(mapped_background),
      background_GO_BP_mapping_rate = length(mapped_background) /
        length(background),
      stringsAsFactors = FALSE
    )

    if (length(genes) < minimum_input_genes) {
      run_log[[log_index]] <- data.frame(
        cancer_project = cancer_project,
        overlap_class = overlap_class,
        input_genes_n = length(genes),
        mapped_input_genes_n = length(mapped_genes),
        background_genes_n = length(background),
        mapped_background_genes_n = length(mapped_background),
        tested_terms_n = 0L,
        significant_within_list_n = 0L,
        status = paste0("SKIPPED: fewer than ", minimum_input_genes, " genes"),
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )
      log_index <- log_index + 1L
      next
    }

    go_result <- tryCatch(
      clusterProfiler::enrichGO(
        gene = genes,
        universe = background,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        qvalueCutoff = 1,
        minGSSize = minimum_go_set_size,
        maxGSSize = maximum_go_set_size,
        readable = FALSE
      ),
      error = function(error) error
    )

    if (inherits(go_result, "error")) {
      run_log[[log_index]] <- data.frame(
        cancer_project = cancer_project,
        overlap_class = overlap_class,
        input_genes_n = length(genes),
        mapped_input_genes_n = length(mapped_genes),
        background_genes_n = length(background),
        mapped_background_genes_n = length(mapped_background),
        tested_terms_n = 0L,
        significant_within_list_n = 0L,
        status = "FAILED",
        error_message = conditionMessage(go_result),
        stringsAsFactors = FALSE
      )
      log_index <- log_index + 1L
      next
    }

    results <- as.data.frame(go_result)

    if (nrow(results) == 0) {
      run_log[[log_index]] <- data.frame(
        cancer_project = cancer_project,
        overlap_class = overlap_class,
        input_genes_n = length(genes),
        mapped_input_genes_n = length(mapped_genes),
        background_genes_n = length(background),
        mapped_background_genes_n = length(mapped_background),
        tested_terms_n = 0L,
        significant_within_list_n = 0L,
        status = "NO GO BP TERMS",
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )
      log_index <- log_index + 1L
      next
    }

    results$GeneRatio_numeric <- ratio_to_numeric(results$GeneRatio)
    results$BgRatio_numeric <- ratio_to_numeric(results$BgRatio)
    results$expected_gene_count <- ratio_denominator(results$GeneRatio) *
      results$BgRatio_numeric
    results$FoldEnrichment <- results$GeneRatio_numeric /
      results$BgRatio_numeric
    results$log2_FoldEnrichment <- log2(results$FoldEnrichment)
    results$p_adjust_within_list <- results$p.adjust
    results$cancer_project <- cancer_project
    results$overlap_class <- overlap_class
    results$input_genes_n <- length(genes)
    results$mapped_input_genes_n <- length(mapped_genes)
    results$background_genes_n <- length(background)
    results$mapped_background_genes_n <- length(mapped_background)
    results$p.adjust <- NULL

    all_results[[result_index]] <- results
    result_index <- result_index + 1L

    run_log[[log_index]] <- data.frame(
      cancer_project = cancer_project,
      overlap_class = overlap_class,
      input_genes_n = length(genes),
      mapped_input_genes_n = length(mapped_genes),
      background_genes_n = length(background),
      mapped_background_genes_n = length(mapped_background),
      tested_terms_n = nrow(results),
      significant_within_list_n = sum(
        results$p_adjust_within_list < bh_cutoff,
        na.rm = TRUE
      ),
      status = "COMPLETED",
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    log_index <- log_index + 1L
  }
}


#### 6. Global BH correction and final tables ####

combined_results <- combine_nonempty(all_results, empty_go_results())
combined_mapping_qc <- combine_nonempty(mapping_qc)
combined_run_log <- combine_nonempty(run_log)

if (nrow(combined_results) > 0) {
  combined_results$p_adjust_global <- p.adjust(
    combined_results$pvalue,
    method = "BH"
  )
  combined_results$significant_within_list_0.05 <-
    combined_results$p_adjust_within_list < bh_cutoff
  combined_results$significant_global_0.05 <-
    combined_results$p_adjust_global < bh_cutoff
} else {
  combined_results$p_adjust_global <- numeric()
  combined_results$significant_within_list_0.05 <- logical()
  combined_results$significant_global_0.05 <- logical()
}

write_csv(
  combined_results,
  file.path(table_output_root, "GO_BP_all_results.csv")
)
write_csv(
  combined_results[
    combined_results$significant_within_list_0.05,
    ,
    drop = FALSE
  ],
  file.path(table_output_root, "GO_BP_significant_within_list_BH_0.05.csv")
)
write_csv(
  combined_results[
    combined_results$significant_global_0.05,
    ,
    drop = FALSE
  ],
  file.path(table_output_root, "GO_BP_significant_global_BH_0.05.csv")
)
write_csv(
  combined_mapping_qc,
  file.path(table_output_root, "GO_BP_mapping_QC.csv")
)

# Add the final global-significance count to each run-log row.
if (nrow(combined_run_log) > 0 && nrow(combined_results) > 0) {
  global_counts <- aggregate(
    significant_global_0.05 ~ cancer_project + overlap_class,
    data = combined_results,
    FUN = sum
  )
  names(global_counts)[[3]] <- "significant_global_n"
  combined_run_log <- merge(
    combined_run_log,
    global_counts,
    by = c("cancer_project", "overlap_class"),
    all.x = TRUE,
    sort = FALSE
  )
  combined_run_log$significant_global_n[
    is.na(combined_run_log$significant_global_n)
  ] <- 0L
} else if (nrow(combined_run_log) > 0) {
  combined_run_log$significant_global_n <- 0L
}
write_csv(combined_run_log, file.path(table_output_root, "GO_BP_run_log.csv"))


#### 7. Cancer/class tables and globally adjusted plots ####

for (cancer_project in unique(combined_results$cancer_project)) {
  cancer_table_dir <- file.path(table_output_root, "tables", cancer_project)
  cancer_figure_dir <- file.path(figure_output_root, cancer_project)

  for (overlap_class in overlap_classes) {
    results <- combined_results[
      combined_results$cancer_project == cancer_project &
        combined_results$overlap_class == overlap_class,
      ,
      drop = FALSE
    ]

    if (nrow(results) == 0) {
      next
    }

    write_csv(
      results,
      file.path(cancer_table_dir, paste0(overlap_class, "_GO_BP_all.csv"))
    )
    write_csv(
      results[results$significant_within_list_0.05, , drop = FALSE],
      file.path(
        cancer_table_dir,
        paste0(overlap_class, "_GO_BP_significant_within_list_BH_0.05.csv")
      )
    )
    globally_significant <- results[
      results$significant_global_0.05,
      ,
      drop = FALSE
    ]
    write_csv(
      globally_significant,
      file.path(
        cancer_table_dir,
        paste0(overlap_class, "_GO_BP_significant_global_BH_0.05.csv")
      )
    )

    if (nrow(globally_significant) == 0) {
      next
    }

    plot_data <- head(
      globally_significant[order(globally_significant$p_adjust_global), ],
      top_terms_to_plot
    )
    plot_data$Description <- factor(
      plot_data$Description,
      levels = rev(plot_data$Description)
    )
    plot_data$minus_log10_global_BH <- -log10(pmax(
      plot_data$p_adjust_global,
      .Machine$double.xmin
    ))

    plot_object <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = FoldEnrichment,
        y = Description,
        size = Count,
        colour = minus_log10_global_BH
      )
    ) +
      ggplot2::geom_point() +
      ggplot2::scale_colour_gradient(
        low = "#2166AC",
        high = "#B2182B",
        name = expression(-log[10]("global BH"))
      ) +
      ggplot2::labs(
        title = paste(cancer_project, overlap_class),
        subtitle = "GO Biological Process enrichment",
        x = "Fold enrichment",
        y = NULL,
        size = "Overlap genes"
      ) +
      ggplot2::theme_bw(base_size = 11)

    ggplot2::ggsave(
      file.path(cancer_figure_dir, paste0(overlap_class, "_GO_BP.png")),
      plot_object,
      width = 10,
      height = 7,
      dpi = 300,
      bg = "white"
    )
    ggplot2::ggsave(
      file.path(cancer_figure_dir, paste0(overlap_class, "_GO_BP.pdf")),
      plot_object,
      width = 10,
      height = 7,
      bg = "white"
    )
  }
}


#### 8. Reproducibility records ####

settings <- data.frame(
  parameter = c(
    "input_root",
    "table_output_root",
    "figure_output_root",
    "ontology",
    "key_type",
    "p_adjust_method",
    "within_list_and_global_BH_cutoff",
    "minimum_input_genes",
    "minimum_GO_set_size",
    "maximum_GO_set_size",
    "top_terms_to_plot",
    "plot_significance_definition"
  ),
  value = c(
    input_root,
    table_output_root,
    figure_output_root,
    "BP",
    "SYMBOL",
    "BH",
    bh_cutoff,
    minimum_input_genes,
    minimum_go_set_size,
    maximum_go_set_size,
    top_terms_to_plot,
    "global BH < 0.05"
  ),
  stringsAsFactors = FALSE
)
write_csv(settings, file.path(table_output_root, "GO_BP_analysis_settings.csv"))

session_lines <- capture.output(sessionInfo())
writeLines(session_lines, file.path(table_output_root, "GO_BP_session_info.txt"))

cat("\nGO BP enrichment complete.\n")
cat("Tables:", table_output_root, "\n")
cat("Figures:", figure_output_root, "\n")

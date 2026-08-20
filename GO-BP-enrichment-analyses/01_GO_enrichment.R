#### Step 1: GO Biological Process over-representation analysis ####

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
  file.path(getwd(), "01_GO_enrichment.R")
}
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))

local_library <- file.path(script_dir, ".Rlib")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))

source(file.path(script_dir, "pipeline_utils.R"))
source(file.path(script_dir, "config.R"))

require_packages(c(
  "AnnotationDbi", "clusterProfiler", "dplyr", "ggplot2", "org.Hs.eg.db"
))

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
table_root <- file.path(output_root, "per_cancer_GO")
plot_root <- file.path(output_root, "per_cancer_plots")
dir.create(table_root, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_root, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(input_root)) {
  stop(
    "Input folder does not exist: ", input_root, "\n",
    "Set GO_PIPELINE_INPUT_ROOT or place TCGA folders in: ",
    file.path(script_dir, "data")
  )
}

cancer_dirs <- sort(list.dirs(input_root, recursive = FALSE, full.names = TRUE))
cancer_dirs <- cancer_dirs[grepl("^TCGA-", basename(cancer_dirs))]
if (!length(cancer_dirs)) stop("No TCGA-* folders found in: ", input_root)

keep_annotated_symbols <- function(genes) {
  mapped <- suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = genes,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  ))
  names(mapped)[!is.na(mapped)]
}

all_results <- list()
run_log <- list()
result_index <- 1L
log_index <- 1L

for (cancer_dir in cancer_dirs) {
  cancer_folder <- basename(cancer_dir)
  cancer <- sub("^TCGA-", "", cancer_folder)
  cancer_table_dir <- file.path(table_root, cancer_folder)
  cancer_plot_dir <- file.path(plot_root, cancer_folder)
  dir.create(cancer_table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cancer_plot_dir, recursive = TRUE, showWarnings = FALSE)

  background_raw <- read_gene_symbols(
    file.path(cancer_dir, "common_tested_genes.csv")
  )
  background <- keep_annotated_symbols(background_raw)

  for (overlap_class in classes) {
    message("GO enrichment: ", cancer, " / ", overlap_class)
    overlap_file <- file.path(cancer_dir, paste0("overlap_", overlap_class, ".csv"))
    genes_raw <- read_gene_symbols(overlap_file)
    genes_in_background <- intersect(genes_raw, background_raw)
    genes <- intersect(keep_annotated_symbols(genes_in_background), background)

    common_log <- data.frame(
      cancer = cancer,
      cancer_folder = cancer_folder,
      overlap_class = overlap_class,
      input_genes_raw = length(genes_raw),
      input_genes_in_background = length(genes_in_background),
      input_genes_mapped = length(genes),
      background_genes_raw = length(background_raw),
      background_genes_mapped = length(background),
      stringsAsFactors = FALSE
    )

    if (length(genes) < minimum_input_genes) {
      run_log[[log_index]] <- dplyr::mutate(
        common_log,
        tested_terms = 0L,
        significant_terms = 0L,
        status = paste0("SKIPPED: fewer than ", minimum_input_genes, " mapped genes")
      )
      log_index <- log_index + 1L
      next
    }

    go <- tryCatch(
      clusterProfiler::enrichGO(
        gene = genes,
        universe = background,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        qvalueCutoff = 1,
        minGSSize = minimum_gene_set_size,
        maxGSSize = maximum_gene_set_size,
        readable = FALSE
      ),
      error = function(error) error
    )

    if (inherits(go, "error")) {
      run_log[[log_index]] <- dplyr::mutate(
        common_log,
        tested_terms = 0L,
        significant_terms = 0L,
        status = paste0("ERROR: ", conditionMessage(go))
      )
      log_index <- log_index + 1L
      next
    }

    results <- as.data.frame(go)
    if (!nrow(results)) {
      run_log[[log_index]] <- dplyr::mutate(
        common_log,
        tested_terms = 0L,
        significant_terms = 0L,
        status = "NO GO TERMS"
      )
      log_index <- log_index + 1L
      next
    }

    results$GeneRatio_numeric <- parse_ratio(results$GeneRatio)
    results$BgRatio_numeric <- parse_ratio(results$BgRatio)
    results$FoldEnrichment <- results$GeneRatio_numeric / results$BgRatio_numeric
    results$log2_FoldEnrichment <- log2(results$FoldEnrichment)
    results$significant <- !is.na(results$p.adjust) & results$p.adjust < bh_cutoff
    results$cancer <- cancer
    results$cancer_folder <- cancer_folder
    results$overlap_class <- overlap_class
    results$input_genes <- length(genes)
    results$background_genes <- length(background)

    significant <- results[results$significant, , drop = FALSE]
    write_csv_utf8(
      results,
      file.path(cancer_table_dir, paste0(overlap_class, "_GO_BP_all.csv"))
    )
    write_csv_utf8(
      significant,
      file.path(cancer_table_dir, paste0(overlap_class, "_GO_BP_BH_0.05.csv"))
    )

    all_results[[result_index]] <- results
    result_index <- result_index + 1L
    run_log[[log_index]] <- dplyr::mutate(
      common_log,
      tested_terms = nrow(results),
      significant_terms = nrow(significant),
      status = "COMPLETED"
    )
    log_index <- log_index + 1L

    if (nrow(significant)) {
      plot_data <- significant |>
        dplyr::arrange(.data$p.adjust, dplyr::desc(.data$log2_FoldEnrichment), .data$ID) |>
        dplyr::slice_head(n = min(14L, nrow(significant))) |>
        dplyr::arrange(.data$GeneRatio_numeric, .data$p.adjust, .data$ID) |>
        dplyr::mutate(
          Description = factor(
            .data$Description,
            levels = unique(.data$Description)
          )
        )

      plot <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
          x = .data$GeneRatio_numeric,
          y = .data$Description,
          size = .data$Count,
          colour = .data$p.adjust
        )
      ) +
        ggplot2::geom_point() +
        ggplot2::scale_colour_gradient(
          low = "#D7191C", high = "#2C7BB6", trans = "log10",
          name = "BH-adjusted p"
        ) +
        ggplot2::labs(
          title = paste(cancer, overlap_class),
          subtitle = "GO Biological Process enrichment",
          x = "Gene ratio", y = NULL, size = "Gene count"
        ) +
        ggplot2::theme_bw(base_size = 11)

      ggplot2::ggsave(
        file.path(cancer_plot_dir, paste0(overlap_class, "_GO_BP.png")),
        plot, width = 10, height = 7, dpi = 300, bg = "white"
      )
    }
  }
}

combined_results <- bind_nonempty(all_results)
combined_log <- bind_nonempty(run_log)
significant_results <- combined_results |>
  dplyr::filter(!is.na(.data$p.adjust), .data$p.adjust < bh_cutoff)

write_csv_utf8(combined_log, file.path(output_root, "00_GO_enrichment_run_log.csv"))
write_csv_utf8(
  combined_results,
  file.path(output_root, "01_GO_enrichment_all_results.csv")
)
write_csv_utf8(
  significant_results,
  file.path(output_root, "02_GO_significant_terms.csv")
)
capture.output(sessionInfo(), file = file.path(output_root, "00_sessionInfo_step1.txt"))

failed <- grepl("^ERROR:", combined_log$status)
if (any(failed)) {
  stop(
    sum(failed), " cancer/class analyses failed. See: ",
    file.path(output_root, "00_GO_enrichment_run_log.csv")
  )
}

message(
  "Step 1 complete: ", nrow(combined_log), " analyses; ",
  nrow(significant_results), " significant cancer/class/GO rows.\n",
  "Results: ", output_root
)

#### GO BP enrichment analysis ####

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)


#### 1. Settings ####

input_root <- "/Users/danabiruk/Documents/all_TCGA_OIS_overlaps"
output_root <- file.path(input_root, "GO_BP_analysis")

classes <- c("UP-UP", "UP-DOWN", "DOWN-UP", "DOWN-DOWN")
bh_cutoff <- 0.05
top_terms <- 14

dir.create(file.path(output_root, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "plots"), recursive = TRUE, showWarnings = FALSE)


#### 2. Find cancer folders ####

cancer_dirs <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
cancer_dirs <- sort(cancer_dirs[grepl("^TCGA-", basename(cancer_dirs))])


#### 3. Read unique gene symbols ####

read_genes <- function(file) {
  x <- read.csv(file, stringsAsFactors = FALSE)
  genes <- toupper(trimws(x$gene_symbol))
  sort(unique(genes[!is.na(genes) & genes != ""]))
}


#### 4. Empty containers for combined results ####

all_results <- list()
significant_results <- list()
run_log <- list()
index <- 1


#### 5. Run every cancer and every overlap class ####

for (cancer_dir in cancer_dirs) {
  cancer <- basename(cancer_dir)
  
  table_dir <- file.path(output_root, "tables", cancer)
  plot_dir <- file.path(output_root, "plots", cancer)
  
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  background <- read_genes(
    file.path(cancer_dir, "common_tested_genes.csv")
  )
  
  for (overlap_class in classes) {
    message("Running ", cancer, " / ", overlap_class)
    
    overlap_file <- file.path(
      cancer_dir,
      paste0("overlap_", overlap_class, ".csv")
    )
    
    genes <- read_genes(overlap_file)
    genes <- intersect(genes, background)
    
    if (length(genes) < 5) {
      run_log[[index]] <- data.frame(
        cancer = cancer,
        overlap_class = overlap_class,
        input_genes = length(genes),
        background_genes = length(background),
        tested_terms = 0,
        significant_terms = 0,
        status = "SKIPPED: fewer than 5 genes"
      )
      
      index <- index + 1
      next
    }
    
    go <- enrichGO(
      gene = genes,
      universe = background,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      qvalueCutoff = 1,
      minGSSize = 10,
      maxGSSize = 500,
      readable = FALSE
    )
    
    results <- as.data.frame(go)
    
    if (nrow(results) == 0) {
      run_log[[index]] <- data.frame(
        cancer = cancer,
        overlap_class = overlap_class,
        input_genes = length(genes),
        background_genes = length(background),
        tested_terms = 0,
        significant_terms = 0,
        status = "NO GO TERMS"
      )
      
      index <- index + 1
      next
    }
    
    # Convert ratios into numbers.
    results$GeneRatio_numeric <- vapply(
      strsplit(results$GeneRatio, "/"),
      function(x) as.numeric(x[1]) / as.numeric(x[2]),
      numeric(1)
    )
    
    results$BgRatio_numeric <- vapply(
      strsplit(results$BgRatio, "/"),
      function(x) as.numeric(x[1]) / as.numeric(x[2]),
      numeric(1)
    )
    
    results$FoldEnrichment <- (
      results$GeneRatio_numeric / results$BgRatio_numeric
    )
    
    results$log2_FoldEnrichment <- log2(results$FoldEnrichment)
    results$significant <- results$p.adjust < bh_cutoff
    results$cancer <- cancer
    results$overlap_class <- overlap_class
    results$input_genes <- length(genes)
    results$background_genes <- length(background)
    
    significant <- results[results$significant, , drop = FALSE]
    
    write.csv(
      results,
      file.path(table_dir, paste0(overlap_class, "_GO_BP_all.csv")),
      row.names = FALSE
    )
    
    write.csv(
      significant,
      file.path(table_dir, paste0(overlap_class, "_GO_BP_BH_0.05.csv")),
      row.names = FALSE
    )
    
    all_results[[index]] <- results
    significant_results[[index]] <- significant
    
    run_log[[index]] <- data.frame(
      cancer = cancer,
      overlap_class = overlap_class,
      input_genes = length(genes),
      background_genes = length(background),
      tested_terms = nrow(results),
      significant_terms = nrow(significant),
      status = "COMPLETED"
    )
    
    if (nrow(significant) > 0) {
      plot_data <- head(significant[order(significant$p.adjust), ], top_terms)
      plot_data$Description <- factor(
        plot_data$Description,
        levels = rev(plot_data$Description)
      )
      
      p <- ggplot(
        plot_data,
        aes(
          GeneRatio_numeric,
          Description,
          size = Count,
          colour = p.adjust
        )
      ) +
        geom_point() +
        scale_colour_gradient(
          low = "red",
          high = "blue",
          trans = "log10",
          name = "BH adjusted p-value"
        ) +
        labs(
          title = paste(cancer, overlap_class),
          subtitle = "GO Biological Process enrichment",
          x = "Gene ratio",
          y = NULL,
          size = "Gene count"
        ) +
        theme_bw(base_size = 11)
      
      print(p)
      
      ggsave(
        file.path(plot_dir, paste0(overlap_class, "_GO_BP.png")),
        p,
        width = 10,
        height = 7,
        dpi = 300,
        bg = "white"
      )
    }
    
    index <- index + 1
  }
}


#### 6. Save combined tables in the main enrichment folder ####

combine_tables <- function(x) {
  x <- x[vapply(x, function(z) !is.null(z) && nrow(z) > 0, logical(1))]
  if (length(x) == 0) return(data.frame())
  do.call(rbind, x)
}

write.csv(
  combine_tables(all_results),
  file.path(output_root, "GO_BP_ALL_CANCERS_all_results.csv"),
  row.names = FALSE
)

write.csv(
  combine_tables(significant_results),
  file.path(output_root, "GO_BP_ALL_CANCERS_significant_BH_0.05.csv"),
  row.names = FALSE
)

write.csv(
  combine_tables(run_log),
  file.path(output_root, "GO_BP_run_log.csv"),
  row.names = FALSE
)

cat("\nFinished. Results saved in:\n", output_root, "\n")

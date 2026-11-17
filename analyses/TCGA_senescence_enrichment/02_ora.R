#### Over-Representation Analysis (step 02) ####
#
# One-sided hypergeometric test of each gene set in three kinds of gene lists:
# - class: the four overlap classes of each cancer x signature (UP-UP,
#   UP-DOWN, DOWN-UP, DOWN-DOWN; cancer first), universe = genes tested in
#   both (common_tested_genes of the overlap module);
# - cancer: tumor up and down DEGs, universe = genes tested in the tumor;
# - signature: signature up and down DEGs, universe = signature background.
# Class results mix the enrichment of each side (a gene set enriched in tumor-up
# genes and in signature-up genes is enriched in UP-UP even without a joint
# effect), so read them next to the cancer and signature results.
#
# For each gene list and collection: universe restricted to annotated genes,
# gene sets with min_gene_set_size-max_gene_set_size genes in that universe,
# BH over all tested gene sets (including those with no hit). Written: rows
# with p < 0.05 (all BH-significant rows included), and a status table.

source(file.path("analyses", "TCGA_senescence_enrichment", "00_config.R"))

#### Gene lists ####

read_symbols <- function(file) {
  x <- read.csv(file, stringsAsFactors = FALSE)
  unique(toupper(x$gene_symbol[!is.na(x$gene_symbol) & nzchar(x$gene_symbol)]))
}

gene_lists <- function() {
  lists <- list()
  add <- function(kind, cancer_project, signature, direction, genes, universe) {
    lists[[length(lists) + 1L]] <<- list(
      kind = kind, cancer_project = cancer_project, signature = signature,
      direction = direction, genes = genes, universe = universe
    )
  }

  for (signature in signatures) {
    folder <- file.path(overlap_table_dir, signature)
    cancers <- list.dirs(folder, recursive = FALSE, full.names = FALSE)
    for (cancer_project in sort(cancers[grepl("^TCGA-", cancers)])) {
      cancer_folder <- file.path(folder, cancer_project)
      universe <- read_symbols(file.path(cancer_folder, "common_tested_genes.csv"))
      for (class_name in overlap_config$class_order) {
        add("class", cancer_project, signature, class_name,
            read_symbols(file.path(cancer_folder, paste0("overlap_", class_name, ".csv"))),
            universe)
      }
    }
    genes <- read.csv(file.path(folder, paste0(signature, "_signature_genes.csv")), stringsAsFactors = FALSE)
    universe <- read_symbols(file.path(folder, paste0(signature, "_signature_background.csv")))
    for (direction in c("UP", "DOWN")) {
      add("signature", NA_character_, signature, direction,
          unique(genes$gene_symbol[genes$regulation == direction]), universe)
    }
  }

  cancer_files <- list.files(
    file.path(overlap_table_dir, "cancer_gene_tables"), "_genes\\.csv$", full.names = TRUE
  )
  for (file in sort(cancer_files)) {
    genes <- read.csv(file, stringsAsFactors = FALSE)
    cancer_project <- sub("_genes\\.csv$", "", basename(file))
    for (direction in c("UP", "DOWN")) {
      add("cancer", cancer_project, NA_character_, direction,
          unique(genes$gene_symbol[genes$regulation == direction]), unique(genes$gene_symbol))
    }
  }
  lists
}

#### ORA ####

# Hypergeometric test of all gene sets of one collection in one gene list.
ora_one <- function(gene_list, gene_sets) {
  annotated <- intersect(gene_list$universe, gene_sets$gene_symbol)
  input <- intersect(gene_list$genes, annotated)
  status <- data.frame(
    n_input_raw = length(gene_list$genes),
    n_universe_raw = length(gene_list$universe),
    n_input = length(input),
    n_universe = length(annotated),
    n_gene_sets_tested = 0L,
    n_significant = 0L,
    status = "completed",
    stringsAsFactors = FALSE
  )
  if (length(input) < ora_min_input_genes) {
    status$status <- "skipped: too few input genes"
    return(list(table = NULL, status = status))
  }

  members <- gene_sets[gene_sets$gene_symbol %in% annotated, , drop = FALSE]
  set_size <- table(members$gene_set)
  set_size <- set_size[set_size >= min_gene_set_size & set_size <= max_gene_set_size]
  if (!length(set_size)) {
    status$status <- "skipped: no gene set of eligible size"
    return(list(table = NULL, status = status))
  }
  members <- members[members$gene_set %in% names(set_size), , drop = FALSE]
  hit_rows <- members$gene_symbol %in% input
  overlap <- table(factor(members$gene_set[hit_rows], levels = names(set_size)))

  k <- as.integer(overlap)
  K <- as.integer(set_size)
  n <- length(input)
  N <- length(annotated)
  p_value <- stats::phyper(k - 1, K, N - K, n, lower.tail = FALSE)
  table <- data.frame(
    gene_set = names(set_size),
    set_size = K,
    overlap = k,
    expected = K * n / N,
    fold_enrichment = (k / n) / (K / N),
    p_value = p_value,
    padj = stats::p.adjust(p_value, method = "BH"),
    stringsAsFactors = FALSE
  )
  status$n_gene_sets_tested <- nrow(table)
  status$n_significant <- sum(table$padj < fdr_threshold)

  table <- table[table$p_value < 0.05, , drop = FALSE]
  hits <- split(members$gene_symbol[hit_rows], members$gene_set[hit_rows])
  table$overlap_genes <- vapply(
    table$gene_set, function(set) paste(sort(hits[[set]]), collapse = ";"), character(1)
  )
  list(table = table[order(table$p_value), , drop = FALSE], status = status)
}

#### Run ####

lists <- gene_lists()
message(length(lists), " gene lists")
list_info <- function(gene_list) {
  data.frame(
    kind = gene_list$kind, cancer_project = gene_list$cancer_project,
    signature = gene_list$signature, direction = gene_list$direction,
    stringsAsFactors = FALSE
  )
}

statuses <- list()
for (collection in names(gene_set_collections)) {
  message("ORA: ", collection)
  gene_sets <- readRDS(file.path(gene_set_dir, paste0(collection, ".rds")))
  gene_sets <- unique(gene_sets[, c("gene_set", "gene_symbol")])
  results <- lapply(lists, function(gene_list) ora_one(gene_list, gene_sets))

  tables <- dplyr::bind_rows(lapply(seq_along(lists), function(i) {
    if (is.null(results[[i]]$table) || !nrow(results[[i]]$table)) return(NULL)
    cbind(list_info(lists[[i]]), collection = collection, results[[i]]$table, row.names = NULL)
  }))
  write.csv(tables, file.path(results_dir, paste0("02_ORA_", collection, ".csv")), row.names = FALSE)

  statuses[[collection]] <- dplyr::bind_rows(lapply(seq_along(lists), function(i) {
    cbind(list_info(lists[[i]]), collection = collection, results[[i]]$status)
  }))
}

status <- dplyr::bind_rows(statuses)
write.csv(status, file.path(results_dir, "02_ORA_status.csv"), row.names = FALSE)
message("ORA results saved to: ", results_dir)
print(table(status$collection, status$status))

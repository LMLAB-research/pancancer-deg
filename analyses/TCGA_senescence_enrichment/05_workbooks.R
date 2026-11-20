#### Excel Workbooks (step 05) ####
#
# Collects the main enrichment results into Excel workbooks for inspection
# (e.g. which genes drive a gene set, issue #14). Optional: all results are
# also in the CSV tables of steps 02-04.
# - one workbook per cancer: summary (step 04), tumor GSEA (significant gene
#   sets), and per signature the gene sets significant on both sides (step 04)
#   and the significant ORA results of the four classes (step 02);
# - one workbook for the signatures: significant GSEA gene sets per signature.

source(file.path("analyses", "TCGA_senescence_enrichment", "00_config.R"))

workbook_dir <- file.path(results_dir, "workbooks")
dir.create(workbook_dir, recursive = TRUE, showWarnings = FALSE)

read_all <- function(prefix) {
  dplyr::bind_rows(lapply(names(gene_set_collections), function(collection) {
    read.csv(file.path(results_dir, paste0(prefix, collection, ".csv")), stringsAsFactors = FALSE)
  }))
}

# Excel cells hold at most 32,767 characters.
excel_safe <- function(x) {
  for (column in names(x)) {
    if (is.character(x[[column]])) {
      long <- !is.na(x[[column]]) & nchar(x[[column]]) > 32000
      x[[column]][long] <- paste0(substr(x[[column]][long], 1, 32000), "...[truncated; see CSV]")
    }
  }
  x
}

add_sheet <- function(workbook, name, x) {
  if (!nrow(x)) {
    x <- data.frame(note = "no significant results")
  }
  openxlsx::addWorksheet(workbook, name)
  openxlsx::writeData(workbook, name, excel_safe(x), withFilter = TRUE)
  openxlsx::freezePane(workbook, name, firstRow = TRUE)
  openxlsx::setColWidths(workbook, name, cols = seq_len(ncol(x)), widths = "auto")
}

significant <- function(padj) !is.na(padj) & padj < fdr_threshold

summary <- read.csv(file.path(results_dir, "04_compare_summary.csv"), stringsAsFactors = FALSE)
compare <- read_all("04_compare_")
gsea <- read_all("03_GSEA_")
ora <- read_all("02_ORA_")

for (cancer_project in sort(unique(summary$cancer_project))) {
  workbook <- openxlsx::createWorkbook()
  add_sheet(workbook, "summary", summary[summary$cancer_project == cancer_project, ])
  tumor <- gsea[gsea$kind == "cancer" & gsea$cancer_project %in% cancer_project & significant(gsea$padj), ]
  add_sheet(workbook, "tumor_GSEA", tumor[order(tumor$collection, -abs(tumor$NES)), ])
  for (signature in signatures) {
    both <- compare[
      compare$cancer_project == cancer_project & compare$signature == signature &
        compare$direction %in% c("same", "opposite"),
    ]
    add_sheet(workbook, paste0(signature, "_GSEA"), both[order(both$collection, both$direction, -abs(both$NES_cancer)), ])
    classes <- ora[
      ora$kind == "class" & ora$cancer_project %in% cancer_project &
        ora$signature %in% signature & significant(ora$padj),
    ]
    add_sheet(workbook, paste0(signature, "_ORA"), classes[order(classes$direction, classes$collection, classes$padj), ])
  }
  openxlsx::saveWorkbook(
    workbook, file.path(workbook_dir, paste0(cancer_project, "_enrichment.xlsx")), overwrite = TRUE
  )
  message("Saved workbook for ", cancer_project)
}

workbook <- openxlsx::createWorkbook()
for (signature in signatures) {
  x <- gsea[gsea$kind == "signature" & gsea$signature %in% signature & significant(gsea$padj), ]
  add_sheet(workbook, paste0(signature, "_GSEA"), x[order(x$collection, -abs(x$NES)), ])
}
openxlsx::saveWorkbook(workbook, file.path(workbook_dir, "signatures_enrichment.xlsx"), overwrite = TRUE)
message("Workbooks saved to: ", workbook_dir)

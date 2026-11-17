# Input preparation retained from TCGA_senescence_overlaps.R (2026-10-18).
classify_regulation <- function(padj, lfc) {
    ifelse(!is.na(padj) & !is.na(lfc) & padj < padj_threshold & 
        lfc > lfc_threshold, "UP", ifelse(!is.na(padj) & !is.na(lfc) & 
        padj < padj_threshold & lfc < -lfc_threshold, "DOWN", 
        "NS"))
}
validate_numeric <- function(x, column, file) {
    y <- suppressWarnings(as.numeric(x))
    bad <- !is.na(x) & trimws(as.character(x)) != "" & is.na(y)
    if (any(bad)) {
        stop("Non-numeric values in ", column, ": ", file)
    }
    if (any(!is.na(y) & !is.finite(y))) 
        stop("Non-finite values in ", column, ": ", file)
    if (column == "padj" && any(y < 0 | y > 1, na.rm = TRUE)) {
        stop("Adjusted p-values must be between 0 and 1: ", file)
    }
    y
}
resolve_symbol_duplicates <- function(x, source_id_column, symbol_column = "gene_symbol") {
    x <- x[!is.na(x[[symbol_column]]) & x[[symbol_column]] != 
        "", , drop = FALSE]
    if (!nrow(x)) {
        stop("No valid gene symbols remain.")
    }
    groups <- split(x, x[[symbol_column]])
    selected <- list()
    audit <- list()
    for (symbol in names(groups)) {
        group <- groups[[symbol]]
        directions <- unique(group$regulation[group$regulation %in% 
            c("UP", "DOWN")])
        conflict <- length(directions) > 1
        row_order <- order(is.na(group$padj), group$padj, -abs(group$log2FoldChange), 
            group[[source_id_column]], na.last = TRUE)
        selected_row <- row_order[[1]]
        representative <- group[selected_row, , drop = FALSE]
        status <- if (nrow(group) == 1) {
            "UNIQUE"
        }
        else if (conflict) {
            "CONFLICT_OPPOSITE_DIRECTIONS"
        }
        else {
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
        group$selected_representative <- seq_len(nrow(group)) == 
            selected_row
        selected[[symbol]] <- representative
        audit[[symbol]] <- group
    }
    list(data = do.call(rbind, selected), audit = do.call(rbind, 
        audit))
}
read_background <- function(file) {
    x <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
    columns <- intersect(c("gene", "gene_symbol", "gene_name"), 
        names(x))
    if (!length(columns) && ncol(x) == 1) 
        columns <- names(x)
    if (!length(columns)) {
        stop("Background needs a gene, gene_symbol, or gene_name column: ", 
            file)
    }
    genes <- unique(toupper(trimws(as.character(x[[columns[1]]]))))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (!length(genes)) 
        stop("Empty background: ", file)
    genes
}
prepare_genes <- function(file, cancer = FALSE) {
    x <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
    required <- c("log2FoldChange", "padj", if (cancer) {
        c("gene_id", "gene_name", "gene_type")
    } else "gene")
    missing <- setdiff(required, names(x))
    if (length(missing)) {
        stop("Missing columns in ", file, ": ", paste(missing, 
            collapse = ", "))
    }
    if (cancer) {
        x <- x[!is.na(x$gene_type) & x$gene_type == "protein_coding", 
            , drop = FALSE]
        x$ensembl_id <- sub("\\.[0-9]+$", "", x$gene_id)
    }
    else {
        x$source_row <- seq_len(nrow(x))
    }
    x$gene_symbol <- toupper(trimws(as.character(x[[if (cancer) "gene_name" else "gene"]])))
    x$gene_symbol[x$gene_symbol == ""] <- NA_character_
    x$padj <- validate_numeric(x$padj, "padj", file)
    x$log2FoldChange <- validate_numeric(x$log2FoldChange, "log2FoldChange", 
        file)
    x$regulation <- classify_regulation(x$padj, x$log2FoldChange)
    resolved <- resolve_symbol_duplicates(x, if (cancer) 
        "ensembl_id"
    else "source_row")
    list(data = resolved$data[!resolved$data$exclude_from_overlap, 
        , drop = FALSE], audit = resolved$audit)
}
prepare_signature <- function(file, background_file) {
    signature <- prepare_genes(file)
    background <- read_background(background_file)
    outside <- setdiff(signature$audit$gene_symbol, background)
    if (length(outside)) {
        stop("Signature genes absent from background: ", paste(head(outside, 
            10), collapse = ", "), " in ", background_file)
    }
    conflicts <- signature$audit$gene_symbol[signature$audit$exclude_from_overlap]
    signature$background <- setdiff(background, conflicts)
    signature
}
prepare_cancer <- function(file) {
    prepare_genes(file, cancer = TRUE)
}

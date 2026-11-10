# GSE45276 (Kennedy 2011): IMR90 transduced with activated AKT, RAS or control, n = 4,
# Affymetrix HG-U133 Plus 2 (RMA). Array data, so limma rather than the house DESeq2
# function; same DEG thresholds as the house (BH < 0.05, |log2FC| > log2(1.5)).
# Probes collapsed to genes by the highest-mean probe. Run from the pancancer-deg root.
source("exploratory/paths.R")   # repo-relative inputs; run from the pancancer-deg repo root
add_lib(BULK_LIB)
suppressMessages(library(limma))
G <- "data/exploratory/geo_oncogenes"
out_dir <- "results/exploratory/ras_vs_ois/contrast_degs"
l <- readLines(gzfile(file.path(G, "GSE45276_series_matrix.txt.gz")))
titles <- gsub('"', "", strsplit(grep("^!Sample_title", l, value = TRUE), "\t")[[1]][-1])
b <- grep("^!series_matrix_table_begin", l); e <- grep("^!series_matrix_table_end", l)
x <- read.delim(text = l[(b + 1):(e - 1)], check.names = FALSE, row.names = 1)
colnames(x) <- titles
ann <- read.delim(gzfile(file.path(G, "GPL570.annot.gz")), skip = grep("^ID\t", readLines(gzfile(file.path(G, "GPL570.annot.gz")), n = 200)) - 1,
                  check.names = FALSE, quote = "")
sym <- ann$`Gene symbol`[match(rownames(x), ann$ID)]
keep <- !is.na(sym) & sym != "" & !grepl("///", sym)
x <- x[keep, ]; sym <- sym[keep]
o <- order(rowMeans(x), decreasing = TRUE)
x <- x[o, ]; sym <- sym[o]
x <- x[!duplicated(sym), ]; rownames(x) <- sym[!duplicated(sym)]
grp <- factor(sub(" [0-9]+$", "", colnames(x)), levels = c("CONTROL", "AKT", "RAS"))
des <- model.matrix(~ 0 + grp); colnames(des) <- levels(grp)
fit <- lmFit(as.matrix(x), des)
for (arm in c("AKT", "RAS")) {
  f2 <- eBayes(contrasts.fit(fit, makeContrasts(contrasts = paste0(arm, "-CONTROL"), levels = des)))
  tt <- topTable(f2, number = Inf, sort.by = "none")
  d <- data.frame(gene = rownames(tt), log2FoldChange = tt$logFC, padj = tt$adj.P.Val,
                  sig = ifelse(tt$adj.P.Val < 0.05 & abs(tt$logFC) > log2(1.5), "y", "n"),
                  contrast = paste0("IMR90_array_", arm))
  write.csv(d, file.path(out_dir, paste0("IMR90_array_", arm, ".csv")), row.names = FALSE)
  cat(sprintf("IMR90_array_%s: genes=%d up=%d down=%d\n", arm, nrow(d),
              sum(d$sig == "y" & d$log2FoldChange > 0), sum(d$sig == "y" & d$log2FoldChange < 0)))
}

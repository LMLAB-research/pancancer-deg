# Non-RAS oncogene arms with the house find_degs_new():
#  - GSE130099 / GSE130100 (Chan 2020): BJ-hTERT myrAKT1 (day 6) vs control; HRAS-G12V (day 14)
#    vs pBabe. Separate series, so AKT and RAS are compared only through their own controls.
#  - GSE185512 (Zimmerli 2022): dox-inducible MYC / cyclin E / empty vector in RPE-1 (p53 wt and
#    p53-/-), dox vs no dox at 120 h (n = 2). QuantSeq 3' counts. Senescence not demonstrated.
# Raw counts from GEO supplementary files. Run from the pancancer-deg repo root.
source("exploratory/ras_vs_ois/house_env.R")
G <- "data/exploratory/geo_oncogenes"
out_dir <- "results/exploratory/ras_vs_ois/contrast_degs"

run <- function(m, g1, g2, name) {
  m <- m[, c(g1, g2)]
  cd <- S4Vectors::DataFrame(external_id = colnames(m),
                             grp = factor(c(rep("G1", length(g1)), rep("G2", length(g2))), levels = c("G2", "G1")),
                             row.names = colnames(m))
  s <- SummarizedExperiment(assays = list(counts = m), colData = cd)
  s <- s[rownames(filter_counts(assay(s, "counts"))), ]
  storage.mode(assay(s, "counts")) <- "integer"
  d <- find_degs_new(rse_obj = s, group_col = "grp", group_1 = "G1", group_2 = "G2")$degs
  d$contrast <- name
  write.csv(d, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  cat(sprintf("%-30s n=%d vs %d up=%d down=%d\n", name, length(g1), length(g2),
              sum(d$sig == "y" & d$log2FoldChange > 0), sum(d$sig == "y" & d$log2FoldChange < 0)))
}
sym_counts <- function(df, sym_col, count_cols) {
  df <- df[!is.na(df[[sym_col]]) & df[[sym_col]] != "", ]
  m <- as.matrix(df[, count_cols]); storage.mode(m) <- "numeric"
  rowsum(round(m), df[[sym_col]])
}

# --- AKT and RAS in BJ-hTERT ---
# The GEO tables have unquoted commas in the description column; counts are fields 3-8,
# so parse only the first 8 fields of each line.
read_first8 <- function(f) {
  l <- readLines(gzfile(f))
  p <- lapply(strsplit(l, ",", fixed = TRUE), `[`, 1:8)
  df <- as.data.frame(do.call(rbind, p[-1]), stringsAsFactors = FALSE)
  colnames(df) <- p[[1]]
  df
}
akt <- read_first8(file.path(G, "GSE130099_AKT_control_all_columnsDEseq2.csv.gz"))
m <- sym_counts(akt, "hgnc_symbol", c("BJ.T1", "BJ.T2", "BJ.T3", "BJ.T.AKT1", "BJ.T.AKT2", "BJ.T.AKT3"))
run(m, c("BJ.T.AKT1", "BJ.T.AKT2", "BJ.T.AKT3"), c("BJ.T1", "BJ.T2", "BJ.T3"), "BJhT_AKT_d6")
ras <- read_first8(file.path(G, "GSE130100_HRAS_pBABE_all_columnsDEseq2.csv.gz"))
m <- sym_counts(ras, "hgnc_symbol", c("pBABE.NTC.1", "pBABE.NTC.2", "pBABE.NTC.3", "HRAS.NTC.1", "HRAS.NTC.2", "HRAS.NTC.3"))
run(m, c("HRAS.NTC.1", "HRAS.NTC.2", "HRAS.NTC.3"), c("pBABE.NTC.1", "pBABE.NTC.2", "pBABE.NTC.3"), "BJhT_HRAS_d14")

# --- MYC / cyclin E / empty in RPE-1 (Ensembl IDs -> symbols via the recount3 gene annotation) ---
x <- read.csv(file.path(G, "GSE185512_Expression_dataset.csv.gz"), check.names = FALSE)
ann <- as.data.frame(SummarizedExperiment::rowData(readRDS("data/exploratory/ras_vs_ois/SRP062871_gene_rse.rds")))
ann$id <- sub("\\..*", "", ann$gene_id)
ann <- ann[ann$gene_type == "protein_coding", ]
x$symbol <- ann$gene_name[match(x$probe, ann$id)]
smx <- readLines(gzfile(file.path(G, "GSE185512_series_matrix.txt.gz")))
titles <- strsplit(grep("^!Sample_title", smx, value = TRUE), "\t")[[1]][-1]
titles <- gsub('"', "", titles)
key <- sub(".*\\[(G_[0-9]+)\\]$", "\\1", titles)
cond <- sub("_[0-9]+ \\[G_[0-9]+\\]$", "", titles)
cond <- gsub("Doxycycline ", "Doxycycline", gsub("RPE-1-\\(P53-/-\\)", "RPE-1 (P53-/-)", cond))
names(cond) <- key
m <- sym_counts(x, "symbol", key)
pick <- function(line, construct, trt, t) key[cond == paste0(line, "-", construct, "_", trt, "_", t)]
for (line in c("RPE-1", "RPE-1 (P53-/-)")) {
  tag <- if (line == "RPE-1") "RPE1_p53wt" else "RPE1_p53ko"
  for (con in c("Myc", "CyclinE", "Empty")) {
    g1 <- pick(line, con, "Doxycycline", "120h"); g2 <- pick(line, con, "NO Dox", "120h")
    if (length(g1) >= 2 && length(g2) >= 2) run(m, g1, g2, paste0(tag, "_", con, "_dox120h"))
    else cat("missing samples for", line, con, "\n")
  }
  g1 <- pick(line, "Empty", "Hydroxyurea", "120h"); g2 <- pick(line, "Empty", "NO Dox", "120h")
  if (length(g1) >= 2 && length(g2) >= 2) run(m, g1, g2, paste0(tag, "_HU120h"))
}

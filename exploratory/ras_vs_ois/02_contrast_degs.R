# RAS vs OIS (exploratory, local only). Within-study DEG contrasts that separate
# oncogenic RAS/BRAF signalling from the senescent arrest, using the house
# find_degs_new() (DESeq2; padj < 0.05, |LFC| > log2(1.5)). Each contrast compares a
# RAS/BRAF (or other) condition with its matched proliferating control in the same
# genetic background. Run from the pancancer-deg repo root.
source("exploratory/ras_vs_ois/house_env.R")
in_dir <- "data/exploratory/ras_vs_ois"
out_dir <- "results/exploratory/ras_vs_ois/contrast_degs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# recount3 coverage -> AUC-scaled counts (40M), protein-coding genes by symbol.
prep_study <- function(srp) {
  rse <- readRDS(file.path(in_dir, paste0(srp, "_gene_rse.rds")))
  assays(rse)$counts <- round(recount3::transform_counts(rse, by = "auc", targetSize = 40e6))
  rd <- as.data.frame(rowData(rse))
  rse <- rse[rd$gene_type == "protein_coding" & !is.na(rd$gene_name), ]
  m <- rowsum(assay(rse, "counts"), rowData(rse)$gene_name)
  cd <- colData(rse)
  se <- SummarizedExperiment(assays = list(counts = m), colData = cd)
  se$external_id <- colnames(se)
  se$title <- se$sra.sample_title
  se
}

run_contrast <- function(se, name, g1_idx, g2_idx, batch = NULL) {
  s <- se[, c(which(g1_idx), which(g2_idx))]
  s$grp <- factor(c(rep("G1", sum(g1_idx)), rep("G2", sum(g2_idx))), levels = c("G2", "G1"))
  if (!is.null(batch)) s$batch <- factor(batch[c(which(g1_idx), which(g2_idx))])
  # house CPM filter: CPM >= 1 in >= 30% of the samples in this contrast
  s <- s[rownames(filter_counts(assay(s, "counts"))), ]
  assay(s, "counts") <- as.matrix(assay(s, "counts"))
  storage.mode(assay(s, "counts")) <- "integer"
  res <- find_degs_new(rse_obj = s, group_col = "grp", group_1 = "G1", group_2 = "G2",
                       batch_col = if (is.null(batch)) NULL else "batch")
  d <- res$degs
  d$contrast <- name
  d$n_g1 <- sum(g1_idx); d$n_g2 <- sum(g2_idx)
  write.csv(d, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  cat(sprintf("%-22s n=%d vs %d  genes=%d  up=%d down=%d\n", name, sum(g1_idx), sum(g2_idx),
              nrow(d), sum(d$sig == "y" & d$log2FoldChange > 0), sum(d$sig == "y" & d$log2FoldChange < 0)))
  invisible(d)
}

# --- SRP065206: IMR90 HRAS-G12V senescence under knockdowns / BRD4 inhibition ---
se <- prep_study("SRP065206")
base <- sub("_[0-9]+$", "", se$title)  # e.g. "Senescent IMR90 shp53/Rb"
for (kd in c("shRen", "shp53", "shp16/p21", "shp53/Rb", "shp65", "shBrd4", "DMSO", "JQ1")) {
  tag <- gsub("[^A-Za-z0-9]", "", kd)
  run_contrast(se, paste0("IMR90_RAS_", tag),
               base == paste("Senescent IMR90", kd), base == paste("Proliferating IMR90", kd))
}
run_contrast(se, "IMR90_quiescent_shRen", base == "Quiescent IMR90 shRen", base == "Proliferating IMR90 shRen")

# --- SRP066947: TIG3 BRAF-V600E OIS vs OIS bypass, day 4 and day 10; serum starvation ---
se <- prep_study("SRP066947"); t <- se$title
for (d in c("4", "10")) {
  cyc <- startsWith(t, paste0("Cycling Day ", d, " "))
  run_contrast(se, paste0("TIG3_BRAF_OIS_d", d), startsWith(t, paste0("BRAF OIS Day ", d, " ")), cyc)
  run_contrast(se, paste0("TIG3_BRAF_bypass_d", d), startsWith(t, paste0("BRAF Senescence Bypass Day ", d, " ")), cyc)
}
run_contrast(se, "TIG3_quiescent_serum", startsWith(t, "Quiescent Serum"), startsWith(t, "Cycling Control Serum"))

# --- SRP172671: BJ-hTERT RAS senescence vs clones that escaped (lines A and B) ---
se <- prep_study("SRP172671"); t <- se$title
line <- substr(t, 1, 1)
parent <- t %in% c("A_BJ", "B_BJ")
run_contrast(se, "BJ_RAS_OIS", t %in% c("A_BJ_ras", "B_BJ_ras"), parent, batch = line)
run_contrast(se, "BJ_RAS_escaped_clones", grepl("_clone_", t), parent, batch = line)

# --- SRP062872: RAS-induced vs etoposide-induced senescence, same ER:RAS IMR90 cells ---
se <- prep_study("SRP062872"); t <- se$title
run_contrast(se, "IMR90_RIS", startsWith(t, "RIS"), startsWith(t, "Control"))
run_contrast(se, "IMR90_DDIS", startsWith(t, "DDIS"), startsWith(t, "Control"))

# --- SRP127037: ER:RAS IMR90-hTERT time course (day 4 onset, day 10), +/- IL-1R knockdown ---
se <- prep_study("SRP127037"); t <- se$title
for (sh in c("shScr", "shIL1R")) {
  d0 <- startsWith(t, sh) & endsWith(t, "d0")
  run_contrast(se, paste0("IMR90hT_RAS_d4_", sh), startsWith(t, sh) & endsWith(t, "d4"), d0)
  run_contrast(se, paste0("IMR90hT_RAS_d10_", sh), startsWith(t, sh) & endsWith(t, "d10"), d0)
}

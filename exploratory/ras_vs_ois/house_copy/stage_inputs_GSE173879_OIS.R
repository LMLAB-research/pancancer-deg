# LOCAL COPY for pancancer-deg exploratory RAS vs OIS (2026-11-02) of senescence_stress
# house_rerun/scripts/stage_inputs.R. Only two changes: (a) output root O is in pancancer-deg
# (data/exploratory/house_GSE173879_OIS/) so nothing is written into senescence_stress; (b) the GSE173879
# plate pick stages the OIS arms (OISON_1 = GSM5282311 'OIS 4OHT', in 4-OHT at harvest; OISOFF_1 = GSM5282312
# 'OIS', 4-OHT-induced then DMSO at harvest) with the same two control plates as the house RS run.
#!/usr/bin/env Rscript
# stage_inputs.R -- INPUT STAGING ONLY for the house preprocessing (scRNA_seq_preprocess_final.R, modules_intact copy).
# WHY. Roberto (2026-10-15): re-run the external senescence datasets through the house scripts; the only allowed changes
#   are (1) dataset inputs staged in the GEO per-plate layout the house script reads, (2) each paper's published QC,
#   (3) the nuclear-integrity keep lists. This script does (1), the published-QC steps the house script cannot express
#   (OIS: HTODemux singlets + DoubletFinder; GSE250041: ADT total <= 20,000), and writes (3) with the matrix rule
#   (MALAT1+NEAT1 share >= own-control p05 and not [> control p99.5 and ribosomal share < control p01]), thresholds from
#   the QC-passing control cells. Loaders, QC values, HTODemux/DoubletFinder code and the integrity rule are copied from
#   ../../external_preprocess.R (CFG); nothing else from that script is used.
# USAGE: Rscript stage_inputs.R <GSE>
# OUT:   house_rerun/<GSE>/staged/<GSM>_<PLATE>_{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}
#        house_rerun/<GSE>/integrity/integrity_keep_<PLATE>.txt, integrity_thresholds.csv
#        house_rerun/<GSE>/staging_counts.csv (barcodes per plate before the house script), staging_map.csv
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(DoubletFinder) })
options(future.globals.maxSize = 16 * 1024^3); set.seed(1)
GSE <- commandArgs(TRUE)[1]
S <- "/Users/ro/Academic/10_Research/senescence_stress/modules/scrna_seq/cs_dataset_survey/"; R <- paste0(S, "data/raw/"); I <- paste0(S, "data/interim/")
O <- "/Users/ro/code/academic/pancancer-deg/data/exploratory/house_GSE173879_OIS/"; ST <- paste0(O, "staged/"); IN <- paste0(O, "integrity/")
dir.create(ST, showWarnings = FALSE, recursive = TRUE); dir.create(IN, showWarnings = FALSE, recursive = TRUE)
# ---- from external_preprocess.R (verbatim) ----
rd <- function(pre) { m <- readMM(paste0(pre, "matrix.mtx.gz")); f <- read.delim(gzfile(paste0(pre, "features.tsv.gz")), header = FALSE)
  colnames(m) <- readLines(gzfile(paste0(pre, "barcodes.tsv.gz"))); list(m = as(m, "CsparseMatrix"), f = f) }
gex <- function(x) { ge <- if (ncol(x$f) >= 3) x$f$V3 == "Gene Expression" else rep(TRUE, nrow(x$f)); m <- x$m[ge, ]; g <- x$f$V2[ge]
  k <- !(g %in% g[duplicated(g)]); m <- m[k, ]; rownames(m) <- g[k]; m }
WECHTER <- list(g = c(300, Inf), u = c(1000, 120000), mt = 12)
CFG <- list(
  GSE245807 = list(qc = list(g = c(-Inf, Inf), u = c(-Inf, Inf), mt = 5, mt_le = TRUE)),
  GSE250041 = list(qc = list(g = c(300, 12000), u = c(2000, 150000), mt = 10, adt = 20000)),
  GSE282425 = list(qc = list(g = c(400, Inf), u = c(1000, Inf), mt = 12, mt_le = TRUE)),
  GSE306635 = list(qc = WECHTER), GSE173879 = list(qc = WECHTER))
q <- CFG[[GSE]]$qc
# ---- helpers (staging only) ----
write_plate <- function(gsm, plate, m, f) {   # m: full feature x barcode matrix (all feature types), f: feature table
  pre <- paste0(ST, gsm, "_", plate, "_")
  # 10x v3 feature files: keep Gene Expression rows, write 2 columns (id, symbol). The house script's GE filter reads
  # features with read.table (whitespace), which splits "Gene Expression" and would drop every feature (ADAPTATIONS.md).
  if (ncol(f) >= 3) { ge <- f$V3 == "Gene Expression"; m <- m[ge, ]; f <- f[ge, 1:2] }
  tmp <- tempfile(); writeMM(m, tmp); system2("gzip", c("-c", tmp), stdout = paste0(pre, "matrix.mtx.gz")); unlink(tmp)
  g <- gzfile(paste0(pre, "features.tsv.gz"), "w"); write.table(f, g, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE); close(g)
  g <- gzfile(paste0(pre, "barcodes.tsv.gz"), "w"); writeLines(colnames(m), g); close(g)
}
stage_counts <- list(); map <- list()
note <- function(plate, step, n) stage_counts[[length(stage_counts) + 1]] <<- data.frame(plate = plate, step = step, n = n)
plates <- list()   # plate -> list(gsm, m (all features), f, source, cell_name_suffix)

if (GSE == "GSE245807") {   # one library, hashtags; each hashtag sample = one plate (ADAPTATIONS.md)
  x <- rd(paste0(R, "GSE245807/GSE245807_")); hto <- x$m[x$f$V3 == "Antibody Capture", ]; rownames(hto) <- x$f$V1[x$f$V3 == "Antibody Capture"]
  o <- CreateSeuratObject(gex(x)); o[["HTO"]] <- CreateAssayObject(counts = hto)
  o <- NormalizeData(o, assay = "HTO", normalization.method = "CLR", verbose = FALSE); o <- HTODemux(o, assay = "HTO", positive.quantile = 0.99, verbose = FALSE)
  n_lib <- ncol(o); o <- subset(o, subset = HTO_classification.global == "Singlet")
  o$sample <- as.character(o$hash.ID); o$condition <- sub("-[0-9]$|_[0-9]$", "", o$sample); o$replicate <- sub(".*[-_]", "", o$sample)
  o <- subset(o, subset = condition %in% c("Grow", "RIS"))
  # published QC then DoubletFinder on QC-passing cells (as in external_preprocess.R; Olan 2024)
  o$percent.mt <- PercentageFeatureSet(o, pattern = "^MT-"); okq <- o$percent.mt <= q$mt
  d <- subset(o, cells = colnames(o)[okq])
  d <- NormalizeData(d, verbose = FALSE); d <- FindVariableFeatures(d, verbose = FALSE); d <- ScaleData(d, verbose = FALSE); d <- RunPCA(d, verbose = FALSE)
  sw <- summarizeSweep(paramSweep(d, PCs = 1:20, sct = FALSE), GT = FALSE); pk <- find.pK(sw); pK <- as.numeric(as.character(pk$pK[which.max(pk$BCmetric)]))
  nexp <- round(0.008 * ncol(d) / 1000 * ncol(d)); d <- doubletFinder(d, PCs = 1:20, pN = 0.25, pK = pK, nExp = nexp, sct = FALSE)
  cl <- d@meta.data[[grep("^DF.classifications", colnames(d@meta.data), value = TRUE)[1]]]
  df_doublets <- colnames(d)[cl == "Doublet"]
  cat(sprintf("DoubletFinder: pK %.3f, nExp %d, doublets %d of %d QC-passing singlets\n", pK, nexp, length(df_doublets), ncol(d)))
  pm <- c("Grow-1" = "CTRL_1", "Grow-2" = "CTRL_2", "RIS-1" = "RIS_1", "RIS-2" = "RIS_2")
  for (s in names(pm)) {
    cells_s <- colnames(o)[o$sample == s]; keep_s <- setdiff(cells_s, df_doublets)
    note(pm[[s]], "library barcodes (whole library)", n_lib); note(pm[[s]], paste0("HTODemux singlets ", s), length(cells_s))
    note(pm[[s]], "after DoubletFinder (doublets among published-QC-passing cells removed)", length(keep_s))
    plates[[pm[[s]]]] <- list(gsm = "GSM0", m = x$m[, keep_s], f = x$f, source = paste0("GSE245807 single library, hashtag ", s))
  }
  thr_note <- sprintf("DoubletFinder pK %.3f nExp %d", pK, nexp)
} else if (GSE == "GSE250041") {   # one library per condition (CITE-seq); ADT total <= 20,000 (Abdelmohsen 2024)
  for (s in c("Proliferating", "Senescent")) {
    x <- rd(paste0(R, "GSE250041/GSE250041_", s, "_")); adt <- Matrix::colSums(x$m[x$f$V3 != "Gene Expression", , drop = FALSE])
    pl <- if (s == "Proliferating") "CTRL_1" else "IR_1"; keep_s <- colnames(x$m)[adt <= q$adt]
    note(pl, "library barcodes", ncol(x$m)); note(pl, "ADT total <= 20,000 (published QC, applied at staging)", length(keep_s))
    plates[[pl]] <- list(gsm = "GSM0", m = x$m[, keep_s], f = x$f, source = paste0("GSE250041_", s))
  }
} else if (GSE == "GSE282425") {   # one library, sample from barcode suffix (external_preprocess.R mapping)
  x <- rd(paste0(R, "GSE282425/GSE282425_")); suf <- sub(".*[-_]", "", colnames(x$m))
  lab <- c("3" = "EP2D", "4" = "EP3D", "1" = "DS2D", "2" = "DS3D")
  for (s in c("EP2D", "DS2D")) { pl <- if (s == "EP2D") "CTRL_1" else "RS_1"; cs <- colnames(x$m)[unname(lab[suf]) == s]
    note(pl, paste0("library barcodes with sample ", s), length(cs))
    plates[[pl]] <- list(gsm = "GSM0", m = x$m[, cs], f = x$f, source = paste0("GSE282425 single library, barcode suffix sample ", s)) }
} else if (GSE == "GSE306635") {
  pick <- c(Proliferating1 = "CTRL_1", Proliferating2 = "CTRL_2", "1Week1" = "ETO_1", "1Week2" = "ETO_2")
  for (s in names(pick)) { pre <- Sys.glob(paste0(I, "GSE306635/GSM*_", s, "_matrix.mtx.gz")); x <- rd(sub("matrix.mtx.gz$", "", pre))
    note(pick[[s]], "library barcodes", ncol(x$m))
    plates[[pick[[s]]]] <- list(gsm = sub("_.*", "", basename(pre)), m = x$m, f = x$f, source = basename(pre)) }
} else if (GSE == "GSE173879") {
  pick <- c(Proliferating_Clone_1_BioRep1 = "CTRL_1:GSM5282307", Proliferating_Clone_2_BioRep1 = "CTRL_2:GSM5282308",
            OIS_4OHT_Clone_1 = "OISON_1:GSM5282311", OIS_Clone_1 = "OISOFF_1:GSM5282312")
  for (s in names(pick)) { d <- paste0(I, "GSE173879/", s, "/"); m <- readMM(paste0(d, "matrix.mtx")); f <- read.delim(paste0(d, "genes.tsv"), header = FALSE)
    colnames(m) <- readLines(paste0(d, "barcodes.tsv")); p <- strsplit(pick[[s]], ":")[[1]]
    note(p[1], "library barcodes", ncol(m))
    plates[[p[1]]] <- list(gsm = p[2], m = as(m, "CsparseMatrix"), f = f, source = s) }
}

# ---- write the GEO per-plate layout ----
for (pl in names(plates)) { P <- plates[[pl]]; write_plate(P$gsm, pl, P$m, P$f)
  map[[length(map) + 1]] <- data.frame(plate = pl, staged_prefix = paste0(P$gsm, "_", pl), source = P$source, n_barcodes = ncol(P$m)) }

# ---- integrity keep lists: matrix rule, thresholds from the QC-passing control (house QC = published QC) ----
ctrl_pl <- grep("^CTRL", names(plates), value = TRUE)
feat <- lapply(names(plates), function(pl) { m <- gex(list(m = plates[[pl]]$m, f = plates[[pl]]$f)); tot <- Matrix::colSums(m)
  nf <- Matrix::colSums(m > 0); mt <- 100 * Matrix::colSums(m[grep("^MT-", rownames(m)), , drop = FALSE]) / tot
  data.frame(plate = pl, barcode = colnames(m), nCount = tot, nFeature = nf, pct_mt = mt,
             nuc = Matrix::colSums(m[intersect(c("MALAT1", "NEAT1"), rownames(m)), , drop = FALSE]) / tot,
             ribo = Matrix::colSums(m[grep("^RP[SL][0-9]", rownames(m)), , drop = FALSE]) / tot) })
feat <- do.call(rbind, feat)
mt_ok <- if (isTRUE(q$mt_le)) feat$pct_mt <= q$mt else feat$pct_mt < q$mt
feat$qc_pass <- feat$nFeature >= q$g[1] & feat$nFeature <= q$g[2] & feat$nCount >= q$u[1] & feat$nCount <= q$u[2] & mt_ok
cc <- feat[feat$plate %in% ctrl_pl & feat$qc_pass, ]
thr <- c(nuc_p05 = unname(quantile(cc$nuc, .05)), nuc_p995 = unname(quantile(cc$nuc, .995)), ribo_p01 = unname(quantile(cc$ribo, .01)))
feat$integrity <- ifelse(feat$nuc < thr["nuc_p05"], "nucleus-poor", ifelse(feat$nuc > thr["nuc_p995"] & feat$ribo < thr["ribo_p01"], "bare nucleus", "intact"))
for (pl in names(plates)) {
  k <- feat$barcode[feat$plate == pl & feat$integrity == "intact"]
  # the house script suffixes control barcodes with _<plate> when there is more than one control library
  if (pl %in% ctrl_pl && length(ctrl_pl) > 1) k <- paste0(k, "_", pl)
  writeLines(k, paste0(IN, "integrity_keep_", pl, ".txt"))
  x <- feat[feat$plate == pl & feat$qc_pass, ]
  note(pl, "[info] published-QC-passing barcodes (staging check)", nrow(x))
  note(pl, "[info] of which intact by matrix rule", sum(x$integrity == "intact"))
}
write.csv(data.frame(t(thr), n_ctrl_qc = nrow(cc), ctrl_plates = paste(ctrl_pl, collapse = ";")), paste0(IN, "integrity_thresholds.csv"), row.names = FALSE)
write.csv(do.call(rbind, stage_counts), paste0(O, "staging_counts.csv"), row.names = FALSE)
write.csv(do.call(rbind, map), paste0(O, "staging_map.csv"), row.names = FALSE)
print(thr); print(do.call(rbind, stage_counts)); cat("staged", GSE, "\n")

# config.R
#
# Single source of truth for where this project reads/writes data on THIS
# machine. Nothing else in bulk_RNA_seq/ should hardcode a path.
#
# Override by setting the env var before running R, e.g.:
#   BULK_RNASEQ_DIR=/path/to/systems_analysis_arrest Rscript meta_analysis/01_build_cs_cq_object.R
# Without it, this assumes you launched R/Rscript with the working directory
# set to bulk_RNA_seq/ itself (i.e. `cd bulk_RNA_seq && Rscript meta_analysis/01_...R`),
# and resolves the project root as one level up from the current working dir.
# Setting BULK_RNASEQ_DIR is the reliable option if you run scripts any other way.

PROJECT_DIR <- Sys.getenv("BULK_RNASEQ_DIR", unset = dirname(getwd()))

DATA_DIR         <- file.path(PROJECT_DIR, "Final")
SAVE_DIR_CSV     <- file.path(DATA_DIR, "SI_tables")
SAVE_DIR_FIGURE  <- file.path(DATA_DIR, "Figures")
SAVE_DIR_FIGURE_SI <- file.path(DATA_DIR, "SI_figures")
ERP021140_DIR    <- file.path(DATA_DIR, "ERP021140")
RERUN_DIR        <- file.path(PROJECT_DIR, "bulk_RNA_seq", "rerun_outputs")

# ---- inputs that are NOT produced by the pipeline and NOT in the repo ------------
# EXTERNAL_DIR holds them; see REPRODUCE.md for where each comes from and its md5.
#   models/  the three tAge clocks (.pkl, Tyshkovskiy et al. - third-party, never commit)
#   geo/     GSE175533_TPM.xlsx (GEO supplementary file, input to meta_analysis/20)
EXTERNAL_DIR     <- Sys.getenv("BULK_EXTERNAL_DIR",
                               unset = file.path(PROJECT_DIR, "bulk_RNA_seq", "external_inputs"))
MODELS_DIR       <- Sys.getenv("TAGE_MODELS_DIR", unset = file.path(EXTERNAL_DIR, "models"))

# The sample sheet meta_analysis/01 builds the object from: the SI table AS SUBMITTED
# (md5 7999c22a..., = archive_original/Final_original/SI_tables/study_info_all.csv).
# meta_analysis/15 later corrected Final/SI_tables/study_info_all.csv in place (tissue,
# cell line, immortalised) and kept the submitted version under this name. Building from
# the corrected table gives the same counts, tAge and DEGs, but meta_analysis/10 then
# records the corrected labels as the "curated" ones, so its audit columns change.
STUDY_INFO_CSV   <- Sys.getenv("BULK_STUDY_INFO",
                               unset = file.path(SAVE_DIR_CSV, "study_info_all.pre_annotation_fix.csv"))

# The Python that reticulate uses for tAge prediction (needs joblib, scikit-learn,
# pandas - see requirements.txt). Must be the SAME interpreter the Python scripts run
# under, so R-side and Python-side predictions come from one environment. Earlier
# versions pointed at the single-cell project's .venv, which had different versions.
PYTHON_BIN       <- Sys.getenv("TAGE_PYTHON",
                               unset = file.path(PROJECT_DIR, "bulk_RNA_seq", ".venv", "bin", "python"))

for (d in c(SAVE_DIR_CSV, SAVE_DIR_FIGURE, SAVE_DIR_FIGURE_SI, ERP021140_DIR, RERUN_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Ensembl biomart host -- NOT used for the protein-coding gene dictionary
# anymore (see get_ensembl_release_pc() in functions.R). Kept for any script
# still calling useMart() directly for other attributes; NULL = biomaRt's
# current live default host.
#
# WHY NOT THE ORIGINAL ARCHIVE: the original analysis pinned biomaRt to
# apr2020.archive.ensembl.org (Ensembl release ~100) so the protein-coding
# gene list stayed fixed across reruns. That interactive archive+biomart
# service is now decommissioned (redirects to a generic "archive retired"
# page; every pre-2021 mirror we tried -- jan2020/apr2019/nov2020/may2021 --
# is dead the same way). BUT Ensembl's plain FTP file archive for the same
# release is still live (ftp.ensembl.org/pub/release-100/), so
# get_ensembl_release_pc() pins the exact Ensembl-100 protein-coding gene set
# by parsing that GTF directly instead -- verified byte-for-byte: 13,681/13,681
# (100%) of the gene symbols in the original cq_samples.rds are reproduced
# exactly this way. No gene-universe gap remains for this part of the pipeline.
ENSEMBL_ARCHIVE_HOST <- NULL  # NULL = biomaRt's current live default host
ENSEMBL_PINNED_RELEASE <- 100L  # the release the original analysis used

message(sprintf("[config] PROJECT_DIR = %s", PROJECT_DIR))
message(sprintf("[config] RERUN_DIR   = %s", RERUN_DIR))
message(sprintf("[config] MODELS_DIR  = %s", MODELS_DIR))

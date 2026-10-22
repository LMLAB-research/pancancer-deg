# Input paths. An optional base directory supports the original two-folder layout.
repo_dir <- dirname(project_dir)
base_dir <- Sys.getenv("TCGA_SENESCENCE_BASE_DIR", unset = "")
senescence_input_dir <- if (nzchar(base_dir)) {
  file.path(base_dir, "meta_analysis")
} else {
  file.path(repo_dir, "data", "senescence_meta_analysis_deg")
}
cancer_input_dir <- if (nzchar(base_dir)) {
  file.path(base_dir, "TCGA_deg")
} else {
  file.path(repo_dir, "results", "tables", "TCGA_deg")
}
results_dir <- Sys.getenv("TCGA_SENESCENCE_RESULTS_DIR", unset = file.path(project_dir, "results"))

expected_signatures <- c("OIS", "RS", "SIPS", "CICQ", "SSCQ")
cancer_file_pattern <- "^TCGA-[A-Za-z0-9-]+_DGE_Raw_Results_Tumor_vs_Normal_Adjusted\\.csv$"

# Thresholds 
padj_threshold <- 0.05
lfc_threshold <- log2(1.5)

# Heatmap settings

cancer_order <- c(
  "HNSC",
  "ESCA",
  "LUAD",
  "LUSC",
  "STAD",
  "COAD",
  "READ",
  "LIHC",
  "CHOL",
  "KICH",
  "KIRC",
  "KIRP",
  "BLCA",
  "THCA",
  "BRCA",
  "PRAD",
  "UCEC"
)

cancer_labels <- c(
  HNSC = "HNSC (head and neck)",
  ESCA = "ESCA (esophagus)",
  LUAD = "LUAD (lung adenocarcinoma)",
  LUSC = "LUSC (lung squamous)",
  STAD = "STAD (stomach)",
  COAD = "COAD (colon)",
  READ = "READ (rectum)",
  LIHC = "LIHC (liver)",
  CHOL = "CHOL (bile duct)",
  KICH = "KICH (kidney chromophobe)",
  KIRC = "KIRC (kidney clear cell)",
  KIRP = "KIRP (kidney papillary)",
  BLCA = "BLCA (bladder)",
  THCA = "THCA (thyroid)",
  BRCA = "BRCA (breast)",
  PRAD = "PRAD (prostate)",
  UCEC = "UCEC (uterus)"
)

class_order <- c( "UP-UP", "UP-DOWN", "DOWN-UP", "DOWN-DOWN" )

# Include the current signature name in the heatmap column labels.
make_class_labels <- function(signature) {
  c(
    "UP-UP" = paste0( "Up in cancer,\nup in ", signature ),
    "UP-DOWN" = paste0( "Up in cancer,\ndown in ", signature ),
    "DOWN-UP" = paste0( "Down in cancer,\nup in ", signature ),
    "DOWN-DOWN" = paste0( "Down in cancer,\ndown in ", signature )
  )
}

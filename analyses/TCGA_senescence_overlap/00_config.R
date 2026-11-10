#### Paths ####
# Scripts are run from the repository root.

analysis_name <- "TCGA_senescence_overlap"

# Senescence/quiescence signatures and their backgrounds.
senescence_input_dir <- file.path("data", "senescence_meta_analysis_deg")
# Full TCGA DESeq2 tables from analyses/TCGA_deg (step 04).
cancer_input_dir <- file.path("results", "tables", "TCGA_deg")

results_dir <- file.path("results", "tables", analysis_name)
figures_dir <- file.path("results", "figures", analysis_name)

expected_signatures <- c("OIS", "RS", "SIPS", "CICQ", "SSCQ")
cancer_file_pattern <- "^TCGA-[A-Za-z0-9-]+_DGE_Raw_Results_Tumor_vs_Normal_Adjusted\\.csv$"

# Thresholds 
padj_threshold <- 0.05
lfc_threshold <- log2(1.5)

# Fisher test sidedness: "two.sided" tests for enrichment and depletion;
# "greater" tests enrichment only (the original GeneOverlap behaviour).
fisher_alternative <- "two.sided"

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

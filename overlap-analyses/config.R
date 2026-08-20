# Input data
cancer_input_dir <- Sys.getenv(
  "TCGA_OIS_INPUT_DIR",
  unset = file.path(script_dir, "data")
)
ois_file <- file.path(cancer_input_dir, "senescence_OIS_genes.csv")

# Output data
results_dir <- Sys.getenv(
  "TCGA_OIS_RESULTS_DIR",
  unset = file.path(script_dir, "results")
)

# Analysis settings
cancer_file_pattern <- paste0(
  "^TCGA-[A-Za-z0-9-]+_DGE_Raw_Results_",
  "Tumor_vs_Normal_Adjusted\\.csv$"
)
padj_threshold <- 0.05
lfc_threshold <- log2(1.5)

# Heatmap settings
cancer_order <- c(
  "HNSC", "ESCA", "LUAD", "LUSC", "STAD", "COAD", "READ", "LIHC",
  "CHOL", "KICH", "KIRC", "KIRP", "BLCA", "THCA", "BRCA", "PRAD", "UCEC"
)

cancer_labels <- c(
  HNSC = "HNSC (head and neck)", ESCA = "ESCA (esophagus)",
  LUAD = "LUAD (lung adenocarcinoma)", LUSC = "LUSC (lung squamous)",
  STAD = "STAD (stomach)", COAD = "COAD (colon)", READ = "READ (rectum)",
  LIHC = "LIHC (liver)", CHOL = "CHOL (bile duct)",
  KICH = "KICH (kidney chromophobe)", KIRC = "KIRC (kidney clear cell)",
  KIRP = "KIRP (kidney papillary)", BLCA = "BLCA (bladder)",
  THCA = "THCA (thyroid)", BRCA = "BRCA (breast)",
  PRAD = "PRAD (prostate)", UCEC = "UCEC (uterus)"
)

class_order <- c("UP-UP", "UP-DOWN", "DOWN-UP", "DOWN-DOWN")

class_labels <- c(
  "UP-UP" = "Up in cancer,\nup in OIS",
  "UP-DOWN" = "Up in cancer,\ndown in OIS",
  "DOWN-UP" = "Down in cancer,\nup in OIS",
  "DOWN-DOWN" = "Down in cancer,\ndown in OIS"
)

library(TCGAbiolinks)
library(dplyr)

# 1. Fetch all TCGA project IDs
tcga_projects <- TCGAbiolinks:::getGDCprojects()$project_id %>% 
  grep("^TCGA", ., value = TRUE) %>% 
  sort()

# Initialize data frame to hold sample type validation statuses
project_status_df <- data.frame(
  cancer       = character(),
  has_tumor    = logical(),
  has_normal   = logical(),
  stringsAsFactors = FALSE
)

message("Scanning cohorts for Tumor and Normal sample types...")

# 2. Check the sample types present in each project
for (proj in tcga_projects) {
  tryCatch({
    query <- GDCquery(
      project       = proj,
      data.category = "Transcriptome Profiling",
      data.type     = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
    
    meta_df <- getResults(query)
    
    # Check if any rows match our simplified groups
    contains_tumor  <- any(meta_df$sample_type %in% c("Primary Tumor", "Metastatic"))
    contains_normal <- any(meta_df$sample_type == "Solid Tissue Normal")
    
    project_status_df <- rbind(project_status_df, data.frame(
      cancer     = proj, 
      has_tumor  = contains_tumor, 
      has_normal = contains_normal
    ))
    
  }, error = function(e) {
    # If a project has no matching RNA-Seq data, mark both as FALSE
    project_status_df <- rbind(project_status_df, data.frame(
      cancer     = proj, 
      has_tumor  = FALSE, 
      has_normal = FALSE
    ))
  })
}

# --- EVALUATE AND FILTER ---

# 3. Print cohorts that DO NOT have normal samples
no_normal_cohorts <- project_status_df %>% 
  filter(has_tumor == TRUE & has_normal == FALSE) %>% 
  pull(cancer)

cat("\n--- Cancers that completely LACK solid tissue normal samples: ---\n")
print(no_normal_cohorts)

# 4. Remove these problem projects from your main tcga_projects vector
# This leaves you with a clean vector containing only cohorts with BOTH types
clean_tcga_projects <- project_status_df %>%
  filter(has_tumor == TRUE & has_normal == TRUE) %>%
  pull(cancer)

cat("\n--- Cleaned list of projects (Safe for differential expression): ---\n")
print(clean_tcga_projects)

# 5. Save the validation log to a CSV file for your project reference
write.csv(project_status_df, file = "TCGA_tumor_normal_checking_summary.csv", row.names = FALSE)
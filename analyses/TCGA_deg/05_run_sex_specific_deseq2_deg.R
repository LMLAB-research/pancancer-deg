# Run sex-specific DESeq2 models ---------------------------------------------
# Load configuration and shared helpers before starting the analysis.

source(file.path("analyses", "TCGA_deg", "00_config.R"))
source(file.path(analysis_dir, "utils.R"))

# Run all TCGA projects ------------------------------------------------------

run_sex_stratified <- function() {
  if (use_lfc_shrinkage && lfc_shrinkage_type != "ashr") {
    stop("Set lfc_shrinkage_type to ashr in 00_config.R")
  }
  
  plan <- readRDS(tcga_design_plan_rds)
  adjusted_projects <- plan$project_id[
    which(plan$eligible & plan$design_type == "adjusted")
  ]
  project_ids <- limit_projects(unique(plan$project_id))
  
  if (length(project_ids) == 0L) {
    stop("No TCGA projects available for analysis")
  }
  
  context <- create_run_context()
  outputs <- list()
  
  # Fit female, male and interaction models independently for each project.
  for (project_id in project_ids) {
    project <- tryCatch(
      load_sex_project(project_id, adjusted_projects),
      error = identity
    )
    
    for (analysis in c("female", "male", "interaction")) {
      message(project_id, " | ", analysis)
      result <- run_sex_model(project_id, analysis, project, context)
      outputs[[length(outputs) + 1L]] <- result
      
      # Preserve each analysis report, including skipped and failed models.
      audit_path <- file.path(
        context$table_dir,
        paste0(project_id, "_", analysis, "_audit_status.rds")
      )
      saveRDS(result, audit_path)
    }
    
    rm(project)
  }
  
  # Write run-wide reports before checking whether any model failed.
  status <- dplyr::bind_rows(lapply(outputs, `[[`, "status"))
  audit <- dplyr::bind_rows(lapply(outputs, `[[`, "audit"))
  
  write.csv(
    status,
    file.path(context$table_dir, "run_status.csv"),
    row.names = FALSE
  )
  write.csv(
    audit,
    file.path(context$table_dir, "sample_audit.csv"),
    row.names = FALSE
  )
  
  print(table(status$analysis_type, status$status))
  message("Reports: ", normalizePath(context$table_dir))
  
  failed <- status[status$status == "failed", , drop = FALSE]
  if (nrow(failed) > 0L) {
    print(failed[, c("project_id", "analysis_type", "reason"), drop = FALSE])
    stop(nrow(failed), " DESeq2 analyses failed. See run_status.csv")
  }
  
  invisible(outputs)
}

run_sex_stratified()

#### Run DESeq2 Differential Expression ####
#
# This script runs projects marked eligible in the design plan, writes raw and
# shrunken DEG tables, and creates first-pass QC plots.

source(file.path("analyses", "TCGA_deg", "00_config.R"))

source(file.path(analysis_dir, "utils.R"))

#### Select Projects From Design Plan ####

design_plan <- readRDS(tcga_design_plan_rds)
eligible_designs <- design_plan |>
  dplyr::filter(eligible)
eligible_projects <- limit_projects(unique(eligible_designs$project_id))
eligible_designs <- eligible_designs |>
  dplyr::filter(project_id %in% eligible_projects) |>
  dplyr::mutate(
    create_plots = design_type == "adjusted" |
      !project_id %in% project_id[design_type == "adjusted"]
  )

#### Run Each Eligible Project ####

run_status <- dplyr::bind_rows(Map(
  run_project_deseq2_safe,
  eligible_designs$project_id,
  eligible_designs$design_type,
  eligible_designs$create_plots
))
write.csv(run_status, tcga_deseq_run_status_file, row.names = FALSE)

#### Run DESeq2 Differential Expression ####
#
# This script runs every design marked eligible in the design plan (all
# samples, female, male, and sex interaction; naive and adjusted), writes raw
# and shrunken DEG tables, and creates first-pass QC plots.

source(file.path("analyses", "TCGA_deg", "00_config.R"))

source(file.path(analysis_dir, "utils.R"))

#### Select Projects From Design Plan ####

design_plan <- readRDS(tcga_design_plan_rds)
eligible_designs <- design_plan |>
  dplyr::filter(eligible)
eligible_projects <- limit_projects(unique(eligible_designs$project_id))
eligible_designs <- eligible_designs |>
  dplyr::filter(project_id %in% eligible_projects) |>
  # QC plots are made once per project, for the all-samples adjusted design
  # (or the naive one if no adjusted design is eligible).
  dplyr::mutate(
    create_plots = design_type == "adjusted" |
      (design_type == "naive" &
        !project_id %in% project_id[design_type == "adjusted"])
  )

#### Run Each Eligible Project ####

run_status <- dplyr::bind_rows(Map(
  run_project_deseq2_safe,
  eligible_designs$project_id,
  eligible_designs$design_type,
  eligible_designs$create_plots
))
write.csv(run_status, tcga_deseq_run_status_file, row.names = FALSE)

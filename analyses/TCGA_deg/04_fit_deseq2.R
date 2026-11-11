#### Fit DESeq2 Models ####
#
# Fits every design marked eligible in the design plan (all samples, female,
# male and sex interaction; naive, adjusted and paired) and saves the fitted
# objects. This is the slow step. Testing is done in step 05, from the fitted
# objects, so that test settings can change without refitting.

source(file.path("analyses", "TCGA_deg", "00_config.R"))

source(file.path(analysis_dir, "utils.R"))

#### Select Designs From Design Plan ####

design_plan <- readRDS(tcga_design_plan_rds)
eligible_designs <- design_plan |>
  dplyr::filter(eligible)
eligible_projects <- limit_projects(unique(eligible_designs$project_id))
eligible_designs <- eligible_designs |>
  dplyr::filter(project_id %in% eligible_projects) |>
  # PCA plots are made once per project, for the all-samples adjusted design
  # (or the naive one if no adjusted design is eligible).
  dplyr::mutate(
    create_plots = design_type == "adjusted" |
      (design_type == "naive" &
        !project_id %in% project_id[design_type == "adjusted"])
  )

#### Fit Each Eligible Design ####

fit_status <- dplyr::bind_rows(Map(
  function(project_id, design_type, create_plots) {
    run_safely(fit_project_deseq2, project_id, design_type, create_plots)
  },
  eligible_designs$project_id,
  eligible_designs$design_type,
  eligible_designs$create_plots
))
write.csv(fit_status, tcga_deseq_fit_status_file, row.names = FALSE)

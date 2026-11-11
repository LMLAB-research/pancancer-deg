#### Test DESeq2 Models ####
#
# Tests the coefficient of every fitted design (step 04) with the test set by
# deseq_test in 00_config.R, shrinks log2 fold changes, and writes raw and
# shrunken result tables, MA and volcano plots, and a test status. Fast: rerun
# from here after changing deseq_test, deseq_alpha or deseq_lfc_threshold.

source(file.path("analyses", "TCGA_deg", "00_config.R"))

source(file.path(analysis_dir, "utils.R"))

#### Select Designs From Design Plan ####

design_plan <- readRDS(tcga_design_plan_rds)
eligible_designs <- design_plan |>
  dplyr::filter(eligible)
eligible_projects <- limit_projects(unique(eligible_designs$project_id))
eligible_designs <- eligible_designs |>
  dplyr::filter(project_id %in% eligible_projects) |>
  # MA and volcano plots for the same design as the PCA plots in step 04.
  dplyr::mutate(
    create_plots = design_type == "adjusted" |
      (design_type == "naive" &
        !project_id %in% project_id[design_type == "adjusted"])
  )

#### Test Each Fitted Design ####

test_status <- dplyr::bind_rows(Map(
  function(project_id, design_type, create_plots) {
    run_safely(test_project_deseq2, project_id, design_type, create_plots)
  },
  eligible_designs$project_id,
  eligible_designs$design_type,
  eligible_designs$create_plots
))
write.csv(test_status, tcga_deseq_test_status_file, row.names = FALSE)

#### Data and results locations ####
#
# All modules read and write under two roots. By default they are the `data/`
# and `results/` folders of the repository (scripts are run from its root).
# To keep data and results elsewhere on a machine, set the environment
# variables, e.g. in ~/.bashrc (for run_pipeline.sh) and ~/.Renviron (for
# interactive R sessions):
#
#   export PANCANCER_DATA_DIR=/path/to/pancancer-deg/data
#   export PANCANCER_RESULTS_DIR=/path/to/pancancer-deg/results

data_root <- Sys.getenv("PANCANCER_DATA_DIR", unset = "data")
results_root <- Sys.getenv("PANCANCER_RESULTS_DIR", unset = "results")

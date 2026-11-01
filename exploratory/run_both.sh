#!/usr/bin/env bash
# Run the two-sided overlap pipeline (worktree code) on the main and the
# proliferation-free signature inputs. Run from the worktree root.
set -u
T=/tmp/claude-501/two
mkdir -p "$T/main_in"
ln -sfn /Users/ro/code/academic/pancancer-deg/data/senescence_meta_analysis_deg "$T/main_in/meta_analysis"
ln -sfn /Users/ro/code/academic/pancancer-deg/results/tables/TCGA_deg "$T/main_in/TCGA_deg"
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
for run in main sens; do
  if [ "$run" = main ]; then B="$T/main_in"; else B=/tmp/claude-501/sens; fi
  for s in TCGA_senescence_overlaps.R heatmap_final.R; do
    TCGA_SENESCENCE_BASE_DIR="$B" TCGA_SENESCENCE_RESULTS_DIR="$T/results_$run" \
      Rscript "TCGA_senescence_overlapping_pipeline/$s" >> "$T/$run.log" 2>&1 \
      || { echo "FAIL $run $s"; tail -15 "$T/$run.log"; exit 1; }
  done
  echo "$run done: $(grep -c COMPLETED "$T/results_$run/batch_run_status.csv") completed"
done

#!/usr/bin/env bash
# Run pipeline steps in order and write one log per step to results/logs/.
#
# Usage (from anywhere in the repository):
#   bash analyses/run_pipeline.sh STEP [STEP ...]
#
# Steps:
#   deg01 deg02 deg03 deg04    analyses/TCGA_deg/01-04
#   overlap01 overlap02        analyses/TCGA_senescence_overlap/01-02
# Groups:
#   deg = deg01-deg04, overlap = overlap01-overlap02, all = deg + overlap
#   from:STEP = STEP and every step after it, e.g. from:deg04
#
# Examples:
#   bash analyses/run_pipeline.sh all
#   bash analyses/run_pipeline.sh from:deg04
#   bash analyses/run_pipeline.sh overlap
#   nohup bash analyses/run_pipeline.sh all > /dev/null 2>&1 &
#
# Steps always run in pipeline order, whatever order they are given in. The run
# stops at the first failing step. Logs: results/logs/<run time>_<step>.log, and
# a summary of the whole run in results/logs/<run time>_run.log.

set -euo pipefail

steps=(deg01 deg02 deg03 deg04 overlap01 overlap02)
declare -A scripts=(
  [deg01]="analyses/TCGA_deg/01_discover_tcga_rnaseq_projects.R"
  [deg02]="analyses/TCGA_deg/02_prepare_tcga_inputs.R"
  [deg03]="analyses/TCGA_deg/03_build_deseq2_designs.R"
  [deg04]="analyses/TCGA_deg/04_run_deseq2_deg.R"
  [overlap01]="analyses/TCGA_senescence_overlap/01_compute_overlaps.R"
  [overlap02]="analyses/TCGA_senescence_overlap/02_plot_heatmaps.R"
)

usage() {
  sed -n '2,22p' "$0" | sed 's/^# \{0,1\}//'
  exit 1
}

[[ $# -eq 0 ]] && usage

# Expand the arguments into a set of selected steps.
declare -A selected=()
for arg in "$@"; do
  case "$arg" in
    all) for s in "${steps[@]}"; do selected[$s]=1; done ;;
    deg) for s in deg01 deg02 deg03 deg04; do selected[$s]=1; done ;;
    overlap) for s in overlap01 overlap02; do selected[$s]=1; done ;;
    from:*)
      start="${arg#from:}"
      [[ -n "${scripts[$start]:-}" ]] || { echo "Unknown step: $start" >&2; usage; }
      on=0
      for s in "${steps[@]}"; do
        [[ "$s" == "$start" ]] && on=1
        [[ $on -eq 1 ]] && selected[$s]=1
      done
      ;;
    -h|--help) usage ;;
    *)
      [[ -n "${scripts[$arg]:-}" ]] || { echo "Unknown step: $arg" >&2; usage; }
      selected[$arg]=1
      ;;
  esac
done

# Scripts expect the repository root as working directory.
cd "$(dirname "$0")/.."
log_dir="results/logs"
mkdir -p "$log_dir"
run_time="$(date +%Y-%m-%d_%H%M%S)"
run_log="$log_dir/${run_time}_run.log"

say() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$run_log"; }

say "Run started; steps: $(for s in "${steps[@]}"; do [[ -n "${selected[$s]:-}" ]] && printf '%s ' "$s"; done)"
say "Commit: $(git rev-parse --short HEAD 2>/dev/null || echo unknown)$(git diff --quiet 2>/dev/null || echo ' (uncommitted changes)')"

for s in "${steps[@]}"; do
  [[ -n "${selected[$s]:-}" ]] || continue
  step_log="$log_dir/${run_time}_${s}.log"
  say "START $s: ${scripts[$s]} (log: $step_log)"
  start_seconds=$SECONDS
  if Rscript "${scripts[$s]}" > "$step_log" 2>&1; then
    say "DONE  $s in $(( (SECONDS - start_seconds) / 60 )) min"
  else
    say "FAILED $s after $(( (SECONDS - start_seconds) / 60 )) min; last lines of $step_log:"
    tail -n 20 "$step_log" | tee -a "$run_log"
    exit 1
  fi
done

say "Run finished"

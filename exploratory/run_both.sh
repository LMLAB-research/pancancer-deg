#!/usr/bin/env bash
# Run the two-sided overlap pipeline on the meta-analysis signatures ("full_signatures") and on the same signatures
# with the 671 Hallmark proliferation genes (exploratory/inputs/other/prolif_genes.txt) removed from every genes and
# background file ("no_proliferation_genes"; DECISIONS B5). Run from the pancancer-deg repo root:
#   bash exploratory/run_both.sh
# Inputs: data/senescence_meta_analysis_deg if present, else exploratory/inputs/senescence_meta_analysis_deg;
#   TCGA DEGs from results/tables/TCGA_deg if present, else the gzipped copies in exploratory/inputs/tcga_degs.
# Output: ${RUN_BOTH_OUT:-results/exploratory/two_sided_rerun}/{full_signatures,no_proliferation_genes}
set -euo pipefail
[ -d exploratory/inputs ] || { echo "run from the pancancer-deg repo root"; exit 1; }
OUT=${RUN_BOTH_OUT:-results/exploratory/two_sided_rerun}
mkdir -p "$OUT"; T="$(cd "$OUT" && pwd)/.inputs_tmp"; rm -rf "$T"   # working copies of the inputs (absolute; removed at the end)
META=data/senescence_meta_analysis_deg; [ -d "$META" ] || META=exploratory/inputs/senescence_meta_analysis_deg
mkdir -p "$OUT" "$T/main/TCGA_deg" "$T/sens/meta_analysis"
cp "$META"/*.csv "$T/sens/meta_analysis/"; ln -s "$(cd "$META" && pwd)" "$T/main/meta_analysis"
if ls results/tables/TCGA_deg/*_DGE_Raw_Results_Tumor_vs_Normal_Adjusted.csv >/dev/null 2>&1; then
  ln -s "$(pwd)"/results/tables/TCGA_deg/*_DGE_Raw_Results_Tumor_vs_Normal_Adjusted.csv "$T/main/TCGA_deg/"
else
  for f in exploratory/inputs/tcga_degs/*.csv.gz; do gunzip -c "$f" > "$T/main/TCGA_deg/$(basename "${f%.gz}")"; done
fi
mkdir -p "$T/sens/TCGA_deg"   # proliferation-free TCGA DEG tables (gene_name, quoted, is column 8)
for f in "$T/main/TCGA_deg"/*.csv; do
  awk -F, 'NR == FNR { p[$1] = 1; next } { g = $8; gsub(/"/, "", g) } FNR == 1 || !(g in p)' \
    exploratory/inputs/other/prolif_genes.txt "$f" > "$T/sens/TCGA_deg/$(basename "$f")"
done
# proliferation-free signatures: drop the proliferation genes (first column) from every genes / background file
for f in "$T/sens/meta_analysis"/*.csv; do
  awk -F, 'NR == FNR { p[$1] = 1; next } FNR == 1 || !($1 in p)' \
    exploratory/inputs/other/prolif_genes.txt "$f" > "$f.tmp" && mv "$f.tmp" "$f"
done
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
for run in ${RUN_BOTH_RUNS:-main sens}; do
  R="$OUT/$([ "$run" = main ] && echo full_signatures || echo no_proliferation_genes)"
  for s in TCGA_senescence_overlaps.R heatmap_final.R; do
    TCGA_SENESCENCE_BASE_DIR="$T/$run" TCGA_SENESCENCE_RESULTS_DIR="$R" \
      Rscript "TCGA_senescence_overlapping_pipeline/$s" >> "$OUT/$run.log" 2>&1 \
      || { echo "FAIL $run $s"; tail -15 "$OUT/$run.log"; exit 1; }
  done
  echo "$run done: $(grep -c COMPLETED "$R/batch_run_status.csv") completed"
done
rm -rf "$T"

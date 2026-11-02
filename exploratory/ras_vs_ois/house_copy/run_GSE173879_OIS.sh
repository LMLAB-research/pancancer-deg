#!/bin/zsh
# Stage GSE173879 OIS arms (local copy of stage_inputs.R) and run the UNCHANGED house one-plate preprocessing
# (senescence_stress house_rerun/scripts/scRNA_seq_preprocess_final_oneplate.R) with the house GSE173879 QC
# (Wechter thresholds, as in run_house_preprocess.sh). All outputs in pancancer-deg/data/exploratory/house_GSE173879_OIS/.
source /Users/ro/Academic/10_Research/senescence_stress/modules_intact/scrna_seq/GSE226225/final_analysis/env/activate.sh
H=/Users/ro/Academic/10_Research/senescence_stress/modules/scrna_seq/cs_dataset_survey/house_rerun
B=${0:A:h}
D=/Users/ro/code/academic/pancancer-deg/data/exploratory/house_GSE173879_OIS
mkdir -p $D/logs
Rscript $B/stage_inputs_GSE173879_OIS.R GSE173879 > $D/logs/stage.log 2>&1 || { echo "staging failed"; tail -20 $D/logs/stage.log; exit 1; }
export DATA_DIR=$D/staged/ OUT_DIR=$D/ INTEGRITY_DIR=$D/integrity/
export BATCHES=OISON_1,OISOFF_1 QC_MT_MAX=12 QC_UMI_MIN=1000 QC_UMI_MAX=120000 QC_GENES_MIN=300 QC_GENES_MAX=Inf
env | grep -E "^(DATA_DIR|OUT_DIR|INTEGRITY_DIR|BATCHES|QC_|MODULES_ROOT)" | sort > $D/run_env.txt
cd $D && Rscript $H/scripts/scRNA_seq_preprocess_final_oneplate.R > $D/logs/preprocess.log 2>&1
echo "preprocess exit $?"

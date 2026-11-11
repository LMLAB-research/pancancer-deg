# pancancer-deg

pancancer-deg contains analysis code for the SenescenceCancer pan-cancer differential expression project. The current codebase focuses on TCGA differential expression analysis, with the repository layout prepared for additional analyses later in the project.

## R setup

This project uses `renv` for reproducible R package management.

```bash
bash setup/setup.sh
```

The setup script installs `renv` if needed, restores packages from `renv.lock`,
and runs a package availability check. To run only the check:

```bash
RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript --vanilla setup/check_setup.R
```

Setup never changes `renv.lock`. When a package is added to or removed from
the code, list it in `setup/install_dependencies.R` and `setup/check_setup.R`,
then update the lockfile and commit it:

```bash
RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript --vanilla setup/install_dependencies.R --snapshot
```

## Project Structure

```text
pancancer-deg/
|-- analyses/
|   |-- TCGA_deg/
|   `-- TCGA_senescence_overlap/
|-- data/
|-- results/
|   |-- figures/
|   |   |-- TCGA_deg/
|   |   `-- TCGA_senescence_overlap/
|   `-- tables/
|       |-- TCGA_deg/
|       `-- TCGA_senescence_overlap/
|-- setup/
|   |-- check_setup.R
|   |-- install_dependencies.R
|   `-- setup.sh
|-- README.md
|-- renv.lock
`-- pancancer_deg.Rproj
```

The `analyses/` directory contains analysis scripts grouped by analysis module. Current modules are `TCGA_deg/` and `TCGA_senescence_overlap/`; future modules can be added as sibling directories under `analyses/`.

The `results/` directory is organized by output type first, then by analysis module. For example, the TCGA differential expression analysis writes figures to `results/figures/TCGA_deg/` and tables to `results/tables/TCGA_deg/`.

The `data/` directory is intended for local input data, downloaded resources, and intermediate analysis objects. Large or generated data files should generally stay out of Git.

## Running the pipeline

`analyses/run_pipeline.sh` runs selected steps in order, stops at the first
failure, and writes one log per step to `results/logs/`:

```bash
bash analyses/run_pipeline.sh all           # everything
bash analyses/run_pipeline.sh from:deg04    # deg04 and all later steps
bash analyses/run_pipeline.sh overlap       # overlap01 and overlap02
nohup bash analyses/run_pipeline.sh all > /dev/null 2>&1 &   # keep running after logout
```

Steps: `deg01`–`deg05` (TCGA_deg) and `overlap01`–`overlap02`
(TCGA_senescence_overlap); groups `deg`, `overlap`, `all`. Run
`bash analyses/run_pipeline.sh` without arguments for help.

## TCGA_deg Analysis

Tumor-versus-normal differential expression (DESeq2) for each TCGA cancer type.
Settings are in `analyses/TCGA_deg/00_config.R`; shared functions in `utils.R`.
Run the scripts in order from the repository root:

```bash
Rscript analyses/TCGA_deg/01_discover_tcga_rnaseq_projects.R
Rscript analyses/TCGA_deg/02_prepare_tcga_inputs.R
Rscript analyses/TCGA_deg/03_build_deseq2_designs.R
Rscript analyses/TCGA_deg/04_run_deseq2_deg.R
Rscript analyses/TCGA_deg/05_compare_designs.R
```

| Step | What it does |
| --- | --- |
| 01 | Lists TCGA projects with at least 5 normal and 5 tumor RNA-seq samples. |
| 02 | Downloads STAR counts and clinical data. Tumor = primary tumor, normal = solid tissue normal; other sample types are excluded. Sums replicate files of the same patient and condition (below). |
| 03 | Builds every DESeq2 design (below) and writes the design plan for review. |
| 04 | Fits every eligible design and writes raw and ashr-shrunken result tables. |
| 05 | Compares naive, adjusted and paired results per cancer and sex model (fold-change correlation, DEG agreement). |

### Replicate files

Some patients have several RNA-seq files for the same condition: the same
sample sequenced twice (e.g. on two plates), or several vials or portions of
the same tissue. In 9 of the 18 eligible cancers this concerns 4 to 88 files
(88 in GBM). Treating them as separate samples would be pseudo-replication: one
tumor would count several times. Step 02 therefore sums their counts into one
sample per patient and condition, as DESeq2's `collapseReplicates` does for
technical replicates. Each patient still counts once, and the summed sample
has a higher sequencing depth, which slightly reduces counting noise for
low-expressed genes. The merged sample keeps the metadata of its first barcode
(vial A before B); `n_files_merged` and `merged_barcodes` record the merge.

### Designs

Each cancer is analysed with four models:

| Model | Samples | Tested coefficient |
| --- | --- | --- |
| all | all | tumor effect $\beta$ |
| female | women only | tumor effect in women $\beta^F$ |
| male | men only | tumor effect in men $\beta^M$ |
| sex interaction | both sexes | $\delta = \beta^M - \beta^F$ |

All effects are log2 fold changes (tumor vs normal).

#### Three tests of the tumor effect

For robustness, each model is fitted in three versions that differ in how
differences between patients are handled:

| Version | Samples | Formula (all-samples model) | Controls for |
| --- | --- | --- | --- |
| naive | all tumors and normals | `~ condition` | nothing |
| adjusted | complete cases on the selected covariates | `~ age + sex + race + ... + condition` | recorded clinical covariates |
| paired | patients with both a tumor and a normal sample | `~ patient_id + condition` | every patient-level factor |

- **Naive** compares tumors and normals from different people. It uses all
  samples but mixes the tumor effect with any patient difference between the
  two groups (normals come from only about 10% of patients).
- **Adjusted** adds clinical covariates (age, sex, race, smoking status),
  selected per design: a covariate is kept if at most 25% of values are missing
  and its smallest group has at least 3 samples (sex only in the all-samples
  model). It corrects only for recorded factors and drops samples with missing
  values, which in some cancers removes 16–34% of the normals.
- **Paired** compares tumor and normal within each patient, so every factor
  that is constant within a patient cancels out, recorded or not (age, sex,
  ancestry, smoking, environment, ...), and no clinical data are needed. About
  97% of normals have a matched tumor. For the interaction, sex is constant
  within patients, so the formula is `~ patient_id + condition + tumor_male`,
  where `tumor_male` = 1 for tumor samples of men; its coefficient is $\delta$.

Conceptually, the paired version estimates the quantity of interest, the
change from normal to tumor tissue within a person, and matches the design of
the senescence signatures (senescent vs proliferating cells of the same
culture). **We expect to keep the paired version as the main analysis**, with
naive as a sensitivity check covering all tumors. The adjusted version is kept
for now. None of the three corrects for tumor purity or for adjacent normal
tissue not being fully healthy.

Step 05 compares the three versions per cancer and sex model. **Status:** the
full chain (steps 03–05) has been tested on one cancer (TCGA-KICH), where the
versions largely agree (all-samples model: Spearman correlation of fold
changes 0.95–0.99, Jaccard index of DEG sets 0.80–0.89). Results for all
cancers are pending a full rerun.

A design is eligible if it has at least 5 normal and 5 tumor samples after
filtering (in each sex for the interaction model). Reasons for ineligible
designs are in `03_TCGA_deseq2_design_plan.csv`.

To find genes whose tumor effect differs between sexes, use the interaction
model ($\delta \neq 0$), not a comparison of the female and male gene lists: a
gene significant in one sex and not in the other does not show that the two
effects differ.

### Outputs

In `results/tables/TCGA_deg/`, per cancer and design:

- `<project>_DGE_Raw_Results_Tumor_vs_Normal_<Naive|Adjusted|Paired>[_Female|_Male].csv`
- `<project>_DGE_Raw_Results_Sex_Interaction_<Naive|Adjusted|Paired>.csv`
- the same names with `DGE_Results` for shrunken log2 fold changes.
- `05_TCGA_design_comparison.csv`: agreement between the naive, adjusted and
  paired versions (step 05).

QC figures (PCA, MA, volcano) are made for the all-samples model only, in
`results/figures/TCGA_deg/`.

## TCGA_senescence_overlap Analysis

Tests whether tumor DEGs (from `TCGA_deg`, all-samples adjusted design) overlap
senescence (OIS, RS, SIPS) and quiescence (CICQ, SSCQ) signatures. Settings are
in `analyses/TCGA_senescence_overlap/00_config.R`. Run from the repository root:

```bash
Rscript analyses/TCGA_senescence_overlap/01_compute_overlaps.R
Rscript analyses/TCGA_senescence_overlap/02_plot_heatmaps.R
```

Script 02 only redraws the heatmaps from the outputs of script 01.

### Inputs

- `data/senescence_meta_analysis_deg/`: one `<group>_<signature>_genes.csv`
  (columns `gene`, `log2FoldChange`, `padj`) and one
  `<group>_<signature>_background.csv` (all genes tested in that signature,
  column `gene`, `gene_symbol` or `gene_name`) per signature, e.g.
  `senescence_OIS_genes.csv`, `quiescence_SSCQ_background.csv`.
- `results/tables/TCGA_deg/<project>_DGE_Raw_Results_Tumor_vs_Normal_Adjusted.csv`.

### Method

A gene is up (down) if $\text{padj} < 0.05$ and $\log_2\text{FC} > \log_2 1.5$
($< -\log_2 1.5$). For each cancer and signature this gives sets $T^\pm$ and
$S^\pm$. The universe $U$ is the TCGA protein-coding genes with a padj,
intersected with the signature background; all sets are restricted to $U$.

Three readouts are computed, each with Fisher's exact test on a 2x2 table:

1. **Class overlaps** (reference). For each class $T^a \cap S^b$,
   $a, b \in \{+, -\}$ (UP-UP, UP-DOWN, DOWN-UP, DOWN-DOWN; cancer first), the
   table is "in $T^a$ or not" x "in $S^b$ or not" over $U$. Two-sided by default
   (`fisher_alternative`). BH across all cancer x class tests within each
   signature.
2. **Sharing.** Over $U$: "DE in cancer or not" x "DE in signature or not".
   Odds ratio > 1: DE genes are shared more than by chance.
3. **Direction.** Over genes DE in both: "up or down in cancer" x "up or down
   in signature". Odds ratio > 1: same direction; < 1: opposite direction.

Readouts 2 and 3 are always two-sided, with BH across cancers within each
signature.

How they relate: the class-overlap enrichment mixes sharing and direction. For
UP-UP, with ratios of fractions,

$$\frac{P(S^+ \mid T^+)}{P(S^+ \mid \text{not } T^+)}
= \frac{P(S\ \text{DE} \mid T^+)}{P(S\ \text{DE} \mid \text{not } T^+)}
\times \frac{P(S^+ \mid T^+, S\ \text{DE})}{P(S^+ \mid \text{not } T^+, S\ \text{DE})}.$$

Because sharing is strong, all four classes can be enriched at once even when
the shared genes move in opposite directions. Use the direction readout to
compare directions.

Duplicate gene symbols with opposite significant directions are excluded; other
duplicates keep one representative and are listed for review.

### Outputs

- `results/tables/TCGA_senescence_overlap/<signature>/<project>/`: shared gene
  lists, class-overlap statistics (`overlap_statistics.csv`), sharing and
  direction statistics (`share_direction_statistics.csv`) and the exact
  universe (`common_tested_genes.csv`).
- `results/tables/TCGA_senescence_overlap/batch_run_status.csv`: one row per
  cancer x signature; check it for failures.
- `results/figures/TCGA_senescence_overlap/TCGA_<signature>_heatmap.png`: log2
  odds ratios of significant class overlaps (blue = enrichment, red = depletion).
- `results/figures/TCGA_senescence_overlap/TCGA_share_heatmap.png` and
  `TCGA_direction_heatmap.png`: cancers x signatures (direction: blue = same,
  red = opposite).
- The plotted values of each heatmap are in
  `results/tables/TCGA_senescence_overlap/TCGA_*_heatmap_data.csv`.

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
bash analyses/run_pipeline.sh from:deg05    # new test settings: deg05 and all later steps, no refit
bash analyses/run_pipeline.sh overlap       # overlap01 and overlap02
nohup bash analyses/run_pipeline.sh all > /dev/null 2>&1 &   # keep running after logout
```

Steps: `deg01`–`deg06` (TCGA_deg) and `overlap01`–`overlap02`
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
Rscript analyses/TCGA_deg/04_fit_deseq2.R
Rscript analyses/TCGA_deg/05_test_deseq2.R
Rscript analyses/TCGA_deg/06_compare_designs.R
```

| Step | What it does |
| --- | --- |
| 01 | Lists TCGA projects with at least 5 normal and 5 tumor RNA-seq samples. |
| 02 | Downloads STAR counts and clinical data. Tumor = primary tumor, normal = solid tissue normal; other sample types are excluded. Sums replicate files of the same patient and condition (below). |
| 03 | Builds every DESeq2 design (below) and writes the design plan for review. |
| 04 | Fits every eligible design (the slow step) and saves the fitted DESeq2 objects. |
| 05 | Tests each fitted design (see *Differential expression test*), shrinks log2 fold changes (ashr), and writes raw and shrunken result tables. |
| 06 | Compares naive, adjusted and paired results per cancer and sex model (fold-change correlation, DEG agreement), and summarises the sex models per cancer and version. |

Each step saves its outputs, so later changes need only a partial rerun:

| To change | Rerun from |
| --- | --- |
| test, FDR level or fold-change threshold | `deg05` (minutes; no refit) |
| designs or covariate rules | `deg03` |
| sample selection or clinical merging | `deg02` (downloads are reused) |

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

Step 06 compares the three versions per cancer and sex model. **Status:** the
full chain (steps 03–06) has been tested on one cancer (TCGA-KICH), where the
versions largely agree (all-samples model, threshold test: Spearman
correlation of fold changes 0.95–0.99, Jaccard index of DEG sets 0.72–0.85,
at most 2 DEGs with opposite directions). Results for all cancers are pending
a full rerun.

A design is eligible if it has at least 5 normal and 5 tumor samples after
filtering (in each sex for the interaction model). Reasons for ineligible
designs are in `03_TCGA_deseq2_design_plan.csv`.

To find genes whose tumor effect differs between sexes, use the interaction
model ($\delta \neq 0$), not a comparison of the female and male gene lists: a
gene significant in one sex and not in the other does not show that the two
effects differ.

### Differential expression test

By default (`deseq_test <- "threshold"` in `00_config.R`), step 05 tests
$H_0: |\beta| \le \log_2 1.5$ against $|\beta| > \log_2 1.5$ (DESeq2
`lfcThreshold` with `altHypothesis = "greaterAbs"`). A gene is a DEG if its
BH-adjusted p-value is below 0.05; its direction is the sign of $\hat\beta$.
The false discovery rate therefore applies to changes larger than 1.5-fold.

The alternative, `deseq_test <- "zero"`, tests $H_0: \beta = 0$ and keeps
genes with $|\hat\beta| > \log_2 1.5$ as a separate filter (the previous
approach). With hundreds of samples nearly every expressed gene rejects
$\beta = 0$, and the p-value does not account for the fold-change filter; in
TCGA-KICH (paired design) the threshold test gives 6,607 DEGs instead of
10,643.

In both cases DEG calls use the raw (unshrunken) $\hat\beta$ and are written
in the `deg` column (`up`, `down`, `not_significant`) of every result table;
steps 06 and the overlap analysis use this column. Shrunken log2 fold changes
are for ranking and plotting.

### Outputs

In `results/tables/TCGA_deg/`, per cancer and design:

- `<project>_DGE_Raw_Results_Tumor_vs_Normal_<Naive|Adjusted|Paired>[_Female|_Male].csv`
- `<project>_DGE_Raw_Results_Sex_Interaction_<Naive|Adjusted|Paired>.csv`
- the same names with `DGE_Results` for shrunken log2 fold changes.
- `03_TCGA_deseq2_design_plan.csv`: every design, eligible or not, with
  sample counts and formula.
- `04_TCGA_deseq2_fit_status.csv` and `05_TCGA_deseq2_test_status.csv`: one
  row per design (status, test used, DEG counts).
- `06_TCGA_design_comparison.csv`: agreement between the naive, adjusted and
  paired versions.
- `06_TCGA_sex_model_summary.csv`: per cancer and version, sample sizes of the
  female and male models, DEG counts of all four models, and the agreement of
  female and male tumor effects (Spearman correlation of log2 fold changes;
  DEGs in both sexes with the same or opposite direction).

Fitted DESeq2 objects are in `data/TCGA_deg/processed/*_fit_dds.rds`.

QC figures (PCA, MA, volcano) are made for the all-samples model only, in
`results/figures/TCGA_deg/`.

### Limitations

- **Adjacent normal tissue.** TCGA normals are taken next to the tumor and can
  already carry tumor-related changes (field effects); they are not healthy
  tissue. Tumor effects may therefore be underestimated.
- **Tumor purity and cell composition.** Tumor and normal samples differ in
  their mix of cell types (tumor cells, stroma, immune cells). Part of each
  tumor effect reflects this mix rather than changes within cells. No design
  corrects for it.
- **Batch effects.** Samples were processed on different plates and in
  different centers. The paired design removes differences between patients
  but not between the batches of a patient's tumor and normal samples.
- **Different power across cancers.** Sample sizes range from 8 to 113 pairs,
  so the number of DEGs partly reflects sample size. Compare cancers by effect
  sizes (e.g. fold-change correlations) rather than by DEG counts alone.
- **GBM.** Its 5 normals have no clinical data and no matched tumor, so GBM
  has only a naive design, based on 5 normals of unknown background.

### To consider

- **False discoveries with large samples.** DESeq2 and edgeR can show inflated
  false discovery rates in large population datasets such as TCGA tumor vs
  normal (Li et al., *Genome Biology* 2022). Planned check: permute the
  tumor/normal labels (within patients for the paired design) in one large
  cancer and count DEGs, which should be close to 0 (issue #29).
- **Tumor purity.** Use purity estimates (e.g. ABSOLUTE, ESTIMATE) as a
  covariate or in a sensitivity analysis (issue #30).

## TCGA_senescence_overlap Analysis

Tests whether tumor DEGs overlap senescence (OIS, RS, SIPS) and quiescence
(CICQ, SSCQ) signatures. Tumor DEGs are the `deg` column of one all-samples
`TCGA_deg` design, set by `cancer_design_type` in
`analyses/TCGA_senescence_overlap/00_config.R`: `"paired"` by default (the
main analysis), or `"naive"` / `"adjusted"`. Run from the repository root:

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
- `results/tables/TCGA_deg/<project>_DGE_Raw_Results_Tumor_vs_Normal_<Paired|Naive|Adjusted>.csv`
  (the design set by `cancer_design_type`) and `05_TCGA_deseq2_test_status.csv`.

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
  universe (`common_tested_genes.csv`). The statistics tables also give the
  cancer's sample sizes (`n_normal`, `n_tumor`, from TCGA_deg step 05) and the
  numbers of cancer and signature DEGs in the universe, so that results can be
  interpreted against power (cancers with more samples have more DEGs).
- `results/tables/TCGA_senescence_overlap/batch_run_status.csv`: one row per
  cancer x signature; check it for failures.
- `results/figures/TCGA_senescence_overlap/TCGA_<signature>_heatmap.png`: log2
  odds ratios of significant class overlaps (blue = enrichment, red = depletion).
- `results/figures/TCGA_senescence_overlap/TCGA_share_heatmap.png` and
  `TCGA_direction_heatmap.png`: cancers x signatures (direction: blue = same,
  red = opposite).
- The plotted values of each heatmap are in
  `results/tables/TCGA_senescence_overlap/TCGA_*_heatmap_data.csv`.

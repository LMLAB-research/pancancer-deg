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

## TCGA_deg Analysis

Tumor-versus-normal differential expression (DESeq2) for each TCGA cancer type.
Settings are in `analyses/TCGA_deg/00_config.R`; shared functions in `utils.R`.
Run the scripts in order from the repository root:

```bash
Rscript analyses/TCGA_deg/01_discover_tcga_rnaseq_projects.R
Rscript analyses/TCGA_deg/02_prepare_tcga_inputs.R
Rscript analyses/TCGA_deg/03_build_deseq2_designs.R
Rscript analyses/TCGA_deg/04_run_deseq2_deg.R
```

| Step | What it does |
| --- | --- |
| 01 | Lists TCGA projects with at least 5 normal and 5 tumor RNA-seq samples. |
| 02 | Downloads STAR counts and clinical data. Tumor = primary tumor, normal = solid tissue normal; other sample types are excluded. |
| 03 | Builds every DESeq2 design (below) and writes the design plan for review. |
| 04 | Fits every eligible design and writes raw and ashr-shrunken result tables. |

### Designs

Each cancer gets a naive and an adjusted version of four models:

| Model | Samples | Formula (naive) | Tested coefficient |
| --- | --- | --- | --- |
| all | all | `~ condition` | tumor effect $\beta$ |
| female | women only | `~ condition` | tumor effect in women $\beta^F$ |
| male | men only | `~ condition` | tumor effect in men $\beta^M$ |
| sex interaction | both sexes | `~ sex * condition` | $\delta = \beta^M - \beta^F$ |

All effects are log2 fold changes (tumor vs normal). The adjusted version adds
covariates in front of these terms, e.g. `~ age + race + condition`. Candidate
covariates (age, sex, race, smoking status) are selected separately for each
design: a covariate is kept if at most 25% of values are missing and its
smallest group has at least 3 samples. Sex is a candidate only in the
all-samples model.

A design is eligible if it has at least 5 normal and 5 tumor samples after
filtering (in each sex for the interaction model). Reasons for ineligible
designs are in `03_TCGA_deseq2_design_plan.csv`.

To find genes whose tumor effect differs between sexes, use the interaction
model ($\delta \neq 0$), not a comparison of the female and male gene lists: a
gene significant in one sex and not in the other does not show that the two
effects differ.

### Outputs

In `results/tables/TCGA_deg/`, per cancer and design:

- `<project>_DGE_Raw_Results_Tumor_vs_Normal_<Naive|Adjusted>[_Female|_Male].csv`
- `<project>_DGE_Raw_Results_Sex_Interaction_<Naive|Adjusted>.csv`
- the same names with `DGE_Results` for shrunken log2 fold changes.

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

Script 02 only redraws the heatmaps from the outputs of script 01. Both need the
Bioconductor package `GeneOverlap`.

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

Four overlaps are tested, $T^a \cap S^b$ for $a, b \in \{+, -\}$ (UP-UP,
UP-DOWN, DOWN-UP, DOWN-DOWN; cancer first), with Fisher's exact test on the
2x2 table of $U$. The test is two-sided by default (`fisher_alternative`), so
both enrichment (odds ratio > 1) and depletion (< 1) are detected. P-values are
BH-adjusted across all cancer x class tests within each signature.

Duplicate gene symbols with opposite significant directions are excluded; other
duplicates keep one representative and are listed for review.

### Outputs

- `results/tables/TCGA_senescence_overlap/<signature>/<project>/`: shared gene
  lists, overlap statistics and the exact universe (`common_tested_genes.csv`).
- `results/tables/TCGA_senescence_overlap/batch_run_status.csv`: one row per
  cancer x signature; check it for failures.
- `results/figures/TCGA_senescence_overlap/TCGA_<signature>_heatmap.png`: log2
  odds ratios of significant overlaps (blue = enrichment, red = depletion);
  the plotted values are in `results/tables/TCGA_senescence_overlap/`.

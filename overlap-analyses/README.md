# Overlap pipeline

This pipeline compares two datasets of differentially expressed genes (DEGs) 
and creates a summary heatmap.

## Folder structure

```text
TCGA_OIS_overlapping_pipeline/
├── .gitignore
├── config.R
├── TCGA_OIS_overlaps.R
├── heatmap_final.R
├── README.md
├── data/       # local input files; not uploaded to GitHub
└── results/    # created automatically; not uploaded to GitHub
```

## Scripts

### `config.R`

Contains the DEG thresholds, filename pattern, cancer order, and heatmap labels.
By default, inputs are read from `data/` and outputs are written to `results/`.

### `TCGA_OIS_overlaps.R`

Main functions:

- `classify_regulation()` assigns `UP`, `DOWN`, or `NS` using the configured
  adjusted p-value and log2 fold-change thresholds.
- `validate_numeric()` stops the analysis if a numeric input column contains
  invalid text values.
- `resolve_symbol_duplicates()` resolves multiple rows mapped to the same gene
  symbol using the rule described below, and produces an audit table.
- `run_overlap_test()` calculates the observed and expected overlap, fold
  enrichment, odds ratio, and overlap p-value.
- `analyze_one_cancer()` prepares one TCGA dataset, runs all four directional
  overlap tests, and saves its result files.

#### Duplicate gene symbol rule

For every gene symbol, `resolve_symbol_duplicates()`:

1. Removes rows with a missing or empty gene symbol.
2. Groups all remaining rows by gene symbol.
3. Checks whether significant rows contain both `UP` and `DOWN` directions.
4. If both directions are present, marks the symbol as
   `CONFLICT_OPPOSITE_DIRECTIONS` and excludes it from the overlap analysis.
5. Otherwise, ranks rows by:
   - non-missing `padj` before missing `padj`;
   - lowest `padj`;
   - largest absolute `log2FoldChange`;
   - source identifier as the final tie-breaker.
6. Keeps the first ranked row as the representative and records the number and
   status of all duplicates in the audit output.

Different gene IDs may represent different transcripts of the same gene, 
and some transcripts may be better biological representatives than others. If all 
duplicates show the same direction of change, we assume that this 
uncertainty is unlikely to have a major effect on the overlap results. 
However, this is still not an ideal solution, and duplicates should 
preferably be resolved before running the pipeline.

The four overlap classes are:

- `UP-UP`: upregulated in cancer and OIS.
- `UP-DOWN`: upregulated in cancer and downregulated in OIS.
- `DOWN-UP`: downregulated in cancer and upregulated in OIS.
- `DOWN-DOWN`: downregulated in cancer and OIS.

The background universe for each cancer is the set of gene symbols tested in
both the cancer and OIS datasets. BH correction in each `overlap_statistics.csv`
is performed across the four tests for that cancer.

### `heatmap_final.R`

Combines the cancer-level `overlap_statistics.csv` files, applies one global BH
correction across all cancer-by-class tests shown in the heatmap, and saves the
final PNG, PDF, and plotted data table. Only globally significant cells are
filled; the color represents `log2(odds ratio)`.

Significance labels:

```text
*     BH < 0.05
**    BH < 0.01
***   BH < 0.001
****  BH < 0.0001
N.S.  BH >= 0.05
```

## Main input: two DEG datasets

The pipeline requires two types of differential-expression results:

1. One OIS differential-expression table (`senescence_OIS_genes.csv`).
2. One or more TCGA cancer differential-expression tables (one file per cancer).

Required columns in the OIS table:

```text
gene, log2FoldChange, padj
```

Required columns in each TCGA cancer table:

```text
gene_id, gene_name, gene_type, log2FoldChange, padj
```

The OIS and TCGA files can be stored in the same input folder. TCGA filenames
must match the pattern configured in `config.R`, for example:

```text
TCGA-BRCA_DGE_Raw_Results_Tumor_vs_Normal_Adjusted.csv
TCGA-LUAD_DGE_Raw_Results_Tumor_vs_Normal_Adjusted.csv
```

## Configuration

The default project-relative paths are:

```r
cancer_input_dir <- file.path(script_dir, "data")
ois_file <- file.path(cancer_input_dir, "senescence_OIS_genes.csv")
results_dir <- file.path(script_dir, "results")

padj_threshold <- 0.05
lfc_threshold <- log2(1.5)
```

Place the input CSV files in `data/`. Personal paths are not stored in the
repository. The `data/` and `results/` folders are excluded by `.gitignore`.

## Required packages

```r
install.packages(c("BiocManager", "ggplot2", "scales"))
BiocManager::install("GeneOverlap")
```

Optional package for higher-quality PNG output:

```r
install.packages("ragg")
```

## Run

Open a terminal in the pipeline folder. Run the overlap analysis first:

```bash
Rscript TCGA_OIS_overlaps.R
```

Check `batch_run_status.csv`. If every cancer has the status `COMPLETED`, run
the heatmap script:

```bash
Rscript heatmap_final.R
```

## Output files

The pipeline creates `results_dir` automatically if it does not exist. It then
creates one subfolder for every TCGA cancer using the project name extracted
from the input filename.

Example output structure:

```text
results_dir/
├── batch_run_status.csv
├── OIS_duplicates_for_review.csv
├── TCGA_OIS_heatmap_cancers_rows_classes_columns.png
├── TCGA_OIS_heatmap_cancers_rows_classes_columns.pdf
├── TCGA_OIS_heatmap_cancers_rows_classes_columns_data.csv
├── TCGA-BRCA/
│   ├── overlap_statistics.csv
│   ├── overlap_counts.csv
│   ├── overlap_all_classes.csv
│   ├── overlap_UP-UP.csv
│   ├── overlap_UP-DOWN.csv
│   ├── overlap_DOWN-UP.csv
│   ├── overlap_DOWN-DOWN.csv
│   ├── common_tested_genes.csv
│   └── CANCER_duplicates_for_review.csv
└── TCGA-LUAD/
    └── ...
```

Files saved directly in `results_dir`:

```text
batch_run_status.csv
OIS_duplicates_for_review.csv
TCGA_OIS_heatmap_cancers_rows_classes_columns.png
TCGA_OIS_heatmap_cancers_rows_classes_columns.pdf
TCGA_OIS_heatmap_cancers_rows_classes_columns_data.csv
```

Each `TCGA-*` subfolder contains:

```text
overlap_statistics.csv
overlap_counts.csv
overlap_all_classes.csv
overlap_UP-UP.csv
overlap_UP-DOWN.csv
overlap_DOWN-UP.csv
overlap_DOWN-DOWN.csv
common_tested_genes.csv
CANCER_duplicates_for_review.csv
```

If the pipeline is run again with the same `results_dir`, files with the same
names are overwritten.

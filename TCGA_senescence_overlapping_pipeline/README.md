# TCGA–senescence overlap analysis

Compare all TCGA cancer datasets with OIS, RS, SIPS, CICQ, and SSCQ signatures.
The pipeline saves shared gene lists and one enrichment heatmap per signature.

## Input files

Set paths and signatures in `config.R`. By default, inputs come from these
repository folders:

- `data/senescence_meta_analysis_deg/`: signature tables and their matching background files.
- `results/tables/TCGA_deg/`: full TCGA differential-expression tables, including non-DEGs.

| File example | Required columns |
| --- | --- |
| `senescence_OIS_genes.csv` | `gene`, `log2FoldChange`, `padj` |
| `senescence_OIS_background.csv` | `gene`, `gene_symbol`, or `gene_name` |
| `TCGA-BRCA_DGE_Raw_Results_Tumor_vs_Normal_Adjusted.csv` | `gene_id`, `gene_name`, `gene_type`, `log2FoldChange`, `padj` |

Each signature needs a matching `_background.csv`, including all eligible genes,
not just DEGs. Quiescence files use the same format, e.g.
`quiescence_SSCQ_genes.csv` and `quiescence_SSCQ_background.csv`.
Use matching gene symbols; every signature gene must be in its background.

## Run

Install packages once in R:

```r
install.packages(c("BiocManager", "ggplot2"))
BiocManager::install("GeneOverlap")
```

From the pipeline folder:

```bash
Rscript TCGA_senescence_overlaps.R
Rscript heatmap_final.R
```

For files in `Downloads/meta_analysis/` and `Downloads/TCGA_deg/`, set these
before running to use the original two-folder input layout:

```bash
export TCGA_SENESCENCE_BASE_DIR="$HOME/Downloads"
export TCGA_SENESCENCE_RESULTS_DIR="$HOME/Downloads/TCGA_senescence_results"
```

To redraw graphs after a successful analysis, run only `heatmap_final.R`.

## Analysis

Genes are classified using `padj < 0.05` and `abs(log2FoldChange) > log2(1.5)`.
Four overlaps are tested: `UP-UP`, `UP-DOWN`, `DOWN-UP`, `DOWN-DOWN`
(cancer first, signature second).

Each one-sided Fisher test uses **TCGA protein-coding genes with non-missing
adjusted p-values and fold changes ∩ signature
background** as its universe. Both DEG sets are restricted to this universe.
Duplicate symbols with opposing significant directions are excluded;
other duplicates use one representative and are exported for review.

Heatmap colors show log2 odds ratios for significant overlaps. Stars use
BH-adjusted p-values across all cancer-by-class tests within each signature.

## Results

The configured results folder (default `results/`) contains gene lists,
counts, statistics, duplicate-review tables, and `common_tested_genes.csv`
with the exact background. Graphs and plotted data are in `heatmaps/`.
Check `batch_run_status.csv` for failed comparisons. Rerunning overwrites
matching output files. Local `data/` and `results/` folders are ignored by Git.

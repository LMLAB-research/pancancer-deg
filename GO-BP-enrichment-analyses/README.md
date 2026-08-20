# GO enrichment of cancer–OIS overlaps

This pipeline performs functional interpretation of genes shared between
cancer and oncogene-induced senescence (OIS) differential-expression results.
It tests GO Biological Process enrichment, groups redundant GO terms by
semantic similarity, and summarizes recurrent biological themes across
cancers.

## Workflow

```text
directional cancer–OIS overlap gene lists
    ↓
GO Biological Process enrichment
    ↓
global semantic clustering of significant GO terms
    ↓
cross-cancer theme ranking
    ↓
TOP5 and TOP10 review tables and figures
```

## Input

The input is the `results/` folder created by the previous overlap pipeline:

```text
TCGA_OIS_overlapping_pipeline/results/
├── TCGA-BLCA/
│   ├── common_tested_genes.csv
│   ├── overlap_UP-UP.csv
│   ├── overlap_UP-DOWN.csv
│   ├── overlap_DOWN-UP.csv
│   └── overlap_DOWN-DOWN.csv
└── TCGA-*/
    └── ...
```

Every input CSV must contain a `gene_symbol` column.

- `common_tested_genes.csv` contains genes tested in both DEG datasets and is
  used as the cancer-specific enrichment background.
- The four `overlap_*.csv` files contain the directional cancer–OIS overlaps.

The upstream overlap pipeline uses adjusted `p < 0.05` and
`|log2FoldChange| > log2(1.5)` for DEG filtering. Overlap enrichment is tested
with `GeneOverlap::testGeneOverlap()`, which uses a one-sided Fisher exact test
from `stats::fisher.test(..., alternative = "greater")`. P-values are corrected
with the Benjamini–Hochberg method.

## Running in RStudio

### 1. Load the launcher

Open `run_pipeline.R` in RStudio and click **Source**.

### 2. Install packages once

Run in the RStudio Console:

```r
run_go_pipeline("install")
```

Packages are installed into the local `.Rlib/` folder.

### 3. Run the complete analysis

If the overlap and GO pipeline folders are next to each other, run:

```r
run_go_pipeline("full")
```

Otherwise, provide the upstream `results/` folder:

```r
run_go_pipeline(
  "full",
  input = "/path/to/TCGA_OIS_overlapping_pipeline/results"
)
```

This creates both TOP5 and TOP10 results.

### 4. Review theme names

Open:

```text
results/TOP5_theme_review.csv
results/TOP10_theme_review.csv
```

Each row contains the automatic representative GO term, all GO terms assigned
to its semantic cluster, enrichment summaries, and two empty columns:

- `manual_theme_name` — reviewed name for the figure;
- `manual_note` — optional explanation.

If `manual_theme_name` is empty, the automatic GO description is used. Manual
names replace only the displayed label; statistics, clustering, and ranking are
preserved.

### 5. Rebuild the figures

After saving the review CSV files, run in the RStudio Console:

```r
run_go_pipeline("plot")
```

## Analysis method

### GO enrichment

`01_GO_enrichment.R` runs `clusterProfiler::enrichGO()` separately for every
cancer and overlap direction.

- Ontology: GO Biological Process.
- Identifier type: gene symbol checked against `org.Hs.eg.db`.
- Background: the corresponding `common_tested_genes.csv`.
- Multiple-testing correction: Benjamini–Hochberg.
- Significance threshold: BH-adjusted `p < 0.05`.
- GO gene-set size: 10–500 genes.

### Global semantic clustering

`02_global_semantic_themes.R` combines the union of significant GO terms from
all enrichment analyses. Pairwise semantic similarity is calculated with
`GOSemSim::mgoSim()` using the Wang method. Redundant terms are grouped with
`rrvgo::reduceSimMatrix()`.

Primary settings:

```text
ontology: GO Biological Process
semantic similarity: Wang
similarity threshold: 0.70
```

Sensitivity outputs are also calculated for Wang thresholds 0.60 and 0.80 and
for the Rel similarity method at 0.70. Cluster agreement is summarized with
the adjusted Rand index.

### Automatic cluster representative

One automatic GO representative is selected for every semantic cluster by:

1. largest number of significant cancer/direction lists;
2. lowest median BH-adjusted p-value;
3. highest median log2 fold enrichment;
4. GO ID as the final tie-break.

The automatic representative provides a reproducible label. Manual review of
all member GO terms is used to assign the final biological name.

### Cross-cancer ranking

Themes are ranked independently within each overlap direction by:

1. number of cancers, decreasing;
2. median log2 fold enrichment, decreasing;
3. median GeneRatio, decreasing;
4. GO ID as the final tie-break.

TOP5 and TOP10 are generated automatically from the same ranking. TOP5 is the
compact main summary. TOP10 provides a broader review but produces a denser
figure.

## Main outputs

```text
results/
├── 00_GO_enrichment_run_log.csv
├── 01_GO_enrichment_all_results.csv
├── 02_GO_significant_terms.csv
├── 03_global_semantic_cluster_mapping.csv
├── 04_theme_per_cancer_class.csv
├── 05_cross_cancer_theme_summary.csv
├── 06_TOP5_cross_cancer_themes_per_class.csv
├── 06_TOP10_cross_cancer_themes_per_class.csv
├── 07_TOP5_cross_cancer_themes_with_manual_names.csv
├── 07_TOP10_cross_cancer_themes_with_manual_names.csv
├── 08_semantic_sensitivity_summary.csv
├── 09_semantic_sensitivity_assignments.csv
├── 10_validation_report.txt
├── TOP5_theme_review.csv
├── TOP10_theme_review.csv
├── per_cancer_GO/
├── per_cancer_plots/
└── figures/
    ├── GO_TOP5_theme_map.png
    ├── GO_TOP5_theme_map_data.csv
    ├── GO_TOP10_theme_map.png
    └── GO_TOP10_theme_map_data.csv
```

Identical final display names are shown once on the y-axis with one point for
each relevant overlap direction.

## Files

- `run_pipeline.R` runs the complete workflow from RStudio.
- `install_packages.R` installs the required packages.
- `config.R` contains analysis parameters and portable paths.
- `pipeline_utils.R` contains shared helper functions.
- `01_GO_enrichment.R` performs GO enrichment.
- `02_global_semantic_themes.R` performs semantic clustering, ranking, and
  sensitivity analysis.
- `03_GO_enrichment_map.R` applies reviewed names and creates both figures.
- `04_validate_results.R` validates the analysis and output tables.

## References

- [GeneOverlap](https://bioconductor.org/packages/GeneOverlap/)
- [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/)
- [Wang et al. (2007)](https://doi.org/10.1093/bioinformatics/btm087)
- [GOSemSim](https://bioconductor.org/packages/GOSemSim/)
- [rrvgo](https://bioconductor.org/packages/rrvgo/)

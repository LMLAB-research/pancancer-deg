# Independent GO BP reanalysis: how to run and audit it

## What this analysis asks

For each TCGA cancer and each directional overlap class, the analysis asks:

> Are genes annotated to a GO Biological Process term represented more often
> in the overlap list than expected from the genes that could have entered that
> cancer-specific overlap?

The four lists are analysed separately:

- UP-UP
- UP-DOWN
- DOWN-UP
- DOWN-DOWN

This is over-representation analysis (ORA), not GSEA.

## Inputs

For every `TCGA-*` directory, the script reads:

- `common_tested_genes.csv`: the cancer-specific statistical background;
- `overlap_UP-UP.csv`;
- `overlap_UP-DOWN.csv`;
- `overlap_DOWN-UP.csv`;
- `overlap_DOWN-DOWN.csv`.

The foreground is forced to be a subset of the corresponding background.

## Statistical test

`clusterProfiler::enrichGO()` performs a GO over-representation test. For each
GO term it compares:

|                         | In GO term | Not in GO term |
|-------------------------|-----------:|---------------:|
| Overlap gene list       | a          | b              |
| Remaining background    | c          | d              |

The over-representation p-value is the upper tail of the hypergeometric
distribution and is equivalent to a one-sided Fisher exact test with
`alternative = "greater"` for this 2 × 2 table.

Benjamini-Hochberg correction is applied within each cancer/class GO analysis.

## Folder structure

The script creates:

```text
GO_BP_independent_reanalysis/
├── GO_BP_all_results.csv
├── GO_BP_significant_results_BH_0.05.csv
├── GO_BP_run_log.csv
├── GO_BP_mapping_QC.csv
├── GO_BP_analysis_settings.csv
├── GO_BP_session_info.txt
├── tables/
│   ├── TCGA-BLCA/
│   │   ├── UP-UP_GO_BP_all_results.csv
│   │   ├── UP-UP_GO_BP_significant_BH_0.05.csv
│   │   └── ...
│   └── ...
├── plots/
│   ├── TCGA-BLCA/
│   │   ├── UP-UP_GO_BP_dotplot.png
│   │   ├── UP-UP_GO_BP_dotplot.pdf
│   │   └── ...
│   └── ...
└── audit/
    ├── TCGA-BLCA/
    │   ├── UP-UP_input_gene_audit.csv
    │   ├── UP-UP_background_audit.csv
    │   ├── UP-UP_GO_BP_annotations.csv
    │   └── ...
    └── ...
```

## How to run it in RStudio

```r
source(
  "/Users/danabiruk/Documents/Codex/2026-08-12/referenced-chatgpt-conversation-this-is-an/outputs/run_GO_BP_independent_reanalysis.R"
)
```

## What each main output means

### `GO_BP_run_log.csv`

One row per cancer/class. Check first that no row has `status = FAILED`.

### `GO_BP_mapping_QC.csv`

Reports:

- input symbols before GO mapping;
- genes inside the cancer-specific background;
- genes with at least one GO BP annotation;
- corresponding mapping rates;
- background size before and after GO annotation filtering.

### `GO_BP_all_results.csv`

Combined table of all returned GO tests across cancers/classes. Important
columns include:

- `ID`, `Description`;
- `GeneRatio`, `BgRatio`;
- numeric versions of both ratios;
- `Count`;
- `expected_gene_count`;
- `FoldEnrichment` and `log2_FoldEnrichment`;
- `pvalue`, `p.adjust`, `qvalue`;
- cancer/class and input/background sizes.

### `GO_BP_significant_results_BH_0.05.csv`

Only rows with BH-adjusted p-value below 0.05. No additional correction across
all cancers is imposed by this file.

## First audit after running

```r
root <- paste0(
  "/Users/danabiruk/Documents/all_TCGA_OIS_overlaps/",
  "GO_BP_independent_reanalysis"
)

run_log <- read.csv(file.path(root, "GO_BP_run_log.csv"))
mapping_qc <- read.csv(file.path(root, "GO_BP_mapping_QC.csv"))
all_go <- read.csv(file.path(root, "GO_BP_all_results.csv"))

table(run_log$status)
summary(mapping_qc$input_GO_BP_mapping_rate)
summary(mapping_qc$background_GO_BP_mapping_rate)
table(all_go$cancer_project, all_go$overlap_class)
```

## Manual check of one GO result

For one row:

```r
x <- all_go[1, ]

# Should reproduce the saved fold enrichment.
x$GeneRatio_numeric / x$BgRatio_numeric

# Expected number of overlap genes under the background model.
x$expected_gene_count

# Observed number.
x$Count
```

If `Count` is much larger than `expected_gene_count`, fold enrichment will be
greater than one. Statistical significance is still determined by the test and
BH correction, not by fold enrichment alone.

## Interpretation boundary

GO ORA demonstrates over-representation of annotations in a selected gene
list. It does not directly measure pathway activity, prove causal mechanism, or
show that every cell in a bulk tumour expresses the program.

# Directional senescence–cancer ORA

R pipeline for directional over-representation analysis using GO, MSigDB and
Reactome, followed by Excel compilation. Input overlaps must already be prepared.

## Install dependencies

Run in R from the repository folder:

```r
dir.create(".Rlib", showWarnings = FALSE)
.libPaths(c(normalizePath(".Rlib"), .libPaths()))
install.packages(c("BiocManager", "openxlsx", "msigdbr"), repos = "https://cloud.r-project.org")
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db", "GO.db", "reactome.db"),
                     ask = FALSE, update = FALSE)
```

## Configure inputs

Edit paths in `config.R`, or override them in an optional ignored `config.local.R`.
Defaults use `data/TCGA_deg`, `data/meta_analysis`, and `data/overlaps`.

- Cancer files: `TCGA-<CANCER>_DGE_Raw_Results_Tumor_vs_Normal_Adjusted.csv`;
  columns `gene_id`, `gene_name`, `gene_type`, `log2FoldChange`, `padj`.
- Signature files: `senescence_<SIGNATURE>_genes.csv`; columns `gene`,
  `log2FoldChange`, `padj`. Signatures: OIS, RS, SIPS, CICQ, SSCQ.
- Each signature needs `senescence_<SIGNATURE>_background.csv`, containing all
  eligible/tested genes in a `gene`, `gene_symbol`, or `gene_name` column.
  Do not use only DEGs as background.
- Prepared overlaps: `data/overlaps/<SIGNATURE>/TCGA-<CANCER>/`, containing
  `common_tested_genes.csv` and `overlap_UP-UP.csv`, `overlap_UP-DOWN.csv`,
  `overlap_DOWN-UP.csv`, `overlap_DOWN-DOWN.csv`, each with `gene_symbol`.

Input filename directions are **cancer–senescence**; workbook directions are
**senescence–cancer**. Thus workbook `UP_DOWN` reads `overlap_DOWN-UP.csv`.

## Run

```sh
Rscript 01_GO.R
Rscript 02_MSigDB.R
Rscript 03_Reactome.R
Rscript 04_Compile.R
```

Run all three analyses before compilation. Missing inputs or failed jobs produce
explicit status records and a nonzero exit. Compilation requires all three RDS
files and checks that their input/configuration fingerprints agree.

Results are saved in `results/`, annotation caches in `cache/`, and workbooks in
`results/excel/`. Keep each workbook's `_full/` CSV folder alongside it: these
files retain content that exceeds Excel limits. After input/config/code changes,
rerun all three analyses. Archive annotation caches with analysis outputs.

## Analysis defaults

Protein-coding cancer genes; strict `padj < 0.05` and
`abs(log2FoldChange) > log2(1.5)`. Conflicting duplicate directions are excluded.
ORA uses the supplied eligible background restricted to each database, pathway
sizes 10–500 and at least 5 eligible input genes. One-sided hypergeometric tests
include zero-hit pathways; BH correction is applied separately per database and
job/direction, with `p_adjust < 0.05`. GO defaults to BP. MSigDB includes H,
C2:CGP and C2:CP excluding Reactome; Reactome is analyzed separately.

Private data, local configuration, dependencies, caches and results are excluded
from Git. Real input data are not bundled; a full biological run is not included.

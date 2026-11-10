# Inputs for the exploratory analyses

The files in this folder were copied from Roberto's machine on 2026-11-09 so the exploratory analyses can be
checked and rerun without access to it. `MANIFEST.tsv` lists every file with its original path, size and md5. Every
copy decompresses to exactly the original, except the two SI-table CSVs, which are extracted sheets. Large tables are
gzipped; R `read.csv()` and pandas `read_csv()` read `.csv.gz` directly.

## How the scripts find their inputs

Run every script from the repository root (`Rscript --vanilla exploratory/...`, `python exploratory/...`). The helpers
`exploratory/paths.R` and `exploratory/paths.py` resolve the paths:
- Inputs from other projects (senescence_stress tables, house code, purity, proliferation list) are always read from
  this folder.
- Regenerable intermediates (TCGA DEG tables, within-study contrast DEGs, meta-analysis signature files) are read from
  `results/` or `data/` when a local run produced them, otherwise from this folder.
- Multi-GB data not in the repo (house single-cell objects, raw GEO data) are read from `SENESCENCE_STRESS_DIR`
  (environment variable) and `data/`.

**Run from the repo alone, in this order:**
1. `bash exploratory/run_both.sh` (two-sided overlap pipeline; full and proliferation-free signatures)
2. `STRICT=1 python exploratory/ras_vs_ois/04_classify_and_project.py` (OIS gene classes; `NO_BRAF=1` for the
   sensitivity)
3. `Rscript --vanilla exploratory/metabolism/00_metabolic_sets.R`
4. Any of:
   - `ras_vs_ois/15_pathway_class_overlap.R`
   - `metabolism/01_metabolic_overlap_tcga.R` (also `NOPROLIF=1`)
   - `hallmarks/01_hallmark_overlap_tcga.R`
   - `stressors/01_inducer_tcga_overlap.R`
   - `breakdown/01_cq_rs_sips_components.R`, then `breakdown/02_shared_layer_drivers.R`
   - Python: `ras_vs_ois/03_sanity_checks.py`, `ras_vs_ois/05_class_examples.py`, `ras_vs_myc/06_nis_class_sharing.py`,
     `ras_vs_myc/09_oncogene_class_sharing.py`, `ois_partition.py`, `summarise.py`, `explore.R`

**Need the large data below:**
- `ras_vs_ois/01, 02` (recount3)
- `ras_vs_myc/01, 02, 07, 08, 10` (TCGA counts, cBioPortal, GEO), and the models that use the tumour scores (`03`,
  `04`, `11`, `metabolism/02`)
- the single-cell scripts (`ras_vs_ois/06–09, 11, 13, 14`, `stats_sc_tests.R`, `figures_sc_umaps.R`)
- `house_copy/`

R packages: the exploratory scripts use the senescence_stress package libraries on Roberto's machine when present
(`.libPaths()`), otherwise the user's own library. They need GeneOverlap, clusterProfiler and the other packages loaded
by `inputs/house_code/bulk_functions.R`, plus msigdbr, dplyr, tidyr, ggplot2, patchwork; the single-cell scripts also
need Seurat and UCell. The overlap pipeline uses the repo's renv lockfile.

## Included (91 files, about 95 MB)

| Folder | Files | What | Source |
|---|---|---|---|
| `senescence_stress/` | 17 | Meta-analysis arrest DEGs (`arrest_degs_final_RERUN`, published `arrest_degs_final`); 21 inducer arms (`degs_nongenotoxic`); GSE60340 p53-mutant contrasts; IR time-course DEGs per cell type; house stress modules v3 (acute, chronic, tested genes); house Hallmark + lysosomal sets (`stress_response_pathways_RERUN`); house Hallmark and metabolic overlap results; SI Tables 18 and 37 (metabolic background and pathways); GSE245807 RIS state labels | senescence_stress project (bulk `rerun_outputs/`, `external_inputs/`, `work/`, `Final/final/si_tables.xlsx`, `cs_dataset_survey/notch_secondary/`) |
| `house_code/` | 7 | House R code the scripts source: bulk `functions.R` (live and HANDOFF_2026-08-31 copies), `config.R`, `17_stress_module_overlap_v3.R` (for `add_corr`), scutils `ucell_parallel.R`, single-cell `lib/scRNA_seq.R` and `lib/general_scripts.R` | senescence_stress, scutils |
| `tcga_degs/` | 17 | TCGA tumour-vs-normal DEGs, covariate-adjusted design (all genes) | this repo's pipeline, scripts 01–04 |
| `contrast_degs/` | 37 | Within-study DEG contrasts (RAS/BRAF bypass and escape, knockdowns, NOTCH1, AKT, MYC, cyclin E, hydroxyurea, quiescence) | `ras_vs_ois/02_contrast_degs.R`, `ras_vs_myc/08_akt_myc_degs.R`, `ras_vs_myc/10_akt_ras_array.R` |
| `senescence_meta_analysis_deg/` | 11 | Meta-analysis signatures split per arrest type (genes and backgrounds), the inputs of the two-sided overlap pipeline | built 2026-10-31 from the published `arrest_degs_final.csv` (see its PROVENANCE.txt) |
| `other/` | 2 | Hallmark proliferation gene list (671 genes); TCGA PanCanAtlas ABSOLUTE purity | `two_sided` run; GDC file 4f277128 |

The OIS gene classes are already tracked as `summary_tables/02_ois_gene_classes_strict.csv`.

## Not included (too large for git)

| Input | Size | How to get it |
|---|---|---|
| TCGA counts (`data/TCGA_deg`, 17 cohorts) | ~45 GB | pipeline scripts 01–04 (GDC) |
| TCGA PAAD / SKCM counts (`data/exploratory/tcga_extra`) | ~3 GB | `ras_vs_myc/07_download_paad_skcm.R` |
| recount3 gene counts for the OIS studies (`data/exploratory/ras_vs_ois`) | 62 MB | `ras_vs_ois/01_recount3_metadata.R` (SRP065206, SRP066947, SRP172671, SRP062872, SRP127037, SRP062871, ...) |
| GEO bulk files (`data/exploratory/geo_oncogenes`) | 48 MB | GEO: GSE130099, GSE130100, GSE185512, GSE45276 (+ GPL570 annotation) |
| Single-cell raw data | ~0.5 GB | GEO: GSE222951, GSE202126, GSE115301 |
| House single-cell objects (GSE245807, GSE173879 house reruns) | ~8.4 GB | senescence_stress house pipeline (`cs_dataset_survey/house_rerun`, `house_copy/` wrapper); on Roberto's machine only |
| cBioPortal driver calls | small | `ras_vs_myc/01_driver_status.py` (downloads) |
| msigdbr gene sets | — | R package msigdbr 26.1 |

Not copied: `senescence_stress/work/results/auc/temporal/all_time_degs.csv`. It is corrupt on disk (33 MB, almost all
NUL bytes); the per-cell-type files here replace it.

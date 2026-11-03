# Figure legends

All figures are exploratory; they are in this folder. Gene classes are the strict classes from
`ras_vs_ois/04_classify_and_project.py` (`summary_tables/02_ois_gene_classes_strict.csv`):

- **RAS-autonomous**: the OIS change is kept in RAS cells that have lost the arrest.
- **Arrest-coupled**: the change is lost there.
- **Intermediate**: neither.
- **NF-κB/BRD4-dependent**: the RAS-senescence change falls to ≤ 0.5 in ≥ 2 of shp65, shBrd4 and JQ1.

UP and DOWN give the direction of the change in OIS vs proliferating cells. "MAPK output" =
DUSP6, SPRY2, SPRY4, ETV4, ETV5, FOSL1, EGR1. Choices behind each figure: `../DECISIONS.md`.

---

**Fig. 1 | Overlap of OIS gene classes with TCGA tumour-vs-normal DEGs.**
`fig1_tcga_overlap_by_class.png`

- **Left:** OIS-up genes of each class against genes up in tumours.
- **Right:** OIS-down genes against genes down in tumours.
- **Points:** one per cancer (17 TCGA cohorts, covariate-adjusted DESeq2 design; DEG = padj < 0.05,
  |log2FC| > log2(1.5), protein-coding).
- **y axis:** log2 odds ratio of a two-sided Fisher test. The universe is tumour genes ∩ genes
  tested in the OIS meta-analysis.
- **Markers:** filled = p < 0.05 (not corrected); bar = median across cancers.
- **Brackets:** paired Wilcoxon signed-rank tests between classes across the 17 cancers (paired by cancer), BH over
  the 6 tests (`summary_tables/07_fig1_paired_class_tests.csv`).
- **Source:** `summary_tables/02_tcga_overlap_by_class_strict.csv`. Script: `figures_summary.R`.

---

**Fig. 2 | Fraction of each OIS gene class reproduced by other oncogene arms.**
`fig2_oncogene_arms.png`

- **Rows:** contrasts of an oncogene arm vs its matched control.
  - IMR90 ER:RAS and NOTCH1-ICD (Hoare 2016, SRP062871; n = 6 per arm, technical runs summed).
  - IMR90 RAS and AKT arrays (GSE45276; n = 4; limma).
  - BJ-hTERT HRAS day 14 and myrAKT1 day 6 (GSE130100 / GSE130099; n = 3; separate series).
  - RPE-1 dox-MYC 120 h, p53 wild-type and −/− (GSE185512; n = 2).
  - The MKI67 log2FC of each contrast is shown in the row label.
- **Columns:** gene classes.
- **Cells:** fraction of the class's genes changed in the same direction and passing the house DEG
  thresholds in that contrast (padj < 0.05, |log2FC| > log2(1.5); house `find_degs_new()`, limma
  for arrays).
- **Source:** `summary_tables/04_oncogene_arms_class_sharing.csv`. Script: `figures_summary.R`.

---

**Fig. 3 | Single-cell, time-course and dose tests of the gene classes.**
`fig3_single_cell_time_dose_tests.png`

Scores are per-cell UCell scores (`scutils::score_ucell_parallel`, RNA counts). "Arrest-coupled"
is the arrest-coupled DOWN set; "RAS-autonomous" and "NF-κB/BRD4-dependent" are UP sets.

- **a. ER:RAS IMR90 time course** (GSE115301, Smart-seq2; the authors' published cells; no house
  state filters).
  - Cells: Growing 100 (two batches), day 2 41, day 4 42, senescent 41 (three batches).
  - Median score per time point, scaled 0–1 within each score.
- **b. RAS dose** (GSE222951, RPE1-hTERT ER:HRAS-G12V, 6 days 4-OHT, sorted by RAS-reporter level).
  - Cells after HTODemux singlets and house default QC: S 1,838, M 1,532, L 1,655, XL 1,148.
  - Median per level, scaled 0–1 within each score.
- **c. KRAS time course** (GSE202126, HPNE parental vs KRAS days 1–5, VASA-seq).
  - Run without the 120,000-UMI cap: parental 150, day 1 171, day 2 170, day 3 171, day 5 130
    cells. The capped run keeps 8–40 cells per group and gives the same direction for every class.
  - Median per group, scaled 0–1 within each score.
- **d. RAS withdrawn after OIS** (GSE173879 clone 1, house one-plate pipeline, house-final cells).
  - OIS with 4-OHT kept (260 cells), OIS with 4-OHT replaced by DMSO before harvest (241 cells),
    proliferating control CTRL_1 (2,173 cells).
  - Bar = (median withdrawn − median control) / (median kept − median control).
  - 1 = no change on withdrawal; 0 = back to the control level.
- **e. RAS cells still cycling** (GSE245807; house published-QC, integrity and control-like filters;
  house LMNB1/G1 subset not applied).
  - Bar = (median RAS-cycling − median growing-cycling) / (median RAS-G1 house-final − median
    growing-G1).
  - Cells: RAS cycling 621, RAS G1 house-final 1,425, growing cycling 3,052, growing G1 2,677.
- **Markers in a–c:** square = reference group; filled = BH < 0.05 vs the reference; open = not significant. Tests:
  house pairwise Wilcoxon on per-cell scores (`stats_sc_tests.R`; full results in the `stats_*` figures below and
  `summary_tables/06_stats_*_wilcox.csv`).
- **Error bars in d–e:** 95% bootstrap CI, 2,000 resamples of cells within each group, seeded per score
  (`summary_tables/06_bootstrap_ratios.csv`).
- **Sources:** `summary_tables/03_*`, `06_*`. Script: `figures_summary.R`.

---

**Fig. 4 | What each OIS gene class tracks across TCGA tumours.**
`fig4_tcga_activity_models.png`

- **Model:** within each of 20 cohorts (primary tumours; SKCM also metastases), a class score
  regressed on ERK output, Hallmark MYC targets (v1 + v2), proliferation (Hallmark E2F + G2M),
  a stromal proxy (10 genes) and ABSOLUTE purity (+ sample type for SKCM). Ordinary least
  squares; all variables standardised within cohort.
- **Scores:** mean within-cohort z of VST expression of the set's genes.
- **Points:** one per cohort and predictor (standardised coefficient). Filled = p < 0.05 (not
  corrected); bar = median.
- **Stars:** Wilcoxon signed-rank test across the 20 cohorts that the coefficient differs from 0, BH over 12 tests.
  **Bracket:** paired Wilcoxon across cohorts, ERK vs MYC coefficient, BH over 4 tests
  (`summary_tables/07_fig4_*.csv`). Stars: * < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001.
- **Caveat:** predictors are correlated, so each coefficient is conditional on the others.
- **Source:** `summary_tables/05_tcga_activity_models_stroma_purity.csv`. Script: `figures_summary.R`.

---

**Fig. 5 | TCGA tumours grouped by driver.**
`fig5_tcga_driver_models.png`

- **Groups** (cBioPortal TCGA PanCanAtlas):
  - RAS-hotspot = KRAS/NRAS/HRAS missense at codons 12/13/61, or BRAF V600.
  - MYC-amplified = GISTIC = 2.
  - Each group is compared with tumours carrying neither; tumours with both are excluded.
  - Cohorts with ≥ 10 tumours in both groups.
- **Effect:** coefficient of the group indicator for a score standardised within cohort, adjusted
  for the stromal proxy and ABSOLUTE purity (+ sample type for SKCM).
- **Colour:** blue = RAS-hotspot vs neither; orange = MYC-amplified vs neither (each cohort row has one point per
  driver that met the group-size rule).
- **Marker:** filled = BH < 0.05 (correction across all cohort × score × driver tests); open = p < 0.05 before
  correction only; cross = p ≥ 0.05.
- **Source:** `summary_tables/05_tcga_driver_models_stroma_purity.csv`. Script: `figures_summary.R`.

---

**UMAP figures.** Point size 3 for < 1,000 cells, 1.6 for 1,000–3,000, 1 above. (`figures_sc_umaps.R`; house plotting calls: `DimPlot` clusters first, no cluster
labels; `customSignaturePlot` for scores, each rescaled 0–10 within the plot, so colours are not
comparable across panels or files).

- **`umap_GSE245807_RAS_final_states.png`**: GSE245807 (Olan 2024) IMR90 ER:RAS, house-final RAS
  cells (1,425; plates RIS_1 and RIS_2).
  - Embedding and clusters from the house run; cluster 3 = state 3 (MAPK output at growing level).
  - Panels: clusters, plate, then the class scores of script 07.
- **`umap_GSE245807_RAS_cycling_vs_G1.png`**: the same RAS cells before the house LMNB1/G1 subset
  (2,062 cells after published QC, integrity and control-like filters).
  - Groups: cycling (S/G2M, 621), G1 house-final (1,425), G1 LMNB1-positive (16).
  - Embedding recomputed here: SCTransform, PCA 50, UMAP dims 1:30, clusters at resolution 0.12,
    seed 1. It is not a house-run embedding, and cluster numbers start at 0.
  - Scores: script 09.
- **`umap_GSE173879_OISON.png` / `umap_GSE173879_OISOFF.png`**: GSE173879 (Evans 2023) LF1
  ER:HRAS clone 1, house-final OIS cells, 4-OHT kept (260) / withdrawn (241).
  - One plate per arm, so no integration and every cluster kept (house one-plate rule).
  - Scores: script 13.

---

**Statistics figures `stats_*.png`** (`stats_sc_tests.R`). House `wilcox_dotplot_wrapper()` from
`senescence_stress/.../lib/scRNA_seq.R`:

- **Plot:** one panel per score; violin + jittered cells.
- **Panel title:** Kruskal–Wallis significance, BH over the 6 scores.
- **Brackets:** all pairwise Wilcoxon rank-sum tests, BH over all pairs × scores in the figure; only significant
  brackets drawn.
- **Tables:** `summary_tables/06_stats_<tag>_{kw,wilcox}.csv`.
- **Unit:** cells are the test unit. Most datasets have 1–2 biological replicates per group, so the p-values
  describe cell-to-cell differences within these cultures, not replicate experiments.

| File | Groups |
|---|---|
| `stats_GSE115301_time_course.png` | Growing, day 2, day 4, senescent |
| `stats_GSE222951_RAS_dose.png` | S, M, L, XL |
| `stats_GSE202126_KRAS_days.png` | parental, days 1, 2, 3, 5 (no UMI cap) |
| `stats_GSE173879_withdrawal.png` | proliferating, OIS 4-OHT kept, OIS 4-OHT withdrawn |
| `stats_GSE245807_cycling.png` | growing G1, growing cycling, RAS cycling, RAS G1 house-final |
| `stats_clusters_GSE245807_final_states.png` | house final states 1–3 |
| `stats_clusters_GSE245807_prefilter.png` | clusters 0, 1 of the recomputed pre-filter embedding |
| `stats_clusters_GSE173879_OISON.png`, `..._OISOFF.png` | house clusters 1, 2 of each arm |

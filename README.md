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
|   |-- TCGA_senescence_overlap/
|   `-- TCGA_senescence_enrichment/
|-- data/
|-- docs/
|-- results/
|   |-- figures/
|   |   |-- TCGA_deg/
|   |   |-- TCGA_senescence_overlap/
|   |   `-- TCGA_senescence_enrichment/
|   `-- tables/
|       |-- TCGA_deg/
|       |-- TCGA_senescence_overlap/
|       `-- TCGA_senescence_enrichment/
|-- setup/
|   |-- check_setup.R
|   |-- install_dependencies.R
|   `-- setup.sh
|-- README.md
|-- renv.lock
`-- pancancer_deg.Rproj
```

The `analyses/` directory contains analysis scripts grouped by analysis module. Current modules are `TCGA_deg/`, `TCGA_senescence_overlap/` and `TCGA_senescence_enrichment/`; future modules can be added as sibling directories under `analyses/`.

The `results/` directory is organized by output type first, then by analysis module. For example, the TCGA differential expression analysis writes figures to `results/figures/TCGA_deg/` and tables to `results/tables/TCGA_deg/`.

The `docs/` directory holds dated notes on findings and follow-up ideas (index in `docs/README.md`), so that this README stays focused on how the analyses work.

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

Steps: `deg01`–`deg08` (TCGA_deg), `overlap01`–`overlap02`
(TCGA_senescence_overlap) and `enrich01`–`enrich05` (TCGA_senescence_enrichment);
groups `deg`, `overlap`, `enrich`, `all`. Run
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
Rscript analyses/TCGA_deg/07_compare_published_degs.R
Rscript analyses/TCGA_deg/08_pca_diagnostics.R
```

| Step | What it does |
| --- | --- |
| 01 | Lists TCGA projects with at least 5 normal and 5 tumor RNA-seq samples. |
| 02 | Downloads STAR counts and clinical data. Tumor = primary tumor, normal = solid tissue normal; other sample types are excluded. Sums replicate files of the same patient and condition (below). |
| 03 | Builds every DESeq2 design (below) and writes the design plan for review. |
| 04 | Fits every eligible design (the slow step) and saves the fitted DESeq2 objects. |
| 05 | Tests each fitted design (see *Differential expression test*), shrinks log2 fold changes (ashr), and writes raw and shrunken result tables. |
| 06 | Compares naive, adjusted and paired results per cancer and sex model (fold-change correlation, DEG agreement), and summarises the sex models per cancer and version. |
| 07 | Compares our results with the published cancer-DEGs of Chatsirisupachai et al. 2019 (see *Comparison with published DEGs*). |
| 08 | PCA diagnostics: which sample variables drive the main axes of variation (see *What drives the main axes of variation*). |

Each step saves its outputs, so later changes need only a partial rerun:

| To change | Rerun from |
| --- | --- |
| test, FDR level or fold-change threshold | `deg05` (minutes; no refit) |
| designs or covariate rules | `deg03` |
| sample selection or clinical merging | `deg02` (downloads are reused) |

Step 04 can be interrupted and resumed: rerunning `from:deg04` reuses every
design fitted since step 03 and fits only the others (set
`reuse_existing_fits <- FALSE` in `00_config.R` to refit everything). Rerunning
step 03 always forces a refit.

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

### Comparison with published DEGs

Step 07 compares our results with the cancer-DEGs of Chatsirisupachai,
Palmer, Ferreira & de Magalhães (2019), *Aging Cell* 18:e13041 (PMC6826163),
Data S2: 10 TCGA projects (BRCA, COAD, ESCA, GBM, LIHC, LUAD, PRAD, STAD, THCA,
UCEC), RSEM counts from FireBrowse (GRCh37), limma-voom, primary tumor vs
solid normal without covariates, DEG = BH-adjusted p < 0.01 and
$|\log_2\text{FC}| > 1$. Data S2 is downloaded automatically from Europe PMC
into `data/published_degs/` (checksum-verified). Genes are matched by Ensembl
ID; Data S2 lists only their DEGs, so every measure is over their DEGs:

- **direction concordance:** share of their DEGs with the same sign of
  log2FC in our results (threshold-free);
- **Spearman correlation** of their and our log2FC;
- **recovery:** share of their up (down) DEGs that are up (down) DEGs in our
  results, with our default calls and with their rule applied to our fits;
- **precision and Jaccard index** (their rule on both sides, protein-coding
  genes).

The comparison uses our naive design (theirs) and the paired design. The two
pipelines also differ in genome build and annotation (GRCh37 RSEM vs GRCh38
STAR), method (limma-voom vs DESeq2), gene filters and sample sets, so perfect
agreement is not expected.

### What drives the main axes of variation

Step 08 runs a PCA per cancer (variance-stabilized counts of the naive fit,
500 most variable genes, PCs 1–10) and relates each PC to sample variables by
the adjusted $R^2$ of PC ~ variable: condition, sex, race, age, library size,
plate and tissue source site (from the TCGA barcode), and the annotations of
the TCGA marker papers (`paper_*` columns: mostly molecular subtypes and
histologies; among tumors only). It repeats the PCA on paired samples only
(one normal and one tumor per patient). Results (full run, 18 cancers):

- **Tumor vs normal is PC1** in 9 cancers: KICH (adjusted $R^2$ 0.90), CHOL
  (0.75), KIRC (0.74), KIRP (0.74), LUSC (0.53), LUAD (0.47), COAD (0.46), and
  more weakly LIHC and READ (0.32).
- **Elsewhere PC1 is tumor heterogeneity**, and tumor vs normal appears only
  on PC2–PC4 ($R^2$ 0.12–0.30). Among tumors, PC1 matches published subtypes
  or histologies: ESCA histology/mRNA cluster (0.86), BLCA mRNA cluster
  (0.78), THCA mRNA cluster (0.71), BRCA mRNA clusters (0.69), UCEC mRNA
  cluster (0.51), PRAD iCluster (0.48); weaker in HNSC (0.31) and STAD (0.26).
- **Why:** PCA ranks directions by total variance. A tumor/normal difference
  $\delta$ adds variance $p(1-p)\,\delta^2$, where $p$ is the fraction of
  normals; with $p \approx 0.05$–$0.10$ this is small compared with the
  differences between tumor subtypes. Differential expression tests are not
  affected. In the paired-only PCA ($p = 0.5$), tumor vs normal becomes PC1
  in most cancers ($R^2$ 0.76–0.93, including BRCA 0.81 and UCEC 0.76), but
  stays weak in THCA (0.47), PRAD (0.24), BLCA (0.20), HNSC (0.19), STAD and
  ESCA (≈ 0), where tumor heterogeneity remains larger than the tumor/normal
  difference.
- **GBM:** PC1 is the sequencing plate (0.77), a batch effect.
- **Clinical and technical variables** explain little of PC1–PC2 elsewhere:
  sex ≤ 0.08, age ≤ 0.11, race ≤ 0.27 (ESCA, through histology), library size
  ≤ 0.36.

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
- `07_TCGA_published_deg_comparison.csv` and
  `results/figures/TCGA_deg/07_TCGA_published_deg_comparison.png`: agreement
  with the published cancer-DEGs (step 07).
- `08_TCGA_pca_summary.csv` (one row per cancer), `08_TCGA_pca_associations.csv`
  (every PC x variable) and `08_TCGA_pca_variance.csv`; figures
  `08_<project>_PCA_condition.png` and `08_<project>_PCA_paper_annotation.png`
  in `results/figures/TCGA_deg/` (step 08).
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
- **Tumor effects average over subtypes.** In heterogeneous cancers, PC1 is
  tumor subtype or histology (step 08), and $\beta$ is the average tumor effect
  over the subtype mix of each cohort; a gene can change in one subtype and
  not in another.
- **ESCA mixes two cancers.** Squamous carcinoma and adenocarcinoma arise from
  different epithelia; PC1 is histology (adjusted $R^2$ 0.86) and even the
  paired PCA shows no tumor/normal separation. Its DEGs average over both
  (STAD shows a similar, weaker pattern).
- **GBM.** Its 5 normals have no clinical data and no matched tumor, so GBM
  has only a naive design, based on 5 normals of unknown background, and its PC1
  is the sequencing plate (a batch effect).

### To consider

- **Outlier samples.** In BLCA, 5 of 19 normals form a separate cluster in the
  paired PCA (step 08). Check BLCA, then flag samples far from their group in
  every cancer.
- **ESCA by histology.** Analyse squamous carcinoma and adenocarcinoma
  separately.
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

Four readouts are computed. The first three use Fisher's exact test on a 2x2
table:

1. **Class overlaps** (reference). For each class $T^a \cap S^b$,
   $a, b \in \{+, -\}$ (UP-UP, UP-DOWN, DOWN-UP, DOWN-DOWN; cancer first), the
   table is "in $T^a$ or not" x "in $S^b$ or not" over $U$. Two-sided by default
   (`fisher_alternative`). BH across all cancer x class tests within each
   signature.
2. **Sharing.** Over $U$: "DE in cancer or not" x "DE in signature or not".
   Odds ratio > 1: DE genes are shared more than by chance.
3. **Direction.** Over genes DE in both: "up or down in cancer" x "up or down
   in signature". Odds ratio > 1: same direction; < 1: opposite direction.

4. **Fold-change correlation** (threshold-free). Spearman correlation $\rho$
   between tumor and signature log2 fold changes over all genes of $U$.
   $\rho > 0$: same overall direction; $\rho < 0$: opposite. It checks the
   direction readout without any DEG cutoff.

Readouts 2–4 are always two-sided, with BH across cancers within each
signature.

How they relate: the class-overlap enrichment mixes sharing and direction. For
UP-UP, with ratios of fractions,

$$\frac{P(S^+ \mid T^+)}{P(S^+ \mid \text{not } T^+)}
= \frac{P(S\ \text{DE} \mid T^+)}{P(S\ \text{DE} \mid \text{not } T^+)}
\times \frac{P(S^+ \mid T^+, S\ \text{DE})}{P(S^+ \mid \text{not } T^+, S\ \text{DE})}.$$

Because sharing is strong, all four classes can be enriched at once even when
the shared genes move in opposite directions. Use the direction readout (and
the fold-change correlation) to compare directions.

**Reading the p-values.** All four tests treat genes as independent, but genes
are co-regulated in modules, so the p-values are too optimistic; with about
13,000 genes, even $\rho = 0.03$ is "significant". The effect sizes (odds
ratios with their confidence intervals, $\rho$) are the informative
quantities.

**DEG definitions differ between the two sides.** Tumor DEGs come from TCGA_deg
step 05 (by default the threshold test $H_0: |\beta| \le \log_2 1.5$).
Signature DEGs are genes with padj < 0.05 and $|\log_2\text{FC}| > \log_2 1.5$
in the meta-analysis tables (identical to their `significant` column); these
tables give no standard errors, so the threshold test cannot be applied to
them.

Genes are matched by symbol. Signature symbols that were renamed since the
signatures were built (e.g. `MRVI1` → `IRAG1`) are first updated: a symbol
absent from TCGA and no longer an official symbol is replaced by its current
symbol when the official aliases (`org.Hs.eg.db`) give exactly one candidate
present in TCGA. This recovers about 30 genes per signature; the renamings are
listed in `<signature>/<signature>_symbol_mapping.csv`.

Duplicate gene symbols with opposite significant directions are excluded; other
duplicates keep one representative and are listed for review.

### Outputs

- `results/tables/TCGA_senescence_overlap/<signature>/<project>/`: shared gene
  lists, class-overlap statistics (`overlap_statistics.csv`), sharing and
  direction statistics (`share_direction_statistics.csv`), fold-change
  correlation (`lfc_correlation_statistics.csv`) and the exact
  universe (`common_tested_genes.csv`). The statistics tables also give the
  cancer's sample sizes (`n_normal`, `n_tumor`, from TCGA_deg step 05) and the
  numbers of cancer and signature DEGs in the universe, so that results can be
  interpreted against power (cancers with more samples have more DEGs).
- `results/tables/TCGA_senescence_overlap/batch_run_status.csv`: one row per
  cancer x signature; check it for failures.
- `results/figures/TCGA_senescence_overlap/TCGA_<signature>_heatmap.png`: log2
  odds ratios of significant class overlaps (blue = enrichment, red = depletion).
- `results/figures/TCGA_senescence_overlap/TCGA_share_heatmap.png`,
  `TCGA_direction_heatmap.png` and `TCGA_lfc_correlation_heatmap.png`: cancers
  x signatures (direction and correlation: blue = same, red = opposite).
- The plotted values of each heatmap are in
  `results/tables/TCGA_senescence_overlap/TCGA_*_heatmap_data.csv`.
- `results/tables/TCGA_senescence_overlap/<signature>/<signature>_signature_genes.csv`
  and `_signature_background.csv`, and `cancer_gene_tables/<project>_genes.csv`:
  prepared signature and tumor gene tables (after symbol mapping and duplicate
  handling), used by the enrichment module.

## TCGA_senescence_enrichment Analysis

Which pathways change in tumors and in senescence or quiescence, and do they
change in the same or opposite direction? Settings are in
`analyses/TCGA_senescence_enrichment/00_config.R`; the TCGA design (paired by
default), signatures and input folder are read from the overlap module's
config. Run from the repository root:

```bash
Rscript analyses/TCGA_senescence_enrichment/01_gene_sets.R
Rscript analyses/TCGA_senescence_enrichment/02_ora.R
Rscript analyses/TCGA_senescence_enrichment/03_gsea.R
Rscript analyses/TCGA_senescence_enrichment/04_compare.R
Rscript analyses/TCGA_senescence_enrichment/05_workbooks.R
```

| Step | What it does |
| --- | --- |
| 01 | Downloads the gene sets from MSigDB (`msigdbr`) and saves them with their version. |
| 02 | Over-representation analysis (ORA) of gene lists. |
| 03 | Preranked gene-set enrichment analysis (GSEA, `fgsea`) of tumor and signature rankings. |
| 04 | Compares tumor and signature GSEA results per gene set. |
| 05 | Excel workbooks of the main results (optional). |

### Gene sets

From MSigDB (`msigdbr` 26.1.1): HALLMARK (H, 50 sets), GO_BP (C5 GO:BP, 7,538),
REACTOME (C2 CP:REACTOME, 1,839), KEGG (C2 CP:KEGG_MEDICUS, 658) and CGP (C2
CGP, 3,555: published perturbation signatures). Gene sets are tested if they
have 10–500 genes in the analysis universe.

### Methods

The module uses only the outputs of the earlier modules (tumor DEG calls,
shrunken log2 fold changes, overlap classes and prepared gene tables), so all
analyses use the same genes. Three complementary approaches:

- **B. Each side separately, then compared (main analysis).** GSEA on the
  ranking of all tumor genes (shrunken log2FC) and of all signature genes
  (log2FC), without DEG cutoffs. NES > 0: the gene set is concentrated among
  up-regulated genes; NES < 0: among down-regulated genes. Step 04 compares
  $\text{NES}_\text{cancer}$ and $\text{NES}_\text{signature}$ per gene set
  (same or opposite direction when significant on both sides) and correlates
  them across gene sets (pathway-level concordance).
- **C. Joint test within gene sets.** For gene sets significant on both sides:
  the Jaccard index of the two leading edges (do the same genes drive both
  enrichments?) and the Spearman correlation $\rho_\text{set}$ of tumor and
  signature log2FC over the genes of the set, next to the genome-wide
  $\rho_\text{genome}$.
- **A. Intersection (descriptive).** ORA (one-sided hypergeometric test) of
  the four overlap classes (UP-UP, UP-DOWN, DOWN-UP, DOWN-DOWN; cancer first;
  universe = genes tested in both), next to ORA of tumor and signature DEGs
  alone. A gene set enriched on each side separately is also enriched in their
  intersection, so class results describe the shared genes but are not
  evidence of a joint effect.

BH correction per gene list or ranking and collection; significant if
BH < 0.05.

GSEA is **unweighted** (`gsea_weight <- 0`, fgsea `gseaParam = 0`): only gene
ranks count. In some tumors the largest changes are almost all losses (95–100%
of the top 100 |log2FC| in KICH and KIRP), which makes the null distribution of
weighted GSEA one-sided, so positively enriched gene sets get no valid p-value.
Unweighted GSEA avoids this and makes tumor and signature NES comparable
(HALLMARK: Spearman 0.92 with weighted NES). Gene sets that the default method
still cannot score are rerun by permutation (10,000; none in the current run).

### First results (paired design, threshold test)

- **Pathways move in opposite directions in tumors and in quiescence and
  senescence:** median correlation of tumor and signature NES across gene sets
  −0.40 to −0.63 for CICQ, SSCQ and SIPS, −0.15 to −0.43 for RS, depending on
  the collection; strongest in LUAD, LUSC and UCEC (down to −0.87).
- **OIS is the exception:** NES correlation about 0 (−0.08 to +0.12),
  slightly positive in kidney cancers, THCA and ESCA, in line with the
  gene-level direction readout.
- **Proliferation is the core of the opposition:** HALLMARK E2F targets and
  G2M checkpoint are opposite in all 85 cancer x signature pairs, mitotic
  spindle in 80, MYC targets in 68.
- **Gene-level support:** among gene sets significant on both sides, the
  leading edges overlap moderately (median Jaccard 0.35); 63% of same-direction
  gene sets have $\rho_\text{set} > \rho_\text{genome}$, but only 50% of
  opposite-direction ones have $\rho_\text{set} < \rho_\text{genome}$: many
  pathways are opposite as a whole without their genes being more opposite
  than the genome average.

### Outputs

In `results/tables/TCGA_senescence_enrichment/`:

- `01_gene_set_versions.csv`; gene sets in `data/TCGA_senescence_enrichment/gene_sets/`.
- `02_ORA_<collection>.csv` (rows with p < 0.05, with the matched genes) and
  `02_ORA_status.csv` (input and universe sizes, gene sets tested per list).
- `03_GSEA_<collection>.csv` (NES, p, padj per ranking and gene set; leading
  edges for padj < 0.05) and `03_GSEA_status.csv`.
- `04_compare_<collection>.csv` (gene sets significant on at least one side,
  with direction, leading-edge Jaccard, $\rho_\text{set}$, $\rho_\text{genome}$)
  and `04_compare_summary.csv` (one row per cancer x signature x collection).
- `workbooks/<project>_enrichment.xlsx` and `workbooks/signatures_enrichment.xlsx`.

Figures in `results/figures/TCGA_senescence_enrichment/`:
`04_HALLMARK_NES_heatmap.png` (NES of 17 tumors and 5 signatures) and
`04_NES_correlation_heatmap.png` (NES correlation per cancer x signature and
collection).

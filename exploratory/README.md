# Exploratory analyses (Roberto Avelar)

Exploratory work on top of the TCGA-senescence overlap pipeline. Nothing here changes
the main pipeline outputs; data and results stay local (`data/` and `results/` are
gitignored), so these scripts document what was run rather than ship results.

## Two-sided overlap tests

- `check_fisher.R`: `fisher.test(alternative = "greater")` reproduces the GeneOverlap
  p-values exactly; the two-sided test also detects depletion.
- `run_both.sh`: runs the overlap pipeline (with `fisher_alternative = "two.sided"`) on
  the full signatures and on signatures with Hallmark proliferation genes removed
  (E2F targets, G2M checkpoint, mitotic spindle, MYC targets v1/v2).
- `summarise.py`: enriched/depleted cancer counts per signature × overlap class.
- `explore.R`: recurrent overlap genes (≥12/17 cancers) and Hallmark content per class.
- `ois_partition.py`: OIS genes split by sharing with RS/SIPS or quiescence, vs tumour DEGs.

## RAS vs OIS (`ras_vs_ois/`)

Asks which OIS genes need the senescent arrest and which follow RAS/BRAF signalling
regardless of arrest, and which of these tumours share.

1. `01_recount3_metadata.R`: recount3 gene counts for the OIS studies of the senescence
   meta-analysis, including samples the meta-analysis excluded (p53/RB/p16 knockdowns,
   BRAF bypass, escaped clones, etoposide-induced senescence).
2. `02_contrast_degs.R`: within-study DEGs with the house `find_degs_new()` (DESeq2,
   padj < 0.05, |log2FC| > log2(1.5); counts AUC-scaled to 40M, protein-coding,
   CPM ≥ 1 in ≥ 30% of samples). `house_env.R` sources the house functions read-only.
3. `03_sanity_checks.py`: ERK targets stay up in every RAS/BRAF arm; MKI67 falls less in
   the bypass arms; escaped clones are fully proliferating.
4. `04_classify_and_project.py`: classifies meta-analysis OIS DEGs as RAS-autonomous
   (change kept in escaped clones and ≥1 partial bypass), arrest-coupled (lost) or
   intermediate; tags NF-κB/BRD4 dependence, RAS- vs DNA-damage-specific, early/late;
   tests TCGA tumour DEG overlap per class. `STRICT=1` anchors on the escaped clones.
5. `05_class_examples.py`: example genes per class shared with tumours.
6. `06_inspect_ris_object.R`, `07_score_ris_states.R`: UCell scores of the classes in the
   GSE245807 RAS-senescent single-cell states (house QC kept); state 3 has MAPK output
   at growing level.
7. `08_inspect_prefilter.R`, `09_cycling_vs_arrested_ras.R`, `10_cycling_dose_model.py`:
   GSE245807 RAS cells still cycling vs in G1, each against phase-matched growing cells. House
   QC, integrity and control-like filters kept; the house LMNB1/G1 subset is reproduced
   (773/652 cells) but not applied for this comparison (agreed 2026-11-02).
8. `14_dose_and_kras.R`: class UCell scores against RAS dose (GSE222951, RPE1-hTERT ER:HRAS-G12V
   sorted S/M/L/XL, HTODemux singlets) and an early KRAS time course (GSE202126, HPNE parental vs
   KRAS days 1-5, VASA-seq; plus a run without the 120k UMI cap). House default QC only; no
   senescent-state filters (dose/time course, not a state definition).

## RAS vs MYC (and other oncogenes) (`ras_vs_myc/`)

Asks whether RAS-induced senescence (RIS) features differ by oncogene, and whether RAS-
vs MYC-driven tumours carry different parts of the RIS programme.

1. `01_driver_status.py`: per-sample KRAS/NRAS/HRAS codon 12/13/61 and BRAF V600 hotspots,
   and MYC high-level amplification (GISTIC = 2), from cBioPortal TCGA PanCancer Atlas.
2. `02_tumour_scores.R`: per-tumour scores (mean within-cohort z of VST expression) for the
   RIS gene classes, ERK output, Hallmark MYC targets, proliferation and a stromal proxy.
3. `03_driver_models.py`: within-cohort score ~ RAS-mutant or MYC-amplified (vs neither) + stroma.
4. `04_activity_models.py`: within-cohort score ~ ERK + MYC targets + proliferation + stroma.
5. `05_nis_vs_ris.R`, `06_nis_class_sharing.py`: NOTCH1-ICD- vs RAS-induced senescence in the
   same IMR90 cells (Hoare 2016, SRP062871; runs summed per biological replicate), and which
   RIS gene classes a non-RAS oncogene's senescence shares. No public MYC-induced senescence
   transcriptome with proliferating controls was found (human or mouse).
6. `07_download_paad_skcm.R`: TCGA-PAAD and TCGA-SKCM counts (tumour-only cohorts; SKCM keeps
   metastases, with sample type as a covariate in 03/04).
7. `08_akt_myc_degs.R`, `09_oncogene_class_sharing.py`: GEO raw counts through the house
   `find_degs_new()`: BJ-hTERT myrAKT1 (GSE130099) and HRAS (GSE130100); RPE-1 dox-MYC,
   cyclin E and empty vector, p53 wt and -/- (GSE185512, 120 h); and RIS-class sharing across
   RAS, NOTCH1, AKT and MYC arms with arrest/ERK/MYC read-outs.

Inputs outside this repository (read-only): the published senescence meta-analysis DEG
table and house functions from the senescence_stress project, and the GSE245807 house
objects. Paths are absolute to the author's machine.

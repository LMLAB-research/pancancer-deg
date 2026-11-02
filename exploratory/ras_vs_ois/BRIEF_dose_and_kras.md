# Brief: RAS-dose and KRAS time-course check of the RIS gene classes

For a Claude Code session started in `~/code/academic/pancancer-deg`. Exploratory work for
Roberto Avelar; nothing here goes to GitHub until he pushes it himself.

## Background

Earlier work on branch `roberto/two-sided-overlaps-ras-vs-ois` split the genes that change in
oncogene-induced (RAS-induced) senescence into classes, using RAS-on cells whose arrest was
bypassed or lost:

- **RAS-autonomous**: the change persists when RAS stays on but the arrest is lost.
- **Arrest-coupled**: the change disappears when the arrest is lost.
- **Intermediate**, plus a tag for genes that depend on NF-κB/BRD4 (SASP-like).

Every test so far agrees. RAS-autonomous genes follow RAS/MAPK activity, appear before the
arrest in an ER:RAS time course and revert when 4-OHT is withdrawn. Arrest-coupled genes follow
the arrest. Two datasets were downloaded but not analysed; this brief covers them.

## Questions

1. **RAS dose (GSE222951, Chan et al. "Titration of RAS alters senescent state").** In RPE1-hTERT
   cells with ER:HRAS-G12V, after 6 days of 4-OHT and flow sorting by RAS-reporter level
   (S < M < L < XL), do RAS-autonomous scores rise steadily with dose? Do arrest-coupled scores
   (and loss of proliferation) appear only at the higher levels? Report MAPK output across the
   levels as the positive control.
2. **Early KRAS time course (GSE202126).** In HPNE cells (hTERT pancreatic duct; parental vs KRAS
   at days 1, 2, 3 and 5; VASA-seq single cells), do RAS-autonomous scores and MAPK output rise
   early relative to parental cells? Do arrest-coupled genes stay flat (these cells may never
   arrest)?

## Inputs (all local; nothing needs downloading)

| What | Path (relative to `~/code/academic/pancancer-deg`) |
|---|---|
| GSE222951 RPE1 single-cell, gene expression | `data/exploratory/GSE222951/raw/GSM7948808_rna_{matrix.mtx,features.tsv,barcodes.tsv}.gz` |
| GSE222951 hashtags (feature names give the sorted group, e.g. `monoculture_XL_d6_…`; ignore the co-culture and `unmapped` tags) | `data/exploratory/GSE222951/raw/GSM7948809_hto_{matrix.mtx,features.tsv,barcodes.tsv}.gz` |
| GSE202126 HPNE parental, HPNE KRAS (gene × cell UMI tables; first two columns are gene id and symbol; cell names carry the group, e.g. `HPNE_Kras_day2_sc14`) | `data/exploratory/GSE202126/GSE202126_HPNE_exon_UMI_count_matrix.txt.gz`, `..._HPNE_KRAS_exon_UMI_count_matrix.txt.gz` |
| Gene classes (columns `gene`, `OIS_direction` UP/DOWN, `arrest_class`, `nfkb_brd4`) | `results/exploratory/ras_vs_ois/strict/ois_gene_classes.csv` |
| Proliferation genes (Hallmark E2F, G2M, mitotic spindle, MYC) | `results/exploratory/two_sided_2026-11-01/prolif_genes.txt` |
| Closest worked example: same scores, same libraries, same output style | `.claude/worktrees/two-sided-fisher/exploratory/ras_vs_ois/13_ras_withdrawal_GSE173879.R` |

Do not use `exploratory/ras_vs_ois/12_dose_and_kras_timecourse.R` in the worktree. It is an
incomplete draft; replace it or move it to the Trash.

## Method (match the earlier analyses)

- **Scores.** Per-cell UCell scores with the fork-based helper `score_ucell_parallel()` from
  `~/code/academic/scutils/R/ucell_parallel.R` (source it; read-only). Use the same score sets as
  script 13:
  - MAPK output: DUSP6, SPRY2, SPRY4, ETV4, ETV5, FOSL1, EGR1.
  - Proliferation, capped at about 1,400 genes because UCell's maxRank is 1,500.
  - RAS-autonomous UP and DOWN, intermediate UP, arrest-coupled UP and DOWN, and NF-κB/BRD4-dependent UP.
  - Also report the fraction of cells with MKI67 and with CDKN1A detected.
- **R packages.** Prepend the senescence_stress single-cell library:
  `.libPaths(c("/Users/ro/Academic/10_Research/senescence_stress/modules/scrna_seq/GSE226225/final_analysis/.Rlib/main", .libPaths()))`.
  It has Seurat and UCell. Run with `Rscript --vanilla`.
- **GSE222951 cell calling.** Use Seurat `HTODemux` singlets with `positive.quantile = 0.99`, as in
  the house staging script, after CLR-normalising the HTO assay. Keep only the S/M/L/XL
  monocultures.
- **QC.** House default (Wechter) cutoffs: at least 300 genes, 1,000–120,000 UMI, and less than
  12% mitochondrial reads. Do **not** apply the house senescent-state filters (control-like kNN,
  LMNB1, G1). This is a dose and time-course comparison, not a definition of a senescent state,
  and those filters would remove the cells of interest. State this in the outputs.
- **Statistics.**
  - Median score per group.
  - Spearman correlation of each score with group order (S→XL; parental→day 5).
  - For the time course, the fraction of the day-5 change already reached at each day.
  - Wilcoxon tests with BH correction, for XL vs S and for day 5 vs parental.
- **Seeds.** `set.seed(1)`. No socket-based parallelism (sandbox): serial, or `mclapply` only.

## Rules

- **senescence_stress is read-only.** It belongs to another, unpublished paper: read its
  libraries and code, never write there.
- **Save scripts in the git worktree** at
  `~/code/academic/pancancer-deg/.claude/worktrees/two-sided-fisher/exploratory/ras_vs_ois/`
  (branch `roberto/two-sided-overlaps-ras-vs-ois`), numbered `14_…`. Run them from the repo
  root, `~/code/academic/pancancer-deg`; the paths in existing scripts assume that.
- **Write outputs** to `results/exploratory/ras_vs_ois/dose_and_kras/` (gitignored, local only).
- **Commit** the new script(s), and a line added to `exploratory/README.md`, to the branch with
  `git -C <worktree> commit`. Do not push and do not merge into main; Roberto pushes.
- **Never permanently delete files.** Moving to the Trash is the most you may do.

## What to report back

- Cells per group after QC.
- A table of median scores per group, with the Spearman rho for each score.
- One paragraph answering each question.
- Caveats: RPE1 and HPNE are epithelial, while the classes were built in fibroblasts; GSE222951
  has no untreated monoculture, so it tests dose within RAS-on cells; HPNE-KRAS may not arrest.
- The commit hash.

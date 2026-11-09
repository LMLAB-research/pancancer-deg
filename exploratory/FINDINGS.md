# Exploratory findings: two-sided overlaps, RAS vs OIS, RAS vs MYC

Roberto Avelar, 2026-10-31 to 2026-11-03. Exploratory; not part of the main pipeline. Script
index: `README.md`. Summary tables behind the numbers below are in `summary_tables/`
(regenerate with `export_summaries.py`); full results and downloaded data are local only
(`results/exploratory/`, `data/exploratory/`; both gitignored). Each choice behind these
results, who made it and any sensitivity check: `DECISIONS.md`. Figures: `figures/` (list below).

## 1. Two-sided overlap tests

**Question.** The overlap pipeline used one-sided (enrichment) Fisher tests, so depleted overlaps
were reported as N.S.

**Change.** `fisher_alternative` in `TCGA_senescence_overlapping_pipeline/config.R`
(`"two.sided"` on this branch). `"greater"` reproduces the GeneOverlap p-values exactly
(`check_fisher.R`). The heatmap uses a diverging scale: blue for enrichment, red for depletion.

**Result** (17 cancers; `results/exploratory/two_sided_2026-11-01/summary.txt`):

| Overlap (tumour vs signature) | OIS, RS, SIPS, CICQ, SSCQ |
|---|---|
| Tumour up / signature down | enriched in 15–17 cancers for every signature |
| Tumour down / signature up | enriched in 14–17 cancers for every signature |
| Same direction, CICQ and SSCQ | depleted in 10–15 cancers |
| Same direction, OIS | enriched in 16–17 cancers |

- The opposite-direction enrichment holds after removing about 670 Hallmark proliferation genes
  (16–17/17).
- OIS genes shared with tumours that are also shared with RS/SIPS overlap tumour-up genes more
  (median log2 OR 0.94) than OIS-only genes (0.57) or genes shared with quiescence (0.26)
  (`ois_partition.py`).

## 2. RAS-autonomous vs arrest-coupled OIS genes (`ras_vs_ois/`)

**Question.** In the meta-analysis, every OIS sample is HRAS-G12V or BRAF-V600E, so RAS signalling
and the senescent arrest are confounded. Which OIS genes need the arrest?

**Method.**
- recount3 samples from the meta-analysis OIS studies, including samples the meta-analysis
  excluded (scripts 01–02).
- Within-study DEGs from the house `find_degs_new()`.
- Each meta-analysis OIS DEG (`arrest_degs_final.csv`) was classified (script 04, strict rule):
  - **RAS-autonomous**: the change is kept in BJ-hTERT RAS clones that escaped senescence (fully
    proliferating), plus ≥1 partial bypass (IMR90 shp53/RB; TIG3 BRAF bypass).
  - **Arrest-coupled**: the change is lost there.
  - **Intermediate**: neither.
- Genes were also tagged for NF-κB/BRD4 dependence (shp65, shBrd4, JQ1 arms of SRP065206).

**Classes.**

| | RAS-autonomous | Intermediate | Arrest-coupled |
|---|---|---|---|
| OIS-up genes | 481 | 683 | 37 |
| OIS-down genes | 511 | 899 | 263 |

ERK targets stay up in every bypass arm.

**Overlap with TCGA tumour DEGs.**

| Class | Tumour overlap (of 17 cancers) |
|---|---|
| RAS-autonomous up | enriched in 14 |
| RAS-autonomous down | enriched in 16 |
| Arrest-coupled up | enriched in 2 |
| Arrest-coupled down | depleted in 15 |
| NF-κB/BRD4-dependent up | enriched in 3 |

**Tests of the classes.**

| Test | Data | Result |
|---|---|---|
| MAPK switched off in arrested cells | GSE245807 house final state 3 (MAPK −2.3 log2) | RAS-autonomous up −0.6 to −0.8; RAS-autonomous down +0.85 to +1.18; arrest-coupled up −0.15 to −0.27 |
| RAS cells still cycling | GSE245807, 621 cells; house G1 filter not applied (agreed 2026-11-02) | MAPK output 34% of G1 RAS cells; RAS-autonomous up retained 0.37; arrest-coupled 0.86–1.5 |
| ER:RAS time course | GSE115301 Smart-seq2 | Day 2 (cycling): MAPK doubled, RAS-autonomous about half of the final change; arrest-coupled and NF-κB/BRD4 genes unchanged until day 4 |
| RAS withdrawn after senescence | GSE173879, house pipeline; 4-OHT kept vs DMSO at harvest | MAPK back to baseline; RAS-autonomous up 30% of the change remains, NF-κB/BRD4 29%; arrest-coupled down 113% |
| RAS dose | GSE222951 RPE1, sorted S < M < L < XL | RAS-autonomous up Spearman rho 0.69 with dose; arrest-coupled down and MKI67 fall only at L and XL |
| Early KRAS | GSE202126 HPNE, days 1–5 | Day 1: RAS-autonomous up at 35% of the day-5 change, arrest-coupled down at 5%; cells slow from day 2 |
| Without the BRAF bypass arm | script 04, `NO_BRAF=1` | 91% of genes keep their class; tumour overlaps unchanged; arrest-coupled up shrinks from 37 to 15 genes |

**Statistics on the single-cell tests** (house pairwise Wilcoxon on cells, BH; `stats_sc_tests.R`; cells are the
test unit):

- **GSE115301:** Growing vs day 2 is significant for MAPK output, RAS-autonomous UP and RAS-autonomous DOWN
  (BH ≤ 10⁻³). It is not significant for arrest-coupled UP, arrest-coupled DOWN or NF-κB/BRD4-dependent UP; those
  first differ at day 4.
- **GSE222951:** 33/36 pairwise tests BH < 0.05. Not significant: arrest-coupled DOWN and UP for S vs M, and
  arrest-coupled UP for L vs XL.
- **Bootstrap 95% CIs:**
  - GSE173879 withdrawal (fraction remaining): RAS-autonomous UP 0.30 (0.25–0.36), NF-κB/BRD4 0.29 (0.25–0.34),
    arrest-coupled DOWN 1.13 (0.92–1.41).
  - GSE245807 cycling/G1: RAS-autonomous UP 0.37 (0.31–0.44), arrest-coupled DOWN 1.51 (1.20–1.91).
- **Clusters:** low-MAPK clusters, each with lower RAS-autonomous UP and NF-κB/BRD4 scores (BH ≤ 3×10⁻⁵ vs the
  other clusters):

  | Cluster | MAPK output (median) |
  |---|---|
  | GSE245807 house final state 3 | 0.019 (others about 0.10) |
  | GSE245807 pre-filter cluster 1 | 0.015 (other 0.095) |
  | GSE173879 4-OHT-kept cluster 2 | 0.004 (other 0.078) |

  After 4-OHT withdrawal, both GSE173879 clusters have median MAPK output 0.
- **Tumour-level paired tests:**
  - Fig. 1: RAS-autonomous vs arrest-coupled overlap, BH p = 0.007 (UP) and 3×10⁻⁵ (DOWN).
  - Fig. 4: ERK minus MYC coefficient, paired across 20 cohorts, BH ≤ 8×10⁻⁶ for RAS-autonomous UP,
    NF-κB/BRD4 UP and RAS-autonomous DOWN; reversed (MYC > ERK) for arrest-coupled UP (BH p = 0.001).

## 3. RAS vs other oncogenes (`ras_vs_myc/`)

**Question.** Is "OIS" specific to RAS?

**In vitro.** Fraction of each class changed in the same direction (script 09):

| Arm | Arrest | Arrest-coupled down | RAS-autonomous up |
|---|---|---|---|
| IMR90 RAS (Hoare 2016) | yes | 0.83 | 0.74 |
| IMR90 NOTCH1-ICD | yes | 0.80 | 0.17 |
| IMR90 RAS (array, GSE45276) | partial | 0.46 | 0.75 |
| IMR90 AKT (array, GSE45276) | no | 0.06 | 0.17 |
| BJ-hTERT HRAS (GSE130100) | no | 0.03 | 0.61 |
| BJ-hTERT myrAKT1 (GSE130099) | partial | 0.59 | 0.31 |
| RPE-1 dox-MYC 120 h (GSE185512) | no | 0.01–0.04 | 0.14–0.20 |

No public MYC-induced senescence transcriptome with proliferating controls was found (human or
mouse; search 2026-11-01).

**Tumours.**
- Driver calls: KRAS/NRAS/HRAS codon 12/13/61 and BRAF V600 hotspots, and MYC GISTIC = 2, from
  cBioPortal TCGA PanCancer Atlas.
- Scores: per-tumour z-scores of the class gene sets in 20 cohorts.
- Models adjusted for a stromal proxy and ABSOLUTE purity (PanCanAtlas).
- RAS-mutant vs neither: RAS-autonomous up is higher in 10/11 cohorts (6 at BH < 0.05; median
  +0.27 SD). PAAD +0.50.
- MYC-amplified vs neither: median +0.14 (1/11 significant).
- Continuous, mutually adjusted coefficients:

| Class | ERK output | MYC targets | Proliferation |
|---|---|---|---|
| RAS-autonomous up | +0.38 (20/20 cohorts) | | |
| NF-κB/BRD4-dependent up | +0.28 (20/20) | | |
| RAS-autonomous down | | −0.45 (20/20) | |
| Arrest-coupled up | | | −0.52 (18/20) |

## 4. TGF-β, cGAS–STING and other pathways (2026-11-07)

**Question.** Do TGF-β, cGAS–STING and other senescence-linked pathways behave like the RAS-autonomous or the
arrest-coupled genes?

**Method.** Pathway sets (`ras_vs_ois/pathway_sets.R`, msigdbr):
- TGF-β: Hallmark TGF-β signalling.
- cGAS–STING: KEGG cytosolic DNA sensing and Reactome STING-mediated induction.
- Type-I IFN output: Hallmark IFN-α.
- Also IFN-γ, TNFα/NF-κB, IL6–JAK–STAT3, p53, NOTCH, and KRAS-up as a positive control.

These sets were used in three ways:
- Which class holds the pathway's genes: house `overlap_function()`, background = classified OIS DEGs of that
  direction (`15_pathway_class_overlap.R`).
- UCell scores in the existing single-cell tests (scripts 07, 09, 11, 13, 14 rerun; existing per-cell scores
  reproduced exactly).
- TCGA: pathway scores as outcomes of the activity/driver models, and as extra covariates in the class models
  (`ras_vs_myc/11_pathway_models.py`). For TCGA, pathway sets exclude stromal-proxy and RIS class genes.

**Results.**
- **Genes.** Few pathway genes are OIS DEGs: 8 TGF-β genes are OIS-up and 4 cGAS–STING genes are OIS-up (IL1B and
  POLR3A RAS-autonomous). None of TGF-β, cGAS–STING, IFN-α/γ, TNFα/NF-κB, IL6–STAT3, p53 or NOTCH prefers a class at BH < 0.05. The
  Hallmarks that do are E2F/G2M/MYC (arrest-coupled down) and EMT, myogenesis and apical junction (RAS-autonomous
  down). As many IFN-α genes are OIS-down (12; 10 intermediate) as OIS-up (12).
- **Single cell** (summary: `summary_tables/08_sc_pathway_summary.csv`):

| Pathway | Time course: senescent vs growing | Fraction remaining after RAS withdrawal | Cycling RAS cells (retained) | RAS dose rho |
|---|---|---|---|---|
| KRAS-up (control) | up, p 2×10⁻¹² | 0.18 | 0.66 | 0.73 |
| TNFα/NF-κB | up, p 2×10⁻¹⁵; 17% by day 2 | 0.21 | 0.64 | 0.71 |
| cGAS–STING (KEGG) | up, p 8×10⁻¹⁶; 23% by day 2, 27% by day 4 | 0.00 | 1.02 | 0.66 |
| p53 | up, p 2×10⁻⁵; 4% by day 2, 94% by day 4 | 0.92 | 0.46 | 0.19 |
| TGF-β | not changed (p 0.08) | (no OIS change) | (no OIS change) | 0.08 |
| IFN-α | not changed (p 0.40) | (no OIS change) | (no OIS change) | 0.24 |
| NOTCH | small (p 0.05) | (no OIS change) | (no OIS change) | 0.02 |

  - cGAS–STING (KEGG; the set includes NF-κB and IL1B/IL6/CCL5 outputs) behaves like a RAS-dependent programme:
    it reverts on withdrawal, tracks dose and is fully present in cycling RAS cells, but it rises late in the time
    course.
  - p53 behaves like an arrest-coupled programme.
  - TGF-β and type-I IFN are not induced cell-autonomously in these OIS datasets.
  - GSE245807 state 3 scores lower than states 1/2 for every pathway, including TGF-β and NOTCH (Wilcoxon
    p ≤ 6×10⁻⁴ vs state 1). So it shows no sign of TGF-β- or NOTCH-driven secondary senescence (open item below).
    Every score is lower in state 3, so part of this may be lower detection.
- **Tumours.**
  - The TGF-β score rises with ERK output (median b +0.39, 20/20 cohorts) and is higher in RAS-hotspot tumours
    (median +0.26 SD; 3/11 cohorts BH < 0.05). cGAS–STING (KEGG) rises with MYC targets (+0.32, 15/20).
  - With TGF-β, cGAS–STING, IFN-α and TNFα/NF-κB added to the class models, RAS-autonomous UP still follows ERK
    (+0.21, 18/20; was +0.38). It also follows TNFα/NF-κB (+0.55, 20/20) and cGAS–STING (+0.22, 19/20).
  - Arrest-coupled UP falls with TGF-β (−0.28, 16/20) and proliferation (−0.48, 18/20). Same with purity.
- **Perturbation data.** No public OIS transcriptome blocks cGAS–STING or TGF-β (search 2026-11-07). The closest
  are Dou 2017 shcGAS in etoposide senescence (GSE99028, recount3) and paracrine OIS co-culture (Acosta 2013,
  GSE41318, arrays). Neither was run.

## 5. Metabolism: senescence vs cancer (`metabolism/`)

**Method.** `01_metabolic_overlap_tcga.R`: a copy of the house `19_metabolic_overlap_v3.R` with new inputs only.
- Same 14 SI Table 37 pathway classes, house `overlap_function()`, SI37 background (SI18 sensitivity) and
  `add_corr()` BH.
- Applied to TCGA tumour DEGs, one call per cancer, and to the 21 per-inducer arms (`degs_nongenotoxic.csv`).
- Senescence results are read from the house output.
- Within metabolic genes, the script also tests whether senescence DEGs overlap tumour DEGs in the same or
  opposite direction.

**Results** (primary background; counts are cancers of 17 at BH < 0.05):

| Pathway class | Senescence (meta-analysis) | Inducer arms | Tumours |
|---|---|---|---|
| Nucleotide synthesis | down in SIPS, CICQ, SSCQ | down in 10 arms (genotoxic, CDK4/6i, quiescence, RAS) | **up in 10** |
| Folate cycle | down in RS | down in IR, bleomycin, palbociclib | **up in 4** |
| Eicosanoids (COX) | up in SIPS, CICQ | up in IR (GSE259295) | **down in 7** |
| Eicosanoids (LOX) | n.s. | n.s. | down in 6 |
| Glycolysis | down in CICQ, SSCQ only | down in doxorubicin 24 h | up in 3 |
| Hexosamine | up in RS, OIS | up in RAS, doxorubicin, etoposide, bleomycin | mixed (2 up, 2 down) |
| Sphingolipids | up in CICQ, SSCQ | up in abemaciclib, palbociclib | n.s. |
| OXPHOS | under-represented in 7 of 10 group × direction sets | under-represented among down DEGs in 13 arms | under-represented in 6 (down) and 9 (up); over in 2 |

- **Opposite to tumours:** nucleotide synthesis, the folate cycle and COX eicosanoids change in opposite
  directions in senescence and tumours.
- **Senescence-specific:** hexosamine (OIS, RS, DNA damage) and sphingolipids (quiescence, CDK4/6i) are not changed
  in tumours.
- **Shared:** OXPHOS genes are rarely DEGs in either.
- **Gene level** (metabolic genes only), matching the genome-wide result: CICQ, SSCQ and SIPS metabolic DEGs overlap
  tumour DEGs in the opposite direction (4–9 cancers per direction), while OIS-down overlaps tumour-down (8 cancers).
  Examples of OIS-up genes also up in tumours: HK2, GPI, ENO2, BPGM, KYNU, CHST11.
- SI18 sensitivity background: same directions (nucleotide synthesis up in 11 cancers; OXPHOS under in 8–9).

**Follow-up checks (2026-11-07).**
- **Proliferation genes removed** (`NOPROLIF=1`; 18 pathway and 85 background genes removed):
  - Nucleotide synthesis and the folate cycle lose their senescence signal; in tumours, nucleotide synthesis drops
    from 10 cancers to 3 and folate from 4 to 1. That axis was the cell cycle.
  - Still present: glycolysis down in CICQ/SSCQ, COX eicosanoids up in SIPS/CICQ, "hexosamine" up in RS/OIS, OXPHOS
    under-represented, and OIS-down / tumour-down metabolic genes (9 cancers).
- **Gene sets.** `00_metabolic_sets.R` splits the SI37 "hexosamine" class into GAG chain synthesis (40 genes: CHST,
  CHPF, HS*ST, XYLT, EXT, NDST, ...) and the hexosamine core (8 genes). These are scored with SI37 glycolysis and
  Hallmark glycolysis.
- **Single cell** (`metabolism/sc_metabolic_summary.csv`):
  - SI37 glycolysis is up in OIS (p 2×10⁻⁷) and early (79% of the change by day 2). It falls below control after RAS
    withdrawal (−0.71) and is partly retained in cycling RAS cells (0.52), as a RAS-dependent programme would be. But it
    does not scale with RAS dose in RPE1 (rho 0.09) and goes **down** with KRAS in HPNE epithelial cells (p 1×10⁻²⁹).
  - GAG synthesis rises early (78% by day 2) and persists after RAS withdrawal (1.41). It does not scale with dose and
    is unchanged in HPNE.
- **Tumours** (`02_tumour_metabolic_purity.py`; ABSOLUTE purity):
  - SI37 glycolysis follows MYC targets (+0.57, 19/20 cohorts), not ERK (+0.04), and is barely higher in RAS-mutant
    tumours (2/11). It is not related to purity or stroma.
  - GAG synthesis is mainly stromal: stroma +0.57 (20/20), purity −0.12 (14/20 negative). It still follows ERK after
    purity adjustment (+0.21, 16/20) but is not higher in RAS-mutant tumours.
- **Reading.**
  - Glycolysis: quiescence switches it off (opposite to tumours); senescence keeps it on; OIS fibroblasts switch it on
    in a RAS-dependent way. Tumours raise the same genes, but there they follow MYC, and epithelial KRAS cells do not
    raise them. The OIS–tumour glycolysis overlap is a shared state reached by different drivers, not a shared RAS
    programme.
  - GAG synthesis: a senescence (not quiescence) feature that persists without RAS. In tumours it comes mostly from
    the stroma.

## 6. Per-stressor senescence vs tumour DEGs (`stressors/`)

**Method.** `01_inducer_tcga_overlap.R`: house `overlap_function()` (two-sided), each senescence DEG set × each
cancer × the 4 direction pairs, background = genes tested in both, BH across all tests. Sets:
- the 5 meta-analysis groups
- 21 inducer arms (`degs_nongenotoxic.csv`)
- the IR time course (3 cell types × 4/10/20 days, AUC per-cell-type files)
- GSE60340 (p53-mutant LFS fibroblasts: Adriamycin, 5-aza, H2O2, quiescence, RS)

**Results** (`figures/fig_stressor_tcga_overlaps.png`; `summary_tables/10_stressor_tcga_summary.csv`):
- **Opposite-direction overlap is near-universal.** Senescence-up/tumour-down is enriched in 11–17 cancers for every
  one of the 40 sets. Senescence-down/tumour-up is enriched in 9–17 cancers for most sets. It is depleted in IR
  melanocytes (4 d, 10 d), doxorubicin 24 h (PAAF), Adriamycin in p53-mutant cells, CICQ and SSCQ. It is weak in
  HDM201 (71 DEGs) and LFS RS.
- **Same-direction overlap depends on the stressor:**
  - Absent, depleted or mixed:
    - CICQ and SSCQ (depleted in 10–15)
    - palbociclib (IMR90)
    - IR in IMR90 (both datasets)
    - quiescence (IMR90)
  - Present in both directions (11–17 cancers each):
    - RAS (IMR90) and OIS (meta)
    - bleomycin (NHDF, NHLF)
    - doxorubicin (BJ 1 d and 8 d; PAAF 24 h)
    - rapamycin (NHLF)
    - abemaciclib 8 d
    - LFS RS
    - IR in keratinocytes (10/20 d) and melanocytes
  - In IR fibroblasts it grows with time: up/up 8 → 15 cancers and down/down 4 → 12 cancers, from 4 d to 20 d.
- **SIPS pools arms that differ.** The SIPS meta group (up/up 13 cancers; down/down 2, depleted in 6) averages
  genotoxic arms that behave differently. Bleomycin and doxorubicin (BJ) share tumour direction; IR in IMR90 does not;
  etoposide (IMR90) shares only up/up (10).

## 7. Components of the CQ, RS and SIPS signatures (`breakdown/`, 2026-11-09)

**Question.** OIS splits into a RAS-autonomous part (shared with tumours) and an arrest-coupled part (lost in tumours).
Do RS, SIPS and CQ split the same way?

**Method.** `01_cq_rs_sips_components.R`; 108 components (direction-specific), each tested against each cancer's
tumour DEGs with house `overlap_function()`. Background = genes tested in the component's source(s) and in that
cancer. BH across all tests; repeated with proliferation genes removed.
- **Arrest-shared vs senescence-specific:**
  - meta level: RS/SIPS/OIS DEGs also changed in CICQ or SSCQ, or not; CQ DEGs also changed in a senescence group,
    or not
  - within GSE162175 (one study, one control): the IR, doxorubicin, etoposide, palbociclib and RAS arms against its
    quiescence arm
- **p53 programme vs the rest:** p53 programme = DEGs in ≥ 2 of 3 MDM2-inhibitor arms. Separately, doxorubicin DEGs in
  ≥ 2 of 4 p53-wild-type arms, kept or lost in the p53-mutant LFS Adriamycin contrast.
- **Timing:** IR 4 d vs 20 d (3 cell types); doxorubicin BJ 1 d vs 8 d.

**Results** (`figures/fig_component_tcga_overlaps.png`; counts are cancers of 17, BH < 0.05):

| Component | Same direction as tumours | Opposite direction |
|---|---|---|
| CQ-specific (CICQ, SSCQ) | depleted (down 15/8, up 16/14) | 16–17 |
| Arrest-shared part of RS / SIPS / OIS (also in CQ) | down: depleted (14–15) or absent; up: mixed | 16–17 |
| Senescence-specific part of RS (not in CQ) | down 16, up 8 | 7–13 |
| Senescence-specific part of SIPS | down 14, up 16 | 6–14 |
| Senescence-specific part of OIS | down 17, up 17 | 5–9 |
| GSE162175 senescence-specific: doxorubicin / etoposide / RAS | down 14 / 10 / 17 | 8–17 |
| GSE162175 senescence-specific: IR / palbociclib | 1–3 | 15–17 |
| p53 programme (MDM2i, 18–50 genes), every group | 0–5 | 5–16 |
| IR fibroblasts: early-only (4 d) | down 0 (4 depleted) | 15–17 |
| IR fibroblasts: late-only (20 d) | down 16, up 15 | down 0, up 16 |
| IR melanocytes: early-only (4 d) | 17 | depleted (12–13) |

- **The OIS pattern holds for RS and SIPS.**
  - The part of each senescence signature that quiescence shares, and the CQ-specific part, are opposite to tumours
    (enriched in 16–17 cancers in the opposite direction; same direction depleted).
  - The senescence-specific part of RS, SIPS and OIS shares tumour direction (down/down in 14–17 cancers).
  - Not a proliferation effect: same results with proliferation genes removed.
- **Within GSE162175 it depends on the stressor:** doxorubicin, etoposide and RAS have a tumour-like
  senescence-specific part; IR and palbociclib in IMR90 do not.
- **The p53 programme is opposite to tumours** in every group (its down genes are tumour-up in 16 cancers) and never
  shares tumour direction. The "not p53" remainder carries the shared part.
- **Doxorubicin, p53-dependent vs independent:** both share tumour direction (9–14 cancers), so the shared part is not
  simply the p53-independent part.
- **Timing:** in IR fibroblasts the early-only changes are opposite to tumours, while changes that appear only by
  day 20 share tumour direction (fibroblasts 15–16 cancers, keratinocytes 14–17). IR melanocytes are tumour-like from
  day 4.
- **What the parts contain** (house Hallmark + lysosomal):
  - Arrest-shared down: E2F, G2M, mitotic spindle, MYC targets.
  - Arrest-shared up: p53 pathway, inflammatory response, IFN, lysosomal genes.
  - Senescence-specific down: EMT, myogenesis, apical junction, KRAS-down.
  - Senescence-specific up: TNFα/NF-κB, KRAS-up, inflammatory response, coagulation (RS: heme metabolism, estrogen
    response, lysosomal).
- **Lysosomal genes** are over-represented among up genes of every arrest type: RS, SIPS, OIS, CICQ and SSCQ, in both
  the arrest-shared and the senescence-specific parts, and outside the p53 programme.

## 8. House stress-response modules vs tumours (`stressors/`, 2026-11-09)

The acute and chronic stress modules (`external_inputs/stress_modules/v3`, "full") are induced (up-regulated) gene
programmes. They were added to `stressors/01_inducer_tcga_overlap.R` as their own BH family; the earlier stressor
results are unchanged.

| Module | Tumour-down (opposite) | Tumour-up (same) |
|---|---|---|
| 2-deoxyglucose (chronic) | 17 (log2 OR 1.33) | 9 (0.39) |
| Glucocorticoid (chronic) | 17 (1.44) | 10 (0.38) |
| Heat shock (senescent / acute) | 16 / 15 | 7 / 9 |
| NF-κB (acute) | 16 | 7 |
| UPR (acute / senescent) | 12 / 14 | 1 / 9 |
| Antimycin A, 5 d (chronic) | 13 (0.60) | 11 (0.37) |
| IL-1β, 14 d (chronic) | 15 (0.88) | 16 (0.96) |
| 21% O2 (chronic) | 16 (0.89) | 17 (0.78) |
| TNF, 24 h (readout) | 11 (0.31) | 15 (0.80) |
| ISR (acute, 133 genes) | 3 | 7 |
| Small acute modules (BPTES, chloramphenicol, rotenone, UVB, DNA damage, H2O2: 14–82 genes) | 0–10 | 0–6 |

- Most stress programmes overlap tumour DEGs in both directions; the median odds ratios show which dominates.
- **Opposite to tumours:** metabolic and glucocorticoid stresses (2-DG, glucocorticoid), heat shock, NF-κB and UPR.
- **Shared with tumours at least as much as opposite:** chronic inflammatory and oxidative programmes (IL-1β 14 d,
  21% O2, TNF 24 h).
- **Antimycin A (mitochondrial complex III):** mixed, leaning opposite.

## 9. House Hallmark + lysosomal sets: senescence vs tumours (`hallmarks/`, 2026-11-09)

**Method.** `01_hallmark_overlap_tcga.R`: a copy of the Hallmark family of the house
`17_stress_module_overlap_v3.R` (49 Hallmark + "Lysosomal Genes"; `overlap_function()`; background = all tested genes;
`add_corr()` BH), with TCGA tumour DEGs and the inducer arms as inputs. Senescence results are read from the house
output. Table: `summary_tables/11_hallmark_senescence_vs_tumour_table.csv`.

**Results** (senescence groups with BH < 0.05; cancers of 17):
- **Proliferation (opposite):** E2F, G2M, MYC targets and mitotic spindle are down in all five groups. In tumours
  they are up in 3–16 cancers (E2F and G2M 16) and depleted among tumour-down genes in 10–17.
- **Up in senescence, mostly down in tumours:** TNFα/NF-κB (up in RS/OIS/SIPS; tumour-down 12, depleted among
  tumour-up 7), KRAS-up (tumour-down 13), estrogen response early (15), myogenesis (13), hypoxia (10), xenobiotic
  metabolism (11). Several of these (myogenesis, adipogenesis, estrogen response) probably reflect normal-tissue loss
  in tumour-vs-normal comparisons.
- **Same direction:** glycolysis (up in RS/OIS; tumour-up 11), angiogenesis (tumour-up 6), and partly IFN-α/γ
  (tumour-up 7, but also tumour-down 2).
- **Unchanged in tumours:** lysosomal genes (up in RS/OIS/SIPS/CICQ; under-represented among tumour DEGs in 5–6
  cancers per direction), p53 pathway (up in RS/OIS/SIPS; 1–4 cancers), OXPHOS and protein secretion.

**House sets in the single-cell RAS tests and tumour models** (`summary_tables/11_sc_house_pathway_summary.csv`).
Rule (Claude): a set is RAS-like if less than half of its OIS change remains after RAS withdrawal and it tracks dose
(|rho| > 0.3); it is arrest-like if ≥ 75% remains and < 30% of the change is reached by day 2.
- **RAS-like:** TNFα/NF-κB, inflammatory response, IL6–STAT3, KRAS-up, IFN-γ, hypoxia, allograft rejection.
- **Arrest-like:** proliferation sets (down), p53, apoptosis, coagulation, EMT, fatty-acid metabolism (down) and
  lysosomal genes. Lysosomal genes rise late (−6% of the change at day 2), are still above the kept-RAS level after
  RAS withdrawal (1.67 of the change), and are partly kept in cycling RAS cells (0.81).
- **No clear OIS change:** TGF-β, IFN-α, NOTCH, Hallmark glycolysis, UPR, OXPHOS.
- **In tumours** (within-cohort models):
  - Lysosomal and p53-pathway scores fall with proliferation (−0.52 and −0.53; 17/20 cohorts) and rise with MYC
    targets. Lysosomal is not linked to ERK or RAS mutation; p53 pathway is higher in RAS-mutant tumours (5/11).
  - TNFα/NF-κB follows ERK (+0.28, 17/20).
  - Same with purity.

## 10. The OIS genes that RS / SIPS share: tumour-likeness and drivers (`breakdown/02_shared_layer_drivers.R`, 2026-11-09)

**Sets.** OIS senescence-specific genes (not in CICQ/SSCQ) and RAS-autonomous genes, each split into "shared with
RS/SIPS" (same-direction DEG in RS or SIPS) and "OIS-only":

| Set | Up | Down |
|---|---|---|
| Senescence-specific, shared with RS/SIPS | 559 | 652 |
| Senescence-specific, OIS-only | 604 | 1,162 |
| RAS-autonomous, shared with RS/SIPS | 247 | 180 |
| RAS-autonomous, OIS-only | 234 | 331 |

**A. Tumours** (background = OIS DEGs of that direction; cancers of 17, BH < 0.05; paired Wilcoxon across cancers):
- **Up genes:** the shared part is the tumour-like part. Shared senescence-specific up genes are tumour-up in 14
  cancers (median log2 OR 0.53); OIS-only up genes in 3 (4 depleted; −0.10); paired BH 0.001. RAS-autonomous: same
  direction (shared 5 vs OIS-only 1), not significant.
- **Down genes:** reversed. OIS-only down genes are tumour-down in 15 cancers (0.82) vs 3 for shared (0.21); paired
  BH 0.001. RAS-autonomous down: both 12.
- So what OIS shares with tumours is genes **switched on** that RS and SIPS also switch on, plus genes **switched
  off** only in OIS.

**B. What reproduces the shared genes outside RAS** (fraction of each set changed the same way in independent
contrasts; shared vs OIS-only compared within the same background):

| Contrast group | Shared up / down | OIS-only up / down |
|---|---|---|
| DNA damage, p53 wild-type (doxorubicin IMR90 and BJ 8 d, etoposide, bleomycin, IR GSE259295) | 0.38–0.77 / 0.30–0.67 | 0.04–0.19 / 0.04–0.17 |
| IR in IMR90 (GSE162175) | 0.15–0.20 / 0.23–0.32 | 0.02–0.03 / 0.00–0.03 |
| IR fibroblasts, day 4 → day 20 | 0.26–0.29 → 0.67–0.77 | 0.02–0.03 → 0.11–0.15 |
| DNA damage, p53 mutant: Adriamycin (LFS) | 0.45–0.51 / 0.21–0.29 (n.s.) | 0.26–0.29 / 0.19–0.31 |
| H2O2 (LFS, p53 mutant) | 0.20–0.30 / 0.04–0.06 | 0.06–0.08 / 0.02–0.03 |
| p53 activation without damage (nutlin, RG7388, HDM201) | 0.01–0.10 / 0.01–0.07 | 0.00–0.05 |
| Arrest without damage: palbociclib / quiescence / rapamycin | 0.12–0.32 / 0.04–0.13 / 0.07–0.08 | 0.00–0.07 |
| Abemaciclib 8 d | 0.45–0.56 | 0.23–0.32 |
| Cytokines (house modules, up only): TNF+IL-1β 3 h / TNF 24 h / IL-1β 14 d | 0.26–0.39 / 0.43–0.52 / 0.25–0.37 | 0.18–0.19 / 0.29–0.34 / 0.17–0.18 |
| AKT (BJ) / NOTCH1-ICD / MYC / cyclin E | 0.38–0.50 / 0.20–0.38 / 0.11–0.20 / 0.00–0.01 | AKT 0.14–0.24; MYC higher in OIS-only |

- Hydroxyurea in RPE1 (GSE185512) has almost no DEGs and is uninformative.
- **Circularity:** "shared" is defined by the meta RS/SIPS signatures, so enrichment in damage contrasts is partly
  expected. The informative comparison is between contrast groups, which are independent studies.
- **The contrast groups:**
  - DNA damage reproduces about half of the shared genes, and in p53-mutant cells still induces about half of the
    shared up genes.
  - p53 activation without damage reproduces almost none; arrest without damage (palbociclib, quiescence, rapamycin)
    little.
  - Cytokines (TNF, IL-1β) induce a third to a half of the shared up genes, more than OIS-only.
  - After IR the shared genes accumulate slowly: about 27% at day 4, about 70% at day 20.

**B2. Dependence inside OIS** (SRP065206 / SRP127037 knockdowns; retention = knockdown change / control change):
- **p65 (NF-κB) knockdown:** almost no effect on either set (median retention about 1; 2–9% dependent).
- **BRD4 knockdown:** the shared genes depend on BRD4 more than OIS-only (down: 49% vs 29% dependent; BH 4×10⁻⁴;
  senescence-specific: 57% vs 43%, BH 3×10⁻⁶).
- **JQ1:** the OIS-only genes depend more (31–41% vs 16–27%).
- **IL-1R knockdown:** shared up genes depend more (23% vs 14% senescence-specific; BH 1×10⁻⁵).
- **p53 and p16/p21 knockdown:** no effect.

**Reading.**
- Outside RAS, the shared layer is switched on by DNA damage, largely without p53, and slowly (accumulating over
  weeks), not by p53 activation or by arrest.
- Inflammatory cytokines induce a large part of it, and in OIS part of it depends on IL-1 signalling and BRD4.
- This fits a damage-triggered, secondary inflammatory / autocrine loop.
- But the "NF-κB" part is not shown directly: p65 knockdown in OIS barely changes these genes (knockdown efficiency
  unknown, possible redundancy with other NF-κB / AP-1 factors).

## Figures (`figures/`)

Full legends (data, n, statistics, scaling): `figures/FIGURE_LEGENDS.md`, which also lists the
`stats_*.png` test figures.

| File | Shows | Script |
|---|---|---|
| `fig1_tcga_overlap_by_class.png` | log2 OR per cancer for each OIS gene class vs tumour DEGs (section 2) | `figures_summary.R` |
| `fig2_oncogene_arms.png` | fraction of each class reproduced by each oncogene arm, with MKI67 log2FC (section 3) | `figures_summary.R` |
| `fig3_single_cell_time_dose_tests.png` | GSE115301 time course, GSE222951 dose, GSE202126 KRAS days (scaled medians); GSE173879 RAS withdrawal and GSE245807 cycling-cell retention (section 2) | `figures_summary.R` |
| `fig4_tcga_activity_models.png` | per-cohort coefficients of ERK, MYC targets and proliferation for each class (section 3) | `figures_summary.R` |
| `fig5_tcga_driver_models.png` | RAS-hotspot and MYC-amplified vs neither, per cohort (section 3) | `figures_summary.R` |
| `umap_GSE245807_RAS_final_states.png` | house final RAS states (clusters, plate) and class scores | `figures_sc_umaps.R` |
| `umap_GSE245807_RAS_cycling_vs_G1.png` | RAS cells before the LMNB1/G1 subset: group (cycling / G1) and class scores | `figures_sc_umaps.R` |
| `umap_GSE173879_OISON.png`, `umap_GSE173879_OISOFF.png` | house final OIS cells with 4-OHT kept / withdrawn: clusters and class scores | `figures_sc_umaps.R` |

UMAPs use the house calls (`DimPlot` clusters first, no cluster labels; `customSignaturePlot` for
scores, scaled 0–10 within each plot). In `umap_GSE245807_RAS_cycling_vs_G1.png`, most cycling RAS
cells form a separate cluster with low MAPK output, low RAS-autonomous UP and high RAS-autonomous
DOWN scores; some house-final G1 cells fall in the same cluster.

## Departures from house practice

- **Metabolism and stressors (sections 5–6):** house functions throughout. The per-cancer and per-arm calls use a
  per-source background (as the house module analyses do). The stressor summary figure is not a house plot.

- **Cycling RAS cells (GSE245807):** the house LMNB1/G1 subset was not applied. Roberto agreed on
  2026-11-02. The rebuild reproduces the house subset exactly (773 and 652 cells).
- **GSE115301 time course:** the authors' published cells were used, without house state filters.
- **GSE173879 staging:** a local copy of `stage_inputs.R` (`ras_vs_ois/house_copy/`) changes only
  the output folder and the plate pick. The house one-plate preprocessing ran unchanged.
- **GSE202126:** the primary run uses the house 120,000-UMI cap, which leaves 8–40 cells per group.
  A run without the cap (all 772 cells) gives the same direction for every class.
  [judgment needed: accept the uncapped run as primary?]
- **GSE45276:** an Affymetrix array, so limma was used (house DEG thresholds).
- **senescence_stress:** files, objects and libraries were read only; nothing was written there.

## Caveats

- All in vitro arms used to build the classes are fibroblasts; RPE1, HPNE and the TCGA tumours are
  epithelial.
- The TIG3 BRAF "bypass" cells also express SV40 small T (GEO GSE76605); their controls do not.
- GSE173879 has one clone and one plate per arm; GEO does not state the 4-OHT withdrawal time.
- The tumour results are bulk correlations. MYC amplification is a proxy for MYC-driven tumours.
- The NOTCH1-, AKT- and MYC-arm results come from different cell types and time points.

## Open

- Push the branch and open a pull request (Roberto).
- Identity of GSE245807 state 3: secondary senescence or RAS non-responders. TGF-β and NOTCH scores are lower in state 3
  than in states 1/2 (section 4), which does not support TGF-β/NOTCH-driven secondary senescence.
- Lysosomal genes: an arrest feature common to senescence and quiescence, lost in proliferating tumours; a candidate
  read-out to add to the arrest-coupled genes for detecting senescence in tumours.
- cGAS dependence: run the Dou 2017 shcGAS etoposide arm (GSE99028, recount3) through the class rule.
- `senescence_stress/work/results/auc/temporal/all_time_degs.csv` is corrupt on disk (33 MB, almost all NUL bytes);
  the per-cell-type `degs.csv` files read fine.
- MYC-induced senescence data would need new experiments or a collaboration.

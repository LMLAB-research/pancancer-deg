# Exploratory findings: two-sided overlaps, RAS vs OIS, RAS vs MYC

Roberto Avelar, 2026-10-31 to 2026-11-03. Exploratory; not part of the main pipeline. Script
index: `README.md`. Summary tables behind the numbers below are in `summary_tables/`
(regenerate with `export_summaries.py`); full results and downloaded data are local only
(`results/exploratory/`, `data/exploratory/`; both gitignored).

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

## Departures from house practice

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
- Identity of GSE245807 state 3: secondary senescence or RAS non-responders.
- MYC-induced senescence data would need new experiments or a collaboration.

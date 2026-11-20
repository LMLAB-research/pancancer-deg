# Pathway enrichment: candidate findings (2026-11-17)

**Status:** leads to investigate, not established results.

**Based on:** enrichment module (`analyses/TCGA_senescence_enrichment/`, run of
2026-11-17 at commit `4f0febb`): unweighted preranked GSEA of 17 TCGA tumors
(paired design, shrunken log2FC) and 5 signatures (OIS, RS, SIPS, CICQ, SSCQ),
MSigDB 26.1.1. Tables: `results/tables/TCGA_senescence_enrichment/03_GSEA_*.csv`
and `04_compare_*.csv`. NES > 0: gene set up; NES < 0: down. "Same" or
"opposite": significant (BH < 0.05) on both sides with the same or opposite
sign.

## Context

Tumor and signature NES are negatively correlated across gene sets for
quiescence (CICQ, SSCQ) and stress-induced and replicative senescence (SIPS,
RS), and close to zero for OIS. Proliferation is the core of this opposition:
HALLMARK E2F targets and G2M checkpoint are opposite in all 85 cancer x
signature pairs.

Senescence and quiescence both arrest proliferation. A pathway that tumors
share with senescence but not with quiescence therefore points to a
senescence-specific programme that tumors also use. The three observations
below come from that contrast.

## 1. OIS and tumors share a "growth without division" programme

A block of biosynthesis pathways is up in OIS and in nearly all tumors, but
down in RS, SIPS and quiescence:

| Process | Gene set | OIS | RS | SIPS | SSCQ | Tumors up |
| --- | --- | --- | --- | --- | --- | --- |
| Ribosome biogenesis | REACTOME rRNA modification in the nucleus and cytosol | +3.2 | −2.0 | −4.0 | −3.9 | 17/17 |
| mRNA processing | REACTOME mRNA splicing | +2.6 | −5.3 | −6.7 | −5.0 | 17/17 |
| Growth signaling | HALLMARK MTORC1 signaling | +1.8 | −0.9 | −1.8 | −2.8 | 16/17 |
| Biosynthetic MYC programme | HALLMARK MYC targets V1 | +2.2 | −4.8 | −7.4 | −7.0 | 17/17 |
| ER protein folding | HALLMARK unfolded protein response | +2.9 | −1.4 | −1.6 | −1.9 | 15/17 |

Leading-edge genes shared by tumors and OIS are translation and folding
machinery: CCT chaperonins (CCT2, CCT5, CCT6A, CCT7), EIF2S2, EIF4G1, HSPA5,
HSP90B1, CALR, and mTORC1 targets such as ACLY and SHMT2.

**Interpretation.** OIS cells are arrested, so their cell-cycle pathways are
opposite to tumors, but they keep a RAS-driven anabolic programme (mTORC1, MYC,
ribosomes, ER load) that other arrest types switch off and that tumors share.
The two layers would cancel each other in the overall comparison, which could
explain why OIS looks neutral (NES correlation ≈ 0). This fits known OIS
biology (senescent cells grow large; mTOR drives the SASP; ER stress in OIS),
but a pan-cancer decoupling of biosynthesis from proliferation, shared
specifically with OIS, would be a new angle. It is the pathway-level
counterpart of the "RAS-autonomous" OIS genes in Roberto's exploratory work
(branch `roberto/two-sided-overlaps-ras-vs-ois`).

**Caveat.** Within these gene sets, tumor and OIS log2FC are not more
correlated than genome-wide ($\rho_\text{set} \le \rho_\text{genome}$ in almost
all pairs): the gene sets shift as a block in both systems, which suggests a
shared upstream driver (MYC/mTORC1) rather than the same gene-by-gene changes.
MYC targets and ribosome biogenesis also rise with proliferation in tumors.

**To test.**
- Score each tumor sample for the anabolic module (ribosome biogenesis,
  translation, mTORC1, MYC biosynthesis, UPR) and for proliferation (E2F,
  G2M): do the two dissociate across tumors?
- Does the anabolic score follow RAS/MAPK activity (mutations, ERK targets)
  rather than proliferation?
- In OIS bypass and escape data (Roberto's datasets): is the anabolic module
  kept when the arrest is lost?

## 2. Glycolysis: a senescence-specific feature shared with tumors

HALLMARK glycolysis is up in 13 of 17 tumors and in all three senescence types
(NES +1.7 to +1.9), but down in quiescence (CICQ −2.0, SSCQ −3.2): it is not a
consequence of arrest. Same direction in 39 cancer x senescence pairs, opposite
to quiescence in 26 pairs. Recurrent shared leading-edge genes are partly
hypoxia/HIF targets (HK2, ERO1A, DDIT4, CXCR4, SPAG4) and partly
glycosaminoglycan synthesis (CHPF, CHPF2, B3GAT3). Gene-level support is
moderate ($\rho_\text{set} > \rho_\text{genome}$ in 54% of pairs). Roberto's
exploratory work independently found glycolysis and GAG synthesis as robust
senescence-vs-quiescence contrasts.

**To test.**
- Separate the HIF part from the GAG-synthesis part of the gene set; check
  HALLMARK hypoxia.
- Control for tumor purity (issue #30).

## 3. NF-κB / SASP splits tumors into two groups

HALLMARK TNFα signaling via NF-κB is strongly up in every senescence signature
(OIS NES +7.5) but divides tumors:

- up: ESCA, HNSC, KIRC, READ, STAD, THCA;
- down: BLCA, BRCA, CHOL, KICH, LIHC, LUAD, LUSC, PRAD, UCEC.

This could define SASP-like and SASP-low cancers, but it is probably confounded
by immune and stromal content of the tumors or of their normals (lung and liver
normals are immune-rich).

**To test.** Repeat after purity adjustment (issue #30) before interpreting.

## Summary

- **Most promising:** an anabolic "growth without division" programme
  (ribosome biogenesis, mRNA processing, mTORC1, MYC biosynthesis, UPR) shared
  by OIS and almost all tumors, and switched off by RS, SIPS and quiescence. It
  is a pathway-level shift, and could explain OIS's neutral overall direction.
- **Senescence-specific glycolysis**, shared with tumors and opposite to
  quiescence, partly through HIF and GAG-synthesis genes.
- **NF-κB/SASP split across cancers**, likely confounded by cell composition.

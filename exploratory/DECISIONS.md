# Decision log

Every choice that shaped the results in `FINDINGS.md`, in the order the analyses were run.

- **Decided by:** "Roberto" means he asked for it or approved it in the session. "Claude" means
  Claude Code chose it and Roberto has not yet reviewed it.
- **Check:** a sensitivity run or sanity check where one exists.
- **Dates:** 2026-10-31 to 2026-11-03.

## A. Inputs and setup

| # | Decision | Alternatives | Reason | Decided by | Check |
|---|---|---|---|---|---|
| A1 | Senescence signatures from the published meta-analysis table `senescence_stress/work/published_reference/arrest_degs_final.csv`, split by `group_1` into OIS, RS, SIPS, CICQ, SSCQ | The later AUC-scaled re-run `work/results/auc/meta_analysis/arrest_degs.csv` | Published version | Claude (offered the AUC version; no reply) | none |
| A2 | Each signature's background = every gene tested for that subtype (about 14.2k) | `common_background.csv` (about 3.4k genes; the shared signature, not a background) | The pipeline needs all tested genes | Claude | none |
| A3 | TCGA DEG pipeline (scripts 01–04) run unchanged; the overlap step reads the covariate-adjusted design | Naive design | Pipeline default (config) | pipeline authors | none |
| A4 | Sandbox proxy fix in `.Rprofile`, active only when `https_proxy` is the credentialed 127.0.0.1 proxy | Home-folder R profile | Keep it project-scoped and off outside the sandbox | Roberto | GDC calls rerun OK; inactive with no proxy or an ordinary proxy |
| A5 | R packages installed with BiocManager into the project renv library at lockfile versions | `renv::restore()` | renv's installer needs a local socket, which the sandbox blocks | Claude | setup check passes; renv itself is 1.3.1, lockfile says 1.2.3 |

## B. Two-sided overlap tests

| # | Decision | Alternatives | Reason | Decided by | Check |
|---|---|---|---|---|---|
| B1 | Two-sided Fisher test (`fisher_alternative`) | One-sided (enrichment only) | Detect depletion | Roberto | `"greater"` reproduces the GeneOverlap p-values exactly (`check_fisher.R`) |
| B2 | Keep GeneOverlap's sample odds ratio; p-value from `fisher.test` | The conditional MLE odds ratio from `fisher.test` | Existing heatmap values stay unchanged | Claude | none |
| B3 | Heatmap on a symmetric diverging scale (blue enrichment, red depletion); OR = 0 drawn at the lower limit | Separate panels | Depletion needs a colour | Claude | none |
| B4 | Kept local, not merged into main | Merge | Roberto: exploratory for now | Roberto | none |
| B5 | Proliferation sensitivity removes the Hallmark E2F, G2M, mitotic spindle, MYC v1 and MYC v2 sets (about 670 genes) from signatures and backgrounds | Other cell-cycle lists | Standard sets in msigdbr | Claude | none |

## C. RAS-autonomous vs arrest-coupled classes (`ras_vs_ois/`)

| # | Decision | Alternatives | Reason | Decided by | Check |
|---|---|---|---|---|---|
| C1 | Use the samples the meta-analysis excluded from its own OIS studies (p53/RB/p16 knockdowns, BRAF bypass, escaped clones, DDIS, time course) | New datasets | Same studies, same cells, matched controls | Claude (plan approved: "do it") | none |
| C2 | Every contrast within one study, against its matched proliferating control in the same genetic background | Pooled controls | Removes study and background effects | Claude | none |
| C3 | DEGs with the house `find_degs_new()` (DESeq2, padj < 0.05, \|log2FC\| > log2(1.5)), sourced read-only from senescence_stress | Other DEG code | House rule: reuse house functions | Roberto (standing rule) | HANDOFF copy identical to `bulk_RNA_seq/R/functions.R` for the four functions used |
| C4 | recount3 counts AUC-scaled to 40M, protein-coding genes; **duplicate gene symbols summed** | House: duplicate symbols dropped | — | Claude | **departure from house practice; not yet tested** |
| C5 | CPM ≥ 1 in ≥ 30% of samples, applied per contrast | House: per condition (kept if it passes in any condition) | Two-group contrasts | Claude | **departure; not yet tested** |
| C6 | SRP062871 technical runs summed per biological replicate | Treat runs as replicates | Runs are not replicates | Claude | 6 per arm; the earlier 12-per-arm run (pseudo-replication) was discarded |
| C7 | Retention = bypass log2FC / reference log2FC; kept if ≥ 0.5, lost if ≤ 0.25; the reference change must exceed log2(1.5) | Other cut-offs | — | Claude | **thresholds not varied** |
| C8 | Escaped clones scored 1 (significant, same direction), 0.5 (same direction, LFC > half threshold) or 0 | Continuous ratio | No matched OIS reference inside SRP172671 | Claude | none |
| C9 | **Strict rule** (used for every result): RAS-autonomous if kept in the escaped clones and ≥ 1 partial bypass; arrest-coupled if lost in the escaped clones and ≥ 1 partial bypass | Lenient rule: ≥ 2 of 3 systems | The lenient rule put 44 proliferation genes in RAS-autonomous, because the knockdown bypasses are partial (MKI67 still −1.6 to −3.5) | Claude | strict rule: 15 proliferation genes in RAS-autonomous; tumour overlaps similar under both rules |
| C10 | NF-κB/BRD4-dependent if retention ≤ 0.5 in ≥ 2 of shp65, shBrd4, JQ1 | — | — | Claude | none |
| C11 | TIG3 BRAF "bypass" arm kept in the main classification | Drop it | Not known at first; GEO later showed hTERT + SV40 small T, with controls lacking small T | Claude | `NO_BRAF=1`: 91% of genes keep their class; arrest-coupled UP 37 → 15 |
| C12 | TCGA projection: tumour DEGs from the adjusted design (padj < 0.05, \|LFC\| > log2(1.5)), protein-coding; universe = tumour genes ∩ OIS-tested genes; two-sided Fisher | — | Matches section B | Claude | none |

## D. Single-cell tests

| # | Decision | Alternatives | Reason | Decided by | Check |
|---|---|---|---|---|---|
| D1 | Scores: UCell via `scutils::score_ucell_parallel` on RNA counts; MAPK output = DUSP6, SPRY2, SPRY4, ETV4, ETV5, FOSL1, EGR1 (the senescence_stress set) | AddModuleScore | House scoring tool | Roberto (house rule) | none |
| D2 | GSE245807 state test on house-final cells; state labels from `notch_secondary/out/step2/cell_scores_GSE245807_RIS.csv.gz` | — | All house filters kept | Claude | none |
| D3 | GSE245807 cycling RAS cells: house LMNB1/G1 subset not applied. Published QC, integrity and control-like filters kept. Phase from the house rebuild (`CellCycleScoring`). Reference: house CTRL_1 cells, phase-matched | Keep the filters | The comparison needs cycling cells | Roberto approved dropping G1 (2026-11-02); LMNB1 also dropped by Claude, with LMNB1+ G1 cells shown as a separate group (16 cells) | rebuild reproduces the house subset exactly (773/652) |
| D4 | Cycling-cell readout: RAS effect in cycling cells / RAS effect in house-final G1 cells, each vs phase-matched growing cells | Model on RAS cells only | The RAS-only model confounds phase with class scores | Claude | both shown; the UMAP shows most cycling RAS cells form a separate low-MAPK cluster |
| D5 | GSE115301 time course: the authors' published (downsampled) cells; no house state filters | House filters | They remove pre-arrest cells by design | Claude | none |
| D6 | GSE173879: OISON = GSM5282311 (in 4-OHT at harvest), OISOFF = GSM5282312 (4-OHT-induced, DMSO at harvest), clone 1; controls CTRL_1/CTRL_2 as in the house RS run | — | GEO sample descriptions | Claude | withdrawal duration not stated on GEO |
| D7 | GSE173879 staging: local copy of `stage_inputs.R` changing only the output root and plate pick; the house one-plate preprocessing ran unchanged | Edit house code | senescence_stress stays read-only | Claude | none |
| D8 | GSE222951 dose: HTODemux singlets (positive.quantile 0.99), house default (Wechter) QC, S/M/L/XL monocultures only, no state filters | — | No untreated monoculture exists, so only a dose comparison within RAS-on cells is possible | Claude (brief) | none |
| D9 | GSE202126: primary run with the house 120,000-UMI cap (8–40 cells per group); uncapped run (772 cells) used in figures | — | VASA-seq median about 156k UMI; the cap targets droplet doublets | second session; [judgment needed: Roberto to approve the uncapped run] | same direction for every class in both runs |

## E. RAS vs other oncogenes

| # | Decision | Alternatives | Reason | Decided by | Check |
|---|---|---|---|---|---|
| E1 | Non-RAS arms: NOTCH1-ICD (SRP062871), AKT (GSE45276 array; GSE130099), MYC dox in RPE-1 at 120 h (GSE185512) | MYC-induced senescence data | None found in human or mouse (search 2026-11-01) | Claude | none |
| E2 | GEO counts parsed from the first 8 fields (descriptions contain unquoted commas) | — | Parsing | Claude | none |
| E3 | GSE45276 (Affymetrix): limma, highest-mean probe per gene, house thresholds | DESeq2 | Array data | Claude | none |
| E4 | GSE185512: dox vs no dox within the same construct at 120 h; the empty-vector dox contrast as a control | 48 h | Latest time point | Claude | empty-vector dox: 0–43 DEGs |

## F. TCGA driver analysis

| # | Decision | Alternatives | Reason | Decided by | Check |
|---|---|---|---|---|---|
| F1 | RAS-pathway = KRAS/NRAS/HRAS missense at codons 12/13/61, or BRAF V600 (cBioPortal PanCanAtlas) | Include NF1 loss, BRAF class II/III, EGFR | Unambiguous hotspots | Claude | ERK output higher in RAS-mutant tumours in 11/11 cohorts (10 at BH < 0.05) |
| F2 | MYC-driven = GISTIC = 2 | Expression-based | — | Claude | MYC targets higher in MYC-amplified tumours in 10/11 cohorts (5 at BH < 0.05) |
| F3 | Groups: RAS-only, MYC-only, neither (both excluded); ≥ 10 tumours per group per cohort | — | — | Claude | none |
| F4 | Scores = mean within-cohort z of VST expression; ten stromal-proxy genes (COL1A1, COL1A2, COL3A1, COL5A1, COL6A3, PDGFRB, FAP, THY1, SPARC, LUM) removed from the class sets and used as a covariate | ssGSEA, singscore | — | Claude | none |
| F5 | Purity covariate: PanCanAtlas ABSOLUTE purity, "called" samples only | ESTIMATE | Genomic estimate | Claude | results with stroma alone and with stroma + purity both reported; same direction |
| F6 | SKCM includes metastases (sample type as covariate); PAAD/SKCM added as tumour-only cohorts | Primary tumours only | Most SKCM samples are metastases | Claude | none |
| F7 | Continuous models: class score ~ ERK output + MYC targets + proliferation + stroma (+ purity); MYC targets = Hallmark MYC v1 + v2; proliferation = E2F + G2M | — | — | Claude | predictors are correlated; coefficients are conditional on each other |

## G. Statistics and figures (2026-11-03)

| # | Decision | Alternatives | Reason | Decided by | Check |
|---|---|---|---|---|---|
| G1 | Single-cell group and cluster comparisons with the house `wilcox_dotplot_wrapper()` (Kruskal–Wallis + all pairwise Wilcoxon, BH) | Pseudobulk by replicate | House function; most datasets have 1–2 replicates, too few for pseudobulk | Roberto asked for statistics; tool by Claude | cells are pseudo-replicates; stated in legends |
| G2 | Ratio read-outs (withdrawal, cycling) get 95% bootstrap CIs: 2,000 resamples of cells within each group, seed 1000 + i / 2000 + i | Delta method | Ratio of medians | Claude | none |
| G3 | Fig. 1: paired Wilcoxon signed-rank between classes across cancers; Fig. 4: one-sample Wilcoxon vs 0 and paired ERK vs MYC across cohorts | Mixed models | Cancer or cohort as the unit | Claude | none |
| G4 | Fig. 5: RAS and MYC in one panel, marker = BH < 0.05 / nominal only / n.s. | Separate rows | Roberto asked | Roberto | none |
| G5 | Pre-filter GSE245807 UMAP recomputed (SCTransform, PCA 50, UMAP dims 1:30, res 0.12, seed 1), not a house-run embedding | House embedding (none exists for these cells) | The house run stores no embedding before the LMNB1/G1 subset | Claude | stated in legend |

## H. Pathways, metabolism, stressors (2026-11-07)

| # | Decision | Alternatives | Reason | Decided by | Check |
|---|---|---|---|---|---|
| H1 | Pathway sets: Hallmark TGF-β; KEGG cytosolic DNA sensing + Reactome STING for cGAS–STING; Hallmark IFN-α as the type-I IFN output; plus IFN-γ, TNFα/NF-κB, IL6–STAT3, p53, NOTCH, KRAS-up (control) | Curated target lists | No Hallmark cGAS–STING set; standard msigdbr sets | Claude (Roberto asked for TGF-β and cGAS–STING "for example") | KRAS-up behaves as RAS-autonomous (positive control) |
| H2 | Class × pathway test: house `overlap_function`, background = classified OIS DEGs of the same direction | All tested genes | Asks which class, not whether OIS | Claude | none |
| H3 | Single-cell scripts rerun with the pathway sets appended | Separate scoring script | One set definition; UCell scores signatures independently | Claude | existing per-cell scores identical (max diff 0) |
| H4 | TCGA pathway scores exclude stromal-proxy and all RIS class genes | Full sets | No shared genes with outcomes or the stroma covariate | Claude | existing model coefficients identical (100/100 rows) |
| H5 | Extended class models add TGF-β, cGAS–STING (KEGG), IFN-α, TNFα/NF-κB as covariates | Pathways as outcomes only | Both reported | Claude | with purity: same direction |
| H6 | Metabolism: copy of house 19_metabolic_overlap_v3.R, inputs only changed; senescence results read from the house output | Recompute | House rule | Claude (Roberto asked to repeat the house analysis) | SI18 sensitivity background gives the same directions |
| H7 | TCGA and inducer-arm metabolic calls: one call per cancer / arm, per-source tested background; BH family = all cancers (or all arms) per bg_type | One pooled background | Tested genes differ by source (as house modules 17/18/19) | Claude | none |
| H8 | Stressor sets: 5 meta groups, 21 nongenotoxic arms, IR time course (AUC per-cell-type files), GSE60340 | AUC `all_time_degs.csv` | That file is corrupt on disk (almost all NUL bytes) | Claude | per-cell-type files have 13.4–13.6k genes and all 3 time points |
| H9 | Stressor × cancer: house `overlap_function`, background = tested in both, BH across all 2,720 tests | BH per set | — | Claude | none |
| H11 | Pathway scores form their own BH family in scripts 03, 07, 09, 14 | One family with the class scores | Existing BH values stay as committed | Claude | every committed summary value identical after the rerun |
| H10 | Stressor summary figure is custom ggplot (median log2 OR + counts), not a house plot | House `summarise_overlaps` | 40 sets × 17 cancers × 4 pairs needs a summary | Claude | full table in results |

| H12 | GAG synthesis vs hexosamine core: gene-name rule on the SI37 "Hexosamine pathway" class (`00_metabolic_sets.R`) | Use the class as is | The class is mostly GAG chain enzymes | Claude | gene lists printed by the script |
| H13 | Metabolic score sets keep RIS class genes in tumour scores (stroma genes removed); own BH family | Remove RIS genes as for the other pathways | Score glycolysis as a whole (HK2 included); not used as covariates | Claude | committed values identical |
| H14 | Proliferation sensitivity for metabolism: the 671-gene list of B5 removed from pathways and backgrounds; senescence arrest results recomputed with the house call | — | Same rule as the main overlaps | Claude (Roberto: "yes do it") | none |
| H15 | Tumour purity check: within-cohort score ~ purity + stroma, ABSOLUTE called samples | ESTIMATE | Same purity source as F5 | Claude | none |

| H16 | House pathway sets (49 Hallmark + Lysosomal Genes, `stress_response_pathways_RERUN.csv`) added as HS_* to the class test, single-cell scores and tumour scores; own BH family | msigdbr Hallmarks only | House rule; lysosomal set was missing | Roberto ("yes do all 3") | class test: same results as msigdbr where sets overlap; committed values identical |
| H17 | Hallmark × tumour DEGs: copy of the Hallmark family of house 17_v3 (`hallmarks/01_hallmark_overlap_tcga.R`), inputs only changed | — | House rule | Claude | none |
| H18 | Components: arrest-shared = same-direction DEG in CICQ or SSCQ (meta) or the GSE162175 quiescence arm; p53 programme = DEG in ≥ 2 of 3 MDM2i arms; doxorubicin p53 dependence = ≥ 2 of 4 wild-type arms vs LFS Adriamycin; timing = 4 d vs 20 d (IR), 1 d vs 8 d (doxorubicin) | Other thresholds | Simplest rules the local data support | Claude | proliferation-removed run: same pattern |
| H19 | Stress modules: "full" version; module genes = DEGs (all up-regulated by construction), background = module tested genes; own BH family | "unique" version | — | Claude | earlier stressor results identical |
| H20 | Single-cell pattern labels (RAS-like / arrest-like) use fixed cut-offs on withdrawal fraction, dose rho and day-2 fraction | Statistical classification | Descriptive summary only | Claude | thresholds not varied |

| H21 | Shared layer = same-direction DEG in meta RS or SIPS; split applied to RAS-autonomous and to OIS senescence-specific genes | RS and SIPS separately | Simplest; both groups pooled | Claude (Roberto asked for the test) | none |
| H22 | Driver test background = the family's OIS genes of that direction tested in the contrast; contrasts from independent studies; etoposide SRP062872 shown but flagged (in meta SIPS) | All genes as background | Compares shared vs OIS-only directly | Claude | none |
| H23 | Dependence = retention ≤ 0.5 (as the NF-κB/BRD4 tag, C7); shared vs OIS-only by Wilcoxon on retention | — | Same rule as the classes | Claude | none |

## Not decided yet

- Accept the GSE202126 uncapped run as primary (D9).
- Test whether summing vs dropping duplicate symbols (C4) and the per-contrast CPM filter (C5)
  change the classes.
- Vary the retention thresholds 0.5 / 0.25 (C7).
- Run the Dou 2017 shcGAS etoposide arm (GSE99028) to tag cGAS-dependent genes.

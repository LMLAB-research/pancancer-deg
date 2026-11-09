"""Copy or derive the summary-level tables behind FINDINGS.md into exploratory/summary_tables/
(tracked in git; the full results stay in the gitignored results/ folder). Per-cell scores,
per-sample tables and full DEG tables are not exported. Run from the pancancer-deg repo root:
    python3 .claude/worktrees/two-sided-fisher/exploratory/export_summaries.py
"""
import csv, os, shutil, statistics
from collections import defaultdict

R = "results/exploratory"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "summary_tables")
os.makedirs(OUT, exist_ok=True)
written = []

def copy(src, dst):
    shutil.copyfile(src, os.path.join(OUT, dst)); written.append(dst)

def write(dst, rows, fields):
    with open(os.path.join(OUT, dst), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
    written.append(dst)

# 1. Two-sided overlaps: per-signature heatmap tables (full signatures and no proliferation genes)
for run, tag in (("full_signatures", "full"), ("no_proliferation_genes", "noprolif")):
    for sig in ("OIS", "RS", "SIPS", "CICQ", "SSCQ"):
        copy(f"{R}/two_sided_2026-11-01/{run}/heatmaps/TCGA_{sig}_heatmap_data.csv",
             f"01_two_sided_overlaps_{tag}_{sig}.csv")
copy(f"{R}/two_sided_2026-11-01/summary.txt", "01_two_sided_overlaps_summary.txt")

# 2. Gene classes and their TCGA overlaps
copy(f"{R}/ras_vs_ois/strict/ois_gene_classes.csv", "02_ois_gene_classes_strict.csv")
copy(f"{R}/ras_vs_ois/strict/tcga_overlap_by_class.csv", "02_tcga_overlap_by_class_strict.csv")
copy(f"{R}/ras_vs_ois/strict_noBRAF/tcga_overlap_by_class.csv", "02_tcga_overlap_by_class_strict_noBRAF.csv")

# 3. Single-cell tests
copy(f"{R}/ras_vs_ois/scrna_GSE245807/state3_vs_states12.csv", "03_GSE245807_state3_vs_states12.csv")
copy(f"{R}/ras_vs_ois/scrna_GSE245807/cycling_vs_arrested_comparisons.csv", "03_GSE245807_cycling_vs_arrested.csv")

def group_medians(src, group_col, dst, skip=()):
    rows = list(csv.DictReader(open(src)))
    scores = [k for k in rows[0] if k not in (group_col, "", "cell", "plate", "Phase", "LMNB1", "cond", "cycling") and k not in skip]
    by = defaultdict(list)
    for r in rows: by[r[group_col]].append(r)
    out = []
    for g, rr in by.items():
        d = {"group": g, "n_cells": len(rr)}
        for s in scores:
            v = [float(x[s]) for x in rr if x[s] not in ("", "NA")]
            d[s] = round(statistics.median(v), 5) if v else ""
        out.append(d)
    write(dst, out, ["group", "n_cells"] + scores)

group_medians(f"{R}/ras_vs_ois/scrna_GSE115301/time_course_cell_scores.csv", "tp", "03_GSE115301_time_course_medians.csv")
group_medians(f"{R}/ras_vs_ois/scrna_GSE173879/withdrawal_cell_scores.csv", "group", "03_GSE173879_ras_withdrawal_medians.csv")
for f in ("GSE222951_RAS_dose_summary.csv", "GSE202126_HPNE_KRAS_timecourse_summary.csv",
          "GSE202126_HPNE_KRAS_timecourse_noUMIcap_summary.csv", "cells_per_group_qc.csv"):
    copy(f"{R}/ras_vs_ois/dose_and_kras/{f}", "03_" + f)

# 4. Oncogene arms: fraction of each class changed in the same direction
C = f"{R}/ras_vs_ois/contrast_degs"
arms = [("IMR90 RAS (Hoare 2016)", "IMR90_N1ICD_RIS"), ("IMR90 NOTCH1-ICD", "IMR90_N1ICD_NIS"),
        ("IMR90 RAS (array)", "IMR90_array_RAS"), ("IMR90 AKT (array)", "IMR90_array_AKT"),
        ("BJ-hTERT HRAS d14", "BJhT_HRAS_d14"), ("BJ-hTERT myrAKT1 d6", "BJhT_AKT_d6"),
        ("RPE-1 p53wt MYC 120h", "RPE1_p53wt_Myc_dox120h"), ("RPE-1 p53-/- MYC 120h", "RPE1_p53ko_Myc_dox120h")]
cls = list(csv.DictReader(open(f"{R}/ras_vs_ois/strict/ois_gene_classes.csv")))
groups = defaultdict(list)
for r in cls:
    if r["arrest_class"] in ("RAS-autonomous", "arrest-coupled", "intermediate"):
        groups[f"{r['OIS_direction']} {r['arrest_class']}"].append(r["gene"])
    if r["OIS_direction"] == "UP" and r["nfkb_brd4"] == "dependent":
        groups["UP NF-kB/BRD4-dependent"].append(r["gene"])
rows = []
for label, name in arms:
    d = {r["gene"].upper(): (float(r["log2FoldChange"]), r["sig"] == "y") for r in csv.DictReader(open(f"{C}/{name}.csv"))}
    row = {"arm": label, "contrast": name}
    for g in ("MKI67", "CDKN1A", "CDKN2A"):
        row[f"log2FC_{g}"] = round(d[g][0], 2) if g in d else ""
    for k, genes in sorted(groups.items()):
        s = 1 if k.startswith("UP") else -1
        gs = [x for x in genes if x in d]
        row[f"frac_same_dir: {k}"] = round(sum(d[x][1] and d[x][0] * s > 0 for x in gs) / max(1, len(gs)), 3)
    rows.append(row)
write("04_oncogene_arms_class_sharing.csv", rows, list(rows[0].keys()))

# 5. Tumour driver and activity models (stroma-adjusted, and stroma + ABSOLUTE purity)
for sub, tag in (("", "stroma"), ("purity_adjusted/", "stroma_purity")):
    copy(f"{R}/ras_vs_myc/{sub}driver_models.csv", f"05_tcga_driver_models_{tag}.csv")
    copy(f"{R}/ras_vs_myc/{sub}activity_models.csv", f"05_tcga_activity_models_{tag}.csv")

# 6. Pathways (TGF-beta, cGAS-STING, IFN, ...; 2026-11-07)
copy(f"{R}/ras_vs_ois/pathways/class_pathway_overlap.csv", "08_ois_class_pathway_overlap.csv")
copy(f"{R}/ras_vs_ois/pathways/pathway_genes_by_class.csv", "08_pathway_genes_by_class.csv")
copy(f"{R}/ras_vs_ois/pathways/sc_pathway_summary.csv", "08_sc_pathway_summary.csv")
for sub, tag in (("", "stroma"), ("purity_adjusted/", "stroma_purity")):
    copy(f"{R}/ras_vs_myc/{sub}pathway_extended_models_summary.csv", f"08_tcga_pathway_extended_models_{tag}.csv")

# 7. Metabolism, senescence vs tumours (2026-11-07)
for f in ("metabolic_overlap_tcga.csv", "metabolic_overlap_inducer_arms.csv", "metabolic_genes_senescence_vs_tumour.csv",
          "metabolic_genes_direction_table.csv"):
    copy(f"{R}/metabolism/{f}", "09_" + f)

for f in ("noprolif/metabolic_overlap_tcga.csv", "noprolif/metabolic_overlap_arrest_noprolif.csv",
          "tumour_metabolic_purity_models.csv", "sc_metabolic_summary.csv", "metabolic_score_sets.csv"):
    copy(f"{R}/metabolism/{f}", "09_" + f.replace("noprolif/", "noprolif_"))

# 9. House Hallmark + lysosomal sets, components, modules (2026-11-09)
copy(f"{R}/ras_vs_ois/pathways/class_pathway_overlap_house.csv", "11_ois_class_house_pathway_overlap.csv")
copy(f"{R}/ras_vs_ois/pathways/sc_house_pathway_summary.csv", "11_sc_house_pathway_summary.csv")
copy(f"{R}/hallmarks/hallmark_senescence_vs_tumour_table.csv", "11_hallmark_senescence_vs_tumour_table.csv")
copy(f"{R}/hallmarks/hallmark_overlap_tcga.csv", "11_hallmark_overlap_tcga.csv")
copy(f"{R}/hallmarks/hallmark_overlap_inducer_arms.csv", "11_hallmark_overlap_inducer_arms.csv")
for f in ("components.csv", "component_tcga_summary.csv", "component_house_pathways.csv"):
    copy(f"{R}/breakdown/{f}", "12_" + f)

for f in ("shared_layer_tcga_summary.csv", "shared_layer_tcga_paired.csv", "shared_layer_perturbations.csv", "shared_layer_ois_dependence.csv"):
    copy(f"{R}/breakdown/shared_layer/{f}", "13_" + f)

# 8. Per-stressor senescence DEGs vs tumour DEGs (2026-11-07)
copy(f"{R}/stressors/stressor_tcga_summary.csv", "10_stressor_tcga_summary.csv")

total = sum(os.path.getsize(os.path.join(OUT, f)) for f in written)
print(f"{len(written)} files, {total / 1e6:.2f} MB -> {OUT}")

"""Continuous version: within each cohort, RIS class score ~ ERK output + MYC targets +
proliferation + stroma (all standardised), across all primary tumours regardless of
mutation calls. Exploratory; run from the pancancer-deg repo root."""
import pandas as pd, numpy as np, statsmodels.formula.api as smf
import sys; sys.path.insert(0, "exploratory"); from paths import *  # repo-relative inputs (paths.py)

OUT = "results/exploratory/ras_vs_myc"
d = pd.read_csv(f"{OUT}/tumour_scores.csv")

import os
PURITY = os.environ.get("PURITY", "0") == "1"
if PURITY:
    # ABSOLUTE purity (TCGA PanCanAtlas mastercalls, GDC file 4f277128-...), called samples only
    pur = pd.read_csv(PURITY_FILE, sep="\t")
    pur = pur[pur["call status"] == "called"][["array", "purity"]].drop_duplicates("array")
    d = d.merge(pur.rename(columns={"array": "sample15"}), on="sample15", how="inner")
    OUT = OUT + "/purity_adjusted"; os.makedirs(OUT, exist_ok=True)
COV = " + Stroma" + (" + purity" if PURITY else "")
targets = ["RIS_RASauto_UP", "RIS_intermediate_UP", "RIS_arrest_UP", "RIS_NFKB_BRD4_UP", "RIS_RASauto_DOWN"]
targets += [c for c in d.columns if c.startswith(("P_", "HS_"))]  # HS_ = house Hallmark + lysosomal (2026-11-09)  # added 2026-11-07: pathway scores as outcomes
preds = ["ERK_output", "MYC_targets", "Proliferation", "Stroma"]
res = []
for coh, g in d.groupby("cohort"):
    g = g.copy()
    for c in targets + preds:
        g[c] = (g[c] - g[c].mean()) / g[c].std()
    for t in targets:
        m = smf.ols(f"{t} ~ ERK_output + MYC_targets + Proliferation" + COV + (" + C(stype)" if g.stype.nunique() > 1 else ""), data=g).fit()
        res.append(dict(cohort=coh, n=len(g), target=t, **{f"b_{p}": m.params[p] for p in preds},
                        **{f"p_{p}": m.pvalues[p] for p in preds}))
r = pd.DataFrame(res)
r.to_csv(f"{OUT}/activity_models.csv", index=False)
print(f"Standardised coefficients (median across {r.cohort.nunique()} cohorts; n cohorts with p<0.05 positive/negative)")
for t in targets:
    s = r[r.target == t]
    line = f"  {t:22}"
    for p in ("ERK_output", "MYC_targets", "Proliferation"):
        b = s[f"b_{p}"]; pv = s[f"p_{p}"]
        line += f" | {p}: {b.median():5.2f} ({((pv < .05) & (b > 0)).sum():2}+/{((pv < .05) & (b < 0)).sum():2}-)"
    print(line)

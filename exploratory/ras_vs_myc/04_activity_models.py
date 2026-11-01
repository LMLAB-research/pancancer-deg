"""Continuous version: within each cohort, RIS class score ~ ERK output + MYC targets +
proliferation + stroma (all standardised), across all primary tumours regardless of
mutation calls. Exploratory; run from the pancancer-deg repo root."""
import pandas as pd, numpy as np, statsmodels.formula.api as smf

OUT = "results/exploratory/ras_vs_myc"
d = pd.read_csv(f"{OUT}/tumour_scores.csv")
targets = ["RIS_RASauto_UP", "RIS_intermediate_UP", "RIS_arrest_UP", "RIS_NFKB_BRD4_UP", "RIS_RASauto_DOWN"]
preds = ["ERK_output", "MYC_targets", "Proliferation", "Stroma"]
res = []
for coh, g in d.groupby("cohort"):
    g = g.copy()
    for c in targets + preds:
        g[c] = (g[c] - g[c].mean()) / g[c].std()
    for t in targets:
        m = smf.ols(f"{t} ~ ERK_output + MYC_targets + Proliferation + Stroma", data=g).fit()
        res.append(dict(cohort=coh, n=len(g), target=t, **{f"b_{p}": m.params[p] for p in preds},
                        **{f"p_{p}": m.pvalues[p] for p in preds}))
r = pd.DataFrame(res)
r.to_csv(f"{OUT}/activity_models.csv", index=False)
print("Standardised coefficients (median across 18 cohorts; n cohorts with p<0.05 positive/negative)")
for t in targets:
    s = r[r.target == t]
    line = f"  {t:22}"
    for p in ("ERK_output", "MYC_targets", "Proliferation"):
        b = s[f"b_{p}"]; pv = s[f"p_{p}"]
        line += f" | {p}: {b.median():5.2f} ({((pv < .05) & (b > 0)).sum():2}+/{((pv < .05) & (b < 0)).sum():2}-)"
    print(line)

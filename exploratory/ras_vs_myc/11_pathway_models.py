"""Do the RIS class scores still follow ERK output once TGF-beta, cGAS-STING, type-I IFN and NF-kB activity are
in the model, and does any class follow those pathways? Within each cohort (all standardised):
  class score ~ ERK + MYC targets + proliferation + TGF-beta + cGAS-STING (KEGG) + IFN-alpha + TNFa/NF-kB + stroma
  (+ purity with PURITY=1; + sample type for SKCM).
Pathway scores come from 02_tumour_scores.R with stromal-proxy and RIS class genes removed. Predictors are
correlated, so each coefficient is conditional on the others. Exploratory; run from the pancancer-deg repo root."""
import os
import sys; sys.path.insert(0, "exploratory"); from paths import *  # repo-relative inputs (paths.py)
import pandas as pd, numpy as np, statsmodels.formula.api as smf
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

OUT = "results/exploratory/ras_vs_myc"
d = pd.read_csv(f"{OUT}/tumour_scores.csv")
PURITY = os.environ.get("PURITY", "0") == "1"
if PURITY:
    pur = pd.read_csv(PURITY_FILE, sep="\t")
    pur = pur[pur["call status"] == "called"][["array", "purity"]].drop_duplicates("array")
    d = d.merge(pur.rename(columns={"array": "sample15"}), on="sample15", how="inner")
    OUT = OUT + "/purity_adjusted"; os.makedirs(OUT, exist_ok=True)
COV = " + Stroma" + (" + purity" if PURITY else "")
targets = ["RIS_RASauto_UP", "RIS_intermediate_UP", "RIS_arrest_UP", "RIS_NFKB_BRD4_UP", "RIS_RASauto_DOWN"]
preds = ["ERK_output", "MYC_targets", "Proliferation", "P_TGFB", "P_cGAS_STING_KEGG", "P_IFNA", "P_TNFA_NFKB", "Stroma"]
res = []
for coh, g in d.groupby("cohort"):
    g = g.copy()
    for c in targets + preds:
        g[c] = (g[c] - g[c].mean()) / g[c].std()
    for t in targets:
        m = smf.ols(f"{t} ~ " + " + ".join(preds[:-1]) + COV + (" + C(stype)" if g.stype.nunique() > 1 else ""), data=g).fit()
        res.append(dict(cohort=coh, n=len(g), target=t, r2=m.rsquared, **{f"b_{p}": m.params[p] for p in preds},
                        **{f"p_{p}": m.pvalues[p] for p in preds}))
r = pd.DataFrame(res)
r.to_csv(f"{OUT}/pathway_extended_models.csv", index=False)

# one-sample Wilcoxon of each coefficient vs 0 across cohorts (cohort = unit), BH across all tests
rows = []
for t in targets:
    s = r[r.target == t]
    for p in preds[:-1]:
        b = s[f"b_{p}"]; pv = s[f"p_{p}"]
        rows.append(dict(target=t, predictor=p, median_b=b.median(), n_pos_sig=int(((pv < .05) & (b > 0)).sum()),
                         n_neg_sig=int(((pv < .05) & (b < 0)).sum()), n=len(b), p_wilcoxon=wilcoxon(b).pvalue))
w = pd.DataFrame(rows); w["padj_BH"] = multipletests(w.p_wilcoxon, method="fdr_bh")[1]
w.to_csv(f"{OUT}/pathway_extended_models_summary.csv", index=False)
pd.set_option("display.width", 250)
print(f"Extended model{' (+ purity)' if PURITY else ''}: median standardised coefficient across {r.cohort.nunique()} cohorts "
      "(n cohorts p<0.05 +/-; * Wilcoxon vs 0 BH<0.05)")
tab = w.assign(cell=lambda x: x.median_b.round(2).astype(str) + " (" + x.n_pos_sig.astype(str) + "+/" +
               x.n_neg_sig.astype(str) + "-)" + np.where(x.padj_BH < .05, "*", "")).pivot(index="target", columns="predictor", values="cell")
print(tab.loc[targets, preds[:-1]].to_string())
print("\nMedian R2:", r.groupby("target").r2.median().round(2).to_dict())

"""Are the tumour glycolysis / GAG-synthesis scores tumour-cell signals or composition? Within each cohort (primary
tumours, ABSOLUTE purity "called" samples), standardised score ~ purity + stroma proxy (+ sample type for SKCM), plus
Spearman rho with purity. A tumour-cell-intrinsic score rises with purity; a stromal one falls with purity / rises with
stroma. Scores from ras_vs_myc/02_tumour_scores.R (metabolic sets keep RIS genes; stroma-proxy genes removed).
The ERK / RAS-driver associations with purity adjustment are in ras_vs_myc/{purity_adjusted/}activity_models.csv and
driver_models.csv (rows P_GLYC_*, P_GAG_SYNTH, P_HEXOSAMINE_CORE). Exploratory, 2026-11-07; run from the repo root."""
import pandas as pd, numpy as np, statsmodels.formula.api as smf
from scipy.stats import spearmanr, wilcoxon

R = "results/exploratory"
d = pd.read_csv(f"{R}/ras_vs_myc/tumour_scores.csv")
pur = pd.read_csv("data/exploratory/tcga_extra/TCGA_mastercalls.abs_tables_JSedit.fixed.txt", sep="\t")
pur = pur[pur["call status"] == "called"][["array", "purity"]].drop_duplicates("array")
d = d.merge(pur.rename(columns={"array": "sample15"}), on="sample15", how="inner")
S = ["P_GLYC_SI37", "P_GLYC_HALLMARK", "P_GAG_SYNTH", "P_HEXOSAMINE_CORE", "Proliferation", "ERK_output"]
rows = []
for coh, g in d.groupby("cohort"):
    g = g.copy()
    for c in S + ["Stroma", "purity"]:
        g[c] = (g[c] - g[c].mean()) / g[c].std()
    for s in S:
        m = smf.ols(f"{s} ~ purity + Stroma" + (" + C(stype)" if g.stype.nunique() > 1 else ""), data=g).fit()
        rows.append(dict(cohort=coh, n=len(g), score=s, rho_purity=spearmanr(g[s], g.purity)[0],
                         b_purity=m.params["purity"], p_purity=m.pvalues["purity"],
                         b_stroma=m.params["Stroma"], p_stroma=m.pvalues["Stroma"]))
r = pd.DataFrame(rows)
r.to_csv(f"{R}/metabolism/tumour_metabolic_purity_models.csv", index=False)
print(f"Within-cohort score ~ purity + stroma ({r.cohort.nunique()} cohorts): median coefficient (n cohorts p<0.05 +/-), "
      "Wilcoxon p across cohorts")
for s in S:
    x = r[r.score == s]
    f = lambda b, p: f"{b.median():+.2f} ({((p < .05) & (b > 0)).sum()}+/{((p < .05) & (b < 0)).sum()}-) W p={wilcoxon(b).pvalue:.1g}"
    print(f"  {s:18} purity {f(x.b_purity, x.p_purity)} | stroma {f(x.b_stroma, x.p_stroma)} | median rho(purity) {x.rho_purity.median():+.2f}")

# ERK / driver associations of the metabolic scores, stroma-only vs stroma + purity
for sub, tag in (("", "stroma"), ("purity_adjusted/", "stroma + purity")):
    a = pd.read_csv(f"{R}/ras_vs_myc/{sub}activity_models.csv"); a = a[a.target.str.match(r"P_(GLYC|GAG|HEXOSAMINE)")]
    dm = pd.read_csv(f"{R}/ras_vs_myc/{sub}driver_models.csv"); dm = dm[dm.score.str.match(r"P_(GLYC|GAG|HEXOSAMINE)")]
    print(f"\n== {tag}: score ~ ERK + MYC + proliferation (median b, n cohorts p<0.05 +/-); RAS / MYC driver effect (median SD, n BH<0.05 +/-)")
    for t, x in a.groupby("target"):
        cell = lambda p: f"{x['b_'+p].median():+.2f} ({((x['p_'+p] < .05) & (x['b_'+p] > 0)).sum()}+/{((x['p_'+p] < .05) & (x['b_'+p] < 0)).sum()}-)"
        dd = {k: f"{v.effect_sd.median():+.2f} ({((v.padj < .05) & (v.effect_sd > 0)).sum()}+/{((v.padj < .05) & (v.effect_sd < 0)).sum()}-, n={len(v)})"
              for k, v in dm[dm.score == t].groupby("driver")}
        print(f"  {t:18} ERK {cell('ERK_output')} MYC {cell('MYC_targets')} prolif {cell('Proliferation')} | RAS {dd.get('RAS')} MYC {dd.get('MYC')}")

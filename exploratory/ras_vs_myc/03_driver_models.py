"""Within-cohort models: score ~ driver group + stroma proxy, groups = RAS-hotspot only,
MYC-amp only, vs neither. Exploratory; run from the pancancer-deg repo root."""
import pandas as pd, numpy as np, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

OUT = "results/exploratory/ras_vs_myc"
d = pd.read_csv(f"{OUT}/tumour_scores.csv")
d = d[d["ras_mut"].notna() & d["myc_amp"].notna()].copy()
d["ras_mut"] = d["ras_mut"].astype(str) == "True"; d["myc_amp"] = d["myc_amp"].astype(str) == "True"
d["group"] = np.select([d.ras_mut & ~d.myc_amp, d.myc_amp & ~d.ras_mut, ~d.ras_mut & ~d.myc_amp],
                       ["RAS", "MYC", "neither"], "both")
scores = ["ERK_output", "MYC_targets", "Proliferation", "RIS_RASauto_UP", "RIS_intermediate_UP",
          "RIS_arrest_UP", "RIS_NFKB_BRD4_UP", "RIS_RASauto_DOWN"]
res = []
for coh, g in d.groupby("cohort"):
    n = g.group.value_counts()
    for drv in ("RAS", "MYC"):
        if n.get(drv, 0) < 10 or n.get("neither", 0) < 10: continue
        s = g[g.group.isin([drv, "neither"])].copy()
        s["is_drv"] = (s.group == drv).astype(int)
        for sc in scores:
            s["y"] = (s[sc] - s[sc].mean()) / s[sc].std()  # standardise within cohort
            m = smf.ols("y ~ is_drv + Stroma", data=s).fit()
            res.append(dict(cohort=coh, driver=drv, n_driver=int(n[drv]), n_neither=int(n["neither"]),
                            score=sc, effect_sd=m.params["is_drv"], p=m.pvalues["is_drv"]))
r = pd.DataFrame(res)
r["padj"] = multipletests(r.p, method="fdr_bh")[1]
r.to_csv(f"{OUT}/driver_models.csv", index=False)

pd.set_option("display.width", 250)
for drv in ("RAS", "MYC"):
    t = r[r.driver == drv]
    tab = t.pivot(index="cohort", columns="score", values="effect_sd")[scores].round(2)
    star = t.pivot(index="cohort", columns="score", values="padj")[scores]
    tab = tab.astype(str) + np.where(star < 0.05, "*", " ")
    ncoh = t.groupby("cohort").n_driver.first()
    tab.insert(0, "n", ncoh)
    print(f"\n=== {drv}-driven vs neither (effect in within-cohort SD, adjusted for stroma; * BH<0.05)")
    print(tab.to_string())
    s = t.groupby("score").apply(lambda x: pd.Series(dict(median=x.effect_sd.median(),
                                                          up=(x.effect_sd > 0).sum() if True else 0,
                                                          sig_up=((x.padj < .05) & (x.effect_sd > 0)).sum(),
                                                          sig_down=((x.padj < .05) & (x.effect_sd < 0)).sum(),
                                                          n=len(x))))
    print(s.loc[scores].round(2).to_string())

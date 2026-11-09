"""Within-cohort models: score ~ driver group + stroma proxy, groups = RAS-hotspot only,
MYC-amp only, vs neither. Exploratory; run from the pancancer-deg repo root."""
import pandas as pd, numpy as np, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

OUT = "results/exploratory/ras_vs_myc"
d = pd.read_csv(f"{OUT}/tumour_scores.csv")

import os
PURITY = os.environ.get("PURITY", "0") == "1"
if PURITY:
    # ABSOLUTE purity (TCGA PanCanAtlas mastercalls, GDC file 4f277128-...), called samples only
    pur = pd.read_csv("data/exploratory/tcga_extra/TCGA_mastercalls.abs_tables_JSedit.fixed.txt", sep="\t")
    pur = pur[pur["call status"] == "called"][["array", "purity"]].drop_duplicates("array")
    d = d.merge(pur.rename(columns={"array": "sample15"}), on="sample15", how="inner")
    OUT = OUT + "/purity_adjusted"; os.makedirs(OUT, exist_ok=True)
COV = " + Stroma" + (" + purity" if PURITY else "")
d = d[d["ras_mut"].notna() & d["myc_amp"].notna()].copy()
d["ras_mut"] = d["ras_mut"].astype(str) == "True"; d["myc_amp"] = d["myc_amp"].astype(str) == "True"
d["group"] = np.select([d.ras_mut & ~d.myc_amp, d.myc_amp & ~d.ras_mut, ~d.ras_mut & ~d.myc_amp],
                       ["RAS", "MYC", "neither"], "both")
scores = ["ERK_output", "MYC_targets", "Proliferation", "RIS_RASauto_UP", "RIS_intermediate_UP",
          "RIS_arrest_UP", "RIS_NFKB_BRD4_UP", "RIS_RASauto_DOWN"]
scores += [c for c in d.columns if c.startswith(("P_", "HS_"))]  # added 2026-11-07: pathway scores
res = []
for coh, g in d.groupby("cohort"):
    n = g.group.value_counts()
    for drv in ("RAS", "MYC"):
        if n.get(drv, 0) < 10 or n.get("neither", 0) < 10: continue
        s = g[g.group.isin([drv, "neither"])].copy()
        s["is_drv"] = (s.group == drv).astype(int)
        for sc in scores:
            s["y"] = (s[sc] - s[sc].mean()) / s[sc].std()  # standardise within cohort
            m = smf.ols("y ~ is_drv" + COV + (" + C(stype)" if s.stype.nunique() > 1 else ""), data=s).fit()
            res.append(dict(cohort=coh, driver=drv, n_driver=int(n[drv]), n_neither=int(n["neither"]),
                            score=sc, effect_sd=m.params["is_drv"], p=m.pvalues["is_drv"]))
r = pd.DataFrame(res)
r["padj"] = r.groupby(np.where(r.score.str.startswith("HS_"), "house", np.where(r.score.str.match(r"P_(GLYC|GAG|HEXOSAMINE)"), "met", np.where(r.score.str.startswith("P_"), "path", "core")))).p.transform(lambda x: multipletests(x, method="fdr_bh")[1])  # pathway scores (2026-11-07): separate BH family
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

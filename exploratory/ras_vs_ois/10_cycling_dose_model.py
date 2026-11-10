"""Within RAS-exposed GSE245807 cells: does each class score follow MAPK output (RAS dose)
or cycling status? score ~ MAPK_output + cycling + plate (all scores standardised).
Run from the pancancer-deg repo root."""
import pandas as pd, statsmodels.formula.api as smf
d = pd.read_csv("results/exploratory/ras_vs_ois/scrna_GSE245807/cycling_vs_arrested_cell_scores.csv")
d = d[(d.cond == "RAS") & (d.group != "RAS_G1_LMNB1pos")].copy()
d["cyc"] = (d.cycling == "cycling").astype(int)
z = lambda s: (s - s.mean()) / s.std()
print(f"RAS cells: {len(d)} ({d.cyc.sum()} cycling). MAPK output, cycling vs G1 (SD units): "
      f"{z(d.MAPK_output)[d.cyc == 1].mean() - z(d.MAPK_output)[d.cyc == 0].mean():.2f}")
print(f"{'score':18} {'b_MAPK':>8} {'p':>9} {'b_cycling':>10} {'p':>9}")
for s in ["RASauto_UP", "intermediate_UP", "arrest_UP", "NFKB_BRD4_dep_UP", "RASauto_DOWN", "arrest_DOWN"]:
    d["y"] = z(d[s]); d["m"] = z(d.MAPK_output)
    f = smf.ols("y ~ m + cyc + C(plate)", data=d).fit()
    print(f"{s:18} {f.params['m']:8.2f} {f.pvalues['m']:9.1e} {f.params['cyc']:10.2f} {f.pvalues['cyc']:9.1e}")

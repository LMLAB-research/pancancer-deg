"""Example genes per class, ranked by how many TCGA cancers share the OIS direction.
Run from the pancancer-deg repo root."""
import csv, glob, os
from collections import Counter
cls = {r["gene"]: r for r in csv.DictReader(open("results/exploratory/ras_vs_ois/strict/ois_gene_classes.csv"))}
OV = "results/exploratory/two_sided_2026-11-01/full_signatures/OIS"
freq = {"UP": Counter(), "DOWN": Counter()}
for d in glob.glob(OV + "/TCGA-*"):
    for k, cl in (("UP", "UP-UP"), ("DOWN", "DOWN-DOWN")):
        for r in csv.DictReader(open(f"{d}/overlap_{cl}.csv")):
            freq[k][r["gene_symbol"]] += 1
for d in ("UP", "DOWN"):
    print(f"\n=== OIS-{d}: genes shared with tumours (n cancers of 17), by class")
    for c in ("RAS-autonomous", "intermediate", "arrest-coupled"):
        g = [(freq[d][x], x) for x, r in cls.items() if r["OIS_direction"] == d and r["arrest_class"] == c and freq[d][x] >= 10]
        g.sort(reverse=True)
        print(f"  {c} ({len(g)} genes in >=10 cancers): " + ", ".join(f"{x}" for n, x in g[:40]))
    if d == "UP":
        g = sorted(((freq[d][x], x) for x, r in cls.items() if r["OIS_direction"] == d and r["nfkb_brd4"] == "dependent"), reverse=True)
        print(f"  NF-kB/BRD4-dependent, top by tumour sharing: " + ", ".join(f"{x}({n})" for n, x in g[:25]))

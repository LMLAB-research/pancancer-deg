"""Which RIS gene classes does a non-RAS oncogene's senescence (NOTCH1-ICD, NIS) share?
Fraction of each class changed in the same direction (house sig) in NIS vs RIS, same cells.
Run from the pancancer-deg repo root."""
import csv
from collections import defaultdict
C = "results/exploratory/ras_vs_ois/contrast_degs"
def load(n):
    return {r["gene"].upper(): (float(r["log2FoldChange"]), r["sig"] == "y") for r in csv.DictReader(open(f"{C}/{n}.csv"))}
nis, ris, rnis = load("IMR90_N1ICD_NIS"), load("IMR90_N1ICD_RIS"), load("IMR90_N1ICD_RNIS")
cls = list(csv.DictReader(open("results/exploratory/ras_vs_ois/strict/ois_gene_classes.csv")))
prolif = set(open("results/exploratory/two_sided_2026-11-01/prolif_genes.txt").read().split())
groups = defaultdict(list)
for r in cls:
    groups[(r["OIS_direction"], r["arrest_class"])].append(r["gene"])
    if r["OIS_direction"] == "UP" and r["nfkb_brd4"] == "dependent":
        groups[("UP", "NF-kB/BRD4-dependent")].append(r["gene"])
groups[("DOWN", "Hallmark proliferation")] = [g for g in prolif]
def frac(d, genes, sign):
    gs = [g for g in genes if g in d]
    return sum(d[g][1] and d[g][0] * sign > 0 for g in gs) / max(1, len(gs)), len(gs)
print(f"{'class':38} {'n':>5} {'in RIS':>7} {'in NIS':>7} {'in RNIS':>8}   (fraction changed same direction, IMR90 Hoare 2016)")
for (d, c), genes in sorted(groups.items()):
    s = 1 if d == "UP" else -1
    fr, n = frac(ris, genes, s); fn, _ = frac(nis, genes, s); frn, _ = frac(rnis, genes, s)
    print(f"{d + ' ' + c:38} {n:5} {fr:7.2f} {fn:7.2f} {frn:8.2f}")
mk = ["MKI67", "CDKN1A", "CDKN2A", "DUSP6", "SPRY2", "ETV4", "HES1", "HEY1", "NOTCH3", "IL1B", "IL6", "CXCL8", "TGFB1", "MMP1", "MMP3"]
print("\nlog2FC  " + " ".join(f"{g:>6}" for g in mk))
for n, d in (("RIS", ris), ("NIS", nis), ("RNIS", rnis)):
    print(f"{n:6}  " + " ".join(f"{d[g][0]:6.1f}" if g in d else "    na" for g in mk))

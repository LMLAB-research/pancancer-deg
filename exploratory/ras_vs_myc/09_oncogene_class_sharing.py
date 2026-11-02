"""Which RAS-induced-senescence (RIS) gene classes does each oncogene arm share?
Fraction of each class changed in the same direction (house sig), plus arrest and pathway
read-outs per contrast. Run from the pancancer-deg repo root."""
import csv, statistics
from collections import defaultdict
C = "results/exploratory/ras_vs_ois/contrast_degs"
def load(n):
    return {r["gene"].upper(): (float(r["log2FoldChange"]), r["sig"] == "y") for r in csv.DictReader(open(f"{C}/{n}.csv"))}
arms = [("IMR90 RAS (Hoare)", "IMR90_N1ICD_RIS"), ("IMR90 NOTCH1-ICD", "IMR90_N1ICD_NIS"),
        ("IMR90 RAS (array)", "IMR90_array_RAS"), ("IMR90 AKT (array)", "IMR90_array_AKT"),
        ("BJ-hTERT HRAS d14", "BJhT_HRAS_d14"), ("BJ-hTERT myrAKT1 d6", "BJhT_AKT_d6"),
        ("RPE-1 p53wt MYC 120h", "RPE1_p53wt_Myc_dox120h"), ("RPE-1 p53-/- MYC 120h", "RPE1_p53ko_Myc_dox120h")]
D = {a: load(n) for a, n in arms}
cls = list(csv.DictReader(open("results/exploratory/ras_vs_ois/strict/ois_gene_classes.csv")))
groups = defaultdict(list)
for r in cls:
    if r["arrest_class"] in ("RAS-autonomous", "arrest-coupled", "intermediate"):
        groups[(r["OIS_direction"], r["arrest_class"])].append(r["gene"])
    if r["OIS_direction"] == "UP" and r["nfkb_brd4"] == "dependent":
        groups[("UP", "NF-kB/BRD4-dependent")].append(r["gene"])
erk = ["DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "FOSL1", "PHLDA1", "SPRED1", "SPRED2", "DUSP4"]
myc_t = ["NPM1", "NCL", "NOP56", "NOP16", "CAD", "ODC1", "SRM", "PA2G4", "CDK4", "LDHA", "POLR1B", "RRP9"]
mk = ["MKI67", "CDKN1A", "CDKN2A", "CCNA2"]
w = max(len(a) for a, _ in arms)
print(f"{'arm':{w}} " + " ".join(f"{g:>7}" for g in mk) + f" {'ERK med':>8} {'MYCt med':>8}")
for a, _ in arms:
    d = D[a]
    f = lambda g: f"{d[g][0]:7.1f}" if g in d else "     na"
    med = lambda gs: statistics.median(d[g][0] for g in gs if g in d)
    print(f"{a:{w}} " + " ".join(f(g) for g in mk) + f" {med(erk):8.2f} {med(myc_t):8.2f}")
print("\nFraction of each RIS class changed in the same direction (house sig):")
print(f"{'class':30}" + "".join(f"{a[:20]:>22}" for a, _ in arms))
for (dr, c), genes in sorted(groups.items()):
    s = 1 if dr == "UP" else -1
    row = f"{dr + ' ' + c:30}"
    for a, _ in arms:
        gs = [g for g in genes if g in D[a]]
        row += f"{sum(D[a][g][1] and D[a][g][0] * s > 0 for g in gs) / max(1, len(gs)):>16.2f} (n={len(gs):4})"[-22:].rjust(22)
    print(row)

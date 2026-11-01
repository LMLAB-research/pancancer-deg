"""Partition OIS DEGs by whether other arrest types share them, then test which
partition the TCGA tumour DEGs overlap. Exploratory; reads the local pipeline outputs."""
import csv, glob, os, math
from collections import defaultdict
SRC = "/Users/ro/Academic/10_Research/senescence_stress/work/published_reference/arrest_degs_final.csv"
TCGA = "/Users/ro/code/academic/pancancer-deg/results/tables/TCGA_deg"
OV = "/Users/ro/code/academic/pancancer-deg/results/exploratory/two_sided_2026-11-01/full_signatures/OIS"
name = {"Oncogene_induced_CS": "OIS", "Replicative_CS": "RS", "Stress_induced_CS": "SIPS",
        "Contact_inhibited_CQ": "CICQ", "Serum_starved_CQ": "SSCQ"}
LFC = math.log2(1.5)
dirn = defaultdict(dict)
for r in csv.DictReader(open(SRC)):
    s = name[r["group_1"]]; g = r["gene"].upper()
    p, l = float(r["padj"]), float(r["log2FoldChange"])
    dirn[s][g] = "UP" if p < .05 and l > LFC else "DOWN" if p < .05 and l < -LFC else "NS"

def fisher_or(a, b, c, d):
    # Haldane-corrected log2 OR and a two-sided normal-approximation z (exploratory)
    a, b, c, d = a + .5, b + .5, c + .5, d + .5
    lor = math.log(a * d / (b * c)); se = math.sqrt(1/a + 1/b + 1/c + 1/d)
    return lor / math.log(2), 2 * (1 - 0.5 * (1 + math.erf(abs(lor / se) / math.sqrt(2))))

for D in ["UP", "DOWN"]:
    part = defaultdict(set)
    for g, v in dirn["OIS"].items():
        if v != D: continue
        cs = any(dirn[s].get(g) == D for s in ("RS", "SIPS"))
        cq = any(dirn[s].get(g) == D for s in ("CICQ", "SSCQ"))
        k = "OIS-only" if not (cs or cq) else "shared CS only" if cs and not cq else "shared with CQ"
        part[k].add(g)
    print(f"\nOIS-{D} genes: " + ", ".join(f"{k} {len(v)}" for k, v in part.items()))
    print(f"  {'cancer':6} " + "  ".join(f"{k:>22}" for k in ["OIS-only", "shared CS only", "shared with CQ"]) +
          f"   (tumour-{D} overlap: log2OR vs rest of universe)")
    agg = defaultdict(list)
    for d in sorted(glob.glob(OV + "/TCGA-*")):
        c = os.path.basename(d)[5:]
        univ = {r[next(iter(r))].upper() for r in csv.DictReader(open(d + "/common_tested_genes.csv"))}
        f = glob.glob(f"{TCGA}/TCGA-{c}_DGE_Raw_Results_Tumor_vs_Normal_Adjusted.csv")[0]
        tum = set()
        for r in csv.DictReader(open(f)):
            try: p, l = float(r["padj"]), float(r["log2FoldChange"])
            except ValueError: continue
            if p < .05 and ((l > LFC) if D == "UP" else (l < -LFC)): tum.add(r["gene_name"].upper())
        tum &= univ
        cells = []
        for k in ["OIS-only", "shared CS only", "shared with CQ"]:
            s = part[k] & univ
            a = len(s & tum); b = len(s) - a; cc = len(tum) - a; dd = len(univ) - a - b - cc
            lor, p = fisher_or(a, b, cc, dd); agg[k].append(lor)
            cells.append(f"{lor:6.2f} ({a:3d}/{len(s):4d}){'*' if p < .05/51 else ' '}")
        print(f"  {c:6} " + "  ".join(f"{x:>22}" for x in cells))
    print("  median " + "  ".join(f"{sorted(v)[len(v)//2]:>22.2f}" for v in agg.values()))

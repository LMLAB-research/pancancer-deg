"""Classify meta-analysis OIS DEGs by whether their change survives when RAS/BRAF stays
on but the arrest is reduced or lost, then test which classes TCGA tumour DEGs overlap.
Exploratory; run from the pancancer-deg repo root."""
import csv, glob, math, os, statistics
import sys; sys.path.insert(0, "exploratory"); from paths import *  # repo-relative inputs (paths.py)
from collections import Counter, defaultdict
from scipy.stats import fisher_exact

CF = contrast_deg_files()
STRICT = os.environ.get("STRICT", "0") == "1"
# NO_BRAF=1 drops the TIG3 BRAF "bypass" arm (hTERT + SV40 small T + BRAF vs cycling cells without small T)
NO_BRAF = os.environ.get("NO_BRAF", "0") == "1"
OUT = "results/exploratory/ras_vs_ois" + ("/strict" if STRICT else "/lenient") + ("_noBRAF" if NO_BRAF else "")
META = ARREST_PUBLISHED
LFC_T = math.log2(1.5)

def load(name):
    return {r["gene"].upper(): (float(r["log2FoldChange"]), r["sig"] == "y")
            for r in csv.DictReader(opent(CF[name]))}
c = {n: load(n) for n in CF}

# Meta-analysis OIS DEGs (published reference) and other arrest types
meta = defaultdict(dict)
lab = {"Oncogene_induced_CS": "OIS", "Replicative_CS": "RS", "Stress_induced_CS": "SIPS",
       "Contact_inhibited_CQ": "CICQ", "Serum_starved_CQ": "SSCQ"}
for r in csv.DictReader(opent(META)):
    p, l = float(r["padj"]), float(r["log2FoldChange"])
    meta[lab[r["group_1"]]][r["gene"].upper()] = "UP" if p < .05 and l > LFC_T else "DOWN" if p < .05 and l < -LFC_T else "NS"

def ratio(g, num, den):
    """Fraction of the reference change retained; None if gene missing or reference change small."""
    if g not in c[num] or g not in c[den]: return None
    d = c[den][g][0]
    if abs(d) < LFC_T: return None
    return c[num][g][0] / d

def retained_escape(g, sign):
    if g not in c["BJ_RAS_escaped_clones"]: return None
    l, s = c["BJ_RAS_escaped_clones"][g]
    return 1.0 if (s and l * sign > 0) else (0.5 if l * sign > LFC_T / 2 else 0.0)

os.makedirs(OUT, exist_ok=True)
rows = []
for g, d in meta["OIS"].items():
    if d == "NS": continue
    sign = 1 if d == "UP" else -1
    r = [ratio(g, "IMR90_RAS_shp53Rb", "IMR90_RAS_shRen"),
         None if NO_BRAF else ratio(g, "TIG3_BRAF_bypass_d10", "TIG3_BRAF_OIS_d10"),
         retained_escape(g, sign)]
    ev = [x for x in r if x is not None]
    keep = sum(x >= 0.5 for x in ev); lost = sum(x <= 0.25 for x in ev)
    if STRICT:
        # anchor on the fully proliferating escaped clones, plus >=1 partial bypass agreeing
        esc = r[2]; part = [x for x in r[:2] if x is not None]
        if esc is None or not part: cls = "unassessable"
        elif esc == 1.0 and any(x >= 0.5 for x in part): cls = "RAS-autonomous"
        elif esc == 0.0 and any(x <= 0.25 for x in part): cls = "arrest-coupled"
        else: cls = "intermediate"
    elif len(ev) < 2: cls = "unassessable"
    elif keep >= 2: cls = "RAS-autonomous"
    elif lost >= 2: cls = "arrest-coupled"
    else: cls = "intermediate"
    # NF-kB / BRD4 dependence of the RAS-senescence change (arrest maintained in these arms)
    dep = [ratio(g, "IMR90_RAS_shp65", "IMR90_RAS_shRen"), ratio(g, "IMR90_RAS_shBrd4", "IMR90_RAS_shRen"),
           ratio(g, "IMR90_RAS_JQ1", "IMR90_RAS_DMSO")]
    dep = [x for x in dep if x is not None]
    nfkb_brd4 = "dependent" if len(dep) >= 2 and sum(x <= 0.5 for x in dep) >= 2 else "independent" if dep else "NA"
    ris = c["IMR90_RIS"].get(g); ddis = c["IMR90_DDIS"].get(g)
    ris_spec = ("RIS-specific" if ris and ris[1] and ris[0] * sign > 0 and not (ddis and ddis[1] and ddis[0] * sign > 0)
                else "shared-with-DDIS" if ris and ddis and ris[1] and ddis[1] and ris[0] * sign > 0 and ddis[0] * sign > 0
                else "other")
    d4 = ratio(g, "IMR90hT_RAS_d4_shScr", "IMR90hT_RAS_d10_shScr")
    timing = "NA" if d4 is None else ("early" if d4 >= 0.5 else "late")
    other_cs = any(meta[s].get(g) == d for s in ("RS", "SIPS"))
    cq = any(meta[s].get(g) == d for s in ("CICQ", "SSCQ"))
    rows.append(dict(gene=g, OIS_direction=d, retention_shp53Rb=r[0], retention_BRAF_bypass=r[1],
                     retention_escape=r[2], arrest_class=cls, nfkb_brd4=nfkb_brd4, ris_vs_ddis=ris_spec,
                     timing=timing, shared_RS_SIPS=other_cs, shared_quiescence=cq))
with open(f"{OUT}/ois_gene_classes.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)

print("OIS DEG classes (meta-analysis OIS DEGs):")
for d in ("UP", "DOWN"):
    sub = [r for r in rows if r["OIS_direction"] == d]
    print(f"  OIS-{d} n={len(sub)}: ", dict(Counter(r["arrest_class"] for r in sub)))
    for k in ("nfkb_brd4", "ris_vs_ddis", "timing"):
        print(f"      {k}: ", dict(Counter(r[k] for r in sub)))
    for cls in ("RAS-autonomous", "arrest-coupled"):
        ss = [r for r in sub if r["arrest_class"] == cls]
        print(f"      {cls}: shared RS/SIPS {sum(r['shared_RS_SIPS'] for r in ss)}, shared quiescence {sum(r['shared_quiescence'] for r in ss)}")

# Hallmark proliferation genes should be arrest-coupled (sanity)
prolif = set(open(PROLIF).read().split())
pc = Counter(r["arrest_class"] for r in rows if r["gene"] in prolif and r["OIS_direction"] == "DOWN")
print("\nSanity: OIS-DOWN proliferation genes by class:", dict(pc))

# Project TCGA tumour DEGs onto the classes
groups = defaultdict(set)
for r in rows:
    groups[(r["OIS_direction"], r["arrest_class"])].add(r["gene"])
    if r["OIS_direction"] == "UP":
        groups[("UP", "nfkb_brd4-" + r["nfkb_brd4"])].add(r["gene"])
        groups[("UP", r["ris_vs_ddis"])].add(r["gene"])
res = []
for f in tcga_deg_files():
    ca = os.path.basename(f).split("_")[0][5:]
    tum = {"UP": set(), "DOWN": set()}; univ = set()
    for r in csv.DictReader(opent(f)):
        if r.get("gene_type") != "protein_coding": continue
        try: p, l = float(r["padj"]), float(r["log2FoldChange"])
        except ValueError: continue
        g = r["gene_name"].upper(); univ.add(g)
        if p < .05 and l > LFC_T: tum["UP"].add(g)
        elif p < .05 and l < -LFC_T: tum["DOWN"].add(g)
    univ &= set(meta["OIS"])
    for (d, cls), gs in groups.items():
        s = gs & univ; t = tum[d] & univ
        a = len(s & t); b = len(s) - a; cc = len(t) - a; dd = len(univ) - a - b - cc
        orr, p = fisher_exact([[a, b], [cc, dd]])
        res.append(dict(cancer=ca, ois_dir=d, cls=cls, n_class=len(s), overlap=a, odds_ratio=orr, p=p))
with open(f"{OUT}/tcga_overlap_by_class.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=res[0].keys()); w.writeheader(); w.writerows(res)

print("\nTumour overlap (same direction as OIS) by class: median log2 OR, n cancers p<0.05/enriched/depleted")
for (d, cls) in sorted(groups):
    rr = [r for r in res if r["ois_dir"] == d and r["cls"] == cls]
    lo = [math.log2(r["odds_ratio"]) if r["odds_ratio"] > 0 else -5 for r in rr]
    enr = sum(r["p"] < .05 and r["odds_ratio"] > 1 for r in rr); dep = sum(r["p"] < .05 and r["odds_ratio"] < 1 for r in rr)
    print(f"  OIS-{d:4} {cls:28} n={rr[0]['n_class']:5}  median log2OR {statistics.median(lo):5.2f}  enriched {enr:2}/17 depleted {dep:2}/17")

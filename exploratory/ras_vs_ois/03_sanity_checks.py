"""Check each contrast: proliferation genes (arrest on/off), ERK targets (RAS on),
and key arrest/SASP genes. Run from the pancancer-deg repo root."""
import csv, glob, os, statistics
import sys; sys.path.insert(0, "exploratory"); from paths import *  # repo-relative inputs (paths.py)
prolif = set(open(PROLIF).read().split())
erk = {"DUSP6", "SPRY2", "SPRY4", "ETV4", "ETV5", "FOSL1", "PHLDA1", "IER3", "SPRED1", "SPRED2", "DUSP4", "EGR1"}
show = ["CDKN2A", "CDKN1A", "TP53", "RB1", "MKI67", "IL1A", "IL1B", "IL6", "CXCL8", "MMP1", "MMP3", "NOTCH1"]
print(f"{'contrast':24} {'prolif med':>10} {'%prolif dn':>10} {'ERK med':>8}  " + " ".join(f"{g:>6}" for g in show))
for f in sorted(contrast_deg_files().values(), key=os.path.basename):
    rows = {r["gene"]: r for r in csv.DictReader(opent(f))}
    lfc = lambda g: float(rows[g]["log2FoldChange"]) if g in rows else float("nan")
    pl = [lfc(g) for g in prolif if g in rows]
    pdn = sum(rows[g]["sig"] == "y" and lfc(g) < 0 for g in prolif if g in rows) / max(1, len(pl))
    el = [lfc(g) for g in erk if g in rows]
    print(f"{os.path.basename(f)[:-4]:24} {statistics.median(pl):10.2f} {100*pdn:9.0f}% {statistics.median(el):8.2f}  "
          + " ".join(f"{lfc(g):6.1f}" for g in show))

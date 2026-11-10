"""Paths for the exploratory Python scripts (2026-11-09); mirrors paths.R. Run scripts from the pancancer-deg repo root.
Inputs from other projects are read from exploratory/inputs/ (originals in inputs/MANIFEST.tsv). Regenerable
intermediates (TCGA DEG tables, contrast DEGs) are read from results/ when present, else from exploratory/inputs/."""
import glob, gzip, os

EXP = "exploratory"
INP = os.path.join(EXP, "inputs")
if not os.path.isdir(INP):
    raise SystemExit(f"Run from the pancancer-deg repo root (exploratory/inputs not found from {os.getcwd()})")


def inp(*p):
    return os.path.join(INP, *p)


def opent(f):
    """Open a text file, transparently decompressing .gz."""
    return gzip.open(f, "rt") if f.endswith(".gz") else open(f)


ARREST_PUBLISHED = inp("senescence_stress", "arrest_degs_final_published.csv.gz")
PROLIF = inp("other", "prolif_genes.txt")
PURITY_FILE = inp("other", "TCGA_mastercalls.abs_tables_JSedit.fixed.txt.gz")


def tcga_deg_files():
    f = sorted(glob.glob("results/tables/TCGA_deg/TCGA-*_DGE_Raw_Results_Tumor_vs_Normal_Adjusted.csv"))
    return f or sorted(glob.glob(inp("tcga_degs", "TCGA-*_DGE_Raw_Results_Tumor_vs_Normal_Adjusted.csv.gz")))


def contrast_deg_files():
    """{contrast name: path}, results/ if regenerated there, else the repo copies."""
    f = glob.glob("results/exploratory/ras_vs_ois/contrast_degs/*.csv")
    if f:
        return {os.path.basename(x)[:-4]: x for x in f}
    return {os.path.basename(x)[:-7]: x for x in glob.glob(inp("contrast_degs", "*.csv.gz"))}


def contrast_deg(name):
    return contrast_deg_files()[name]


def overlap_results():
    """Two-sided overlap pipeline results: the 2026-11-01 run if present, else the output of run_both.sh."""
    for d in ("results/exploratory/two_sided_2026-11-01", os.environ.get("RUN_BOTH_OUT", "results/exploratory/two_sided_rerun")):
        if os.path.isdir(d):
            return d
    raise SystemExit("no overlap results: run bash exploratory/run_both.sh first")

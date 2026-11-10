"""Per-sample RAS-pathway hotspot and MYC amplification calls for TCGA PanCancer Atlas
studies (cBioPortal API). Exploratory; run from the pancancer-deg repo root.
RAS-pathway: KRAS/NRAS/HRAS missense at codons 12/13/61, BRAF V600.
MYC: GISTIC discrete copy number = 2 (high-level amplification), among CNA-profiled samples."""
import csv, json, os, urllib.request, urllib.error

API = "https://www.cbioportal.org/api"
OUT = "results/exploratory/ras_vs_myc"
os.makedirs(OUT, exist_ok=True)
GENES = {"KRAS": 3845, "NRAS": 4893, "HRAS": 3265, "BRAF": 673}
MYC = 4609

def _retry(fn, what, tries=6):
    import time
    for i in range(tries):
        try:
            return fn()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise
            err = e
        except Exception as e:  # truncated responses / proxy resets
            err = e
        time.sleep(2 * (i + 1))
    raise RuntimeError(f"{what} failed after {tries} tries: {err!r}")

# curl rather than urllib: urllib gets truncated chunked responses through some proxies.
import subprocess

def _curl(path, body=None):
    cmd = ["curl", "-sS", "--fail-with-body", "--retry", "3", "-H", "Accept: application/json"]
    if body is not None:
        cmd += ["-H", "Content-Type: application/json", "-X", "POST", "--data", json.dumps(body)]
    p = subprocess.run(cmd + [API + path], capture_output=True, text=True)
    if p.returncode == 22 and '"status":404' in p.stdout:
        raise urllib.error.HTTPError(API + path, 404, "Not Found", None, None)
    if p.returncode != 0:
        raise RuntimeError(f"curl rc={p.returncode}: {p.stderr.strip()} {p.stdout[:200]}")
    return json.loads(p.stdout)

def get(path):
    return _retry(lambda: _curl(path), path)

def post(path, body):
    return _retry(lambda: _curl(path, body), path)

def ras_hotspot(gene, pchange):
    p = (pchange or "").upper()
    if gene == "BRAF":
        return p.startswith("V600")
    return any(p.startswith(x) for x in ("G12", "G13", "Q61"))

studies = sorted(s["studyId"] for s in get("/studies?keyword=pan_can_atlas_2018&projection=SUMMARY"))
rows = []
for st in studies:
    try:
        seq = set(get(f"/sample-lists/{st}_sequenced/sample-ids"))
    except urllib.error.HTTPError:  # 404: study has no sequenced list
        seq = set()
    try:
        cna = set(get(f"/sample-lists/{st}_cna/sample-ids"))
    except urllib.error.HTTPError:
        cna = set()
    muts = post(f"/molecular-profiles/{st}_mutations/mutations/fetch",
                {"entrezGeneIds": list(GENES.values()), "sampleListId": f"{st}_all"}) if seq else []
    inv = {v: k for k, v in GENES.items()}
    ras = {}
    for m in muts:
        g = inv[m["entrezGeneId"]]
        if ras_hotspot(g, m.get("proteinChange")):
            ras.setdefault(m["sampleId"], set()).add(f"{g}_{m['proteinChange']}")
    amp = set()
    if cna:
        for c in post(f"/molecular-profiles/{st}_gistic/discrete-copy-number/fetch?discreteCopyNumberEventType=AMP",
                      {"entrezGeneIds": [MYC], "sampleListId": f"{st}_all"}):
            if c["alteration"] == 2:
                amp.add(c["sampleId"])
    for s in sorted(seq | cna):
        rows.append(dict(study=st.split("_")[0].upper(), sample=s, sequenced=s in seq, cna_profiled=s in cna,
                         ras_hotspot=";".join(sorted(ras.get(s, []))), ras_mut=s in ras,
                         myc_amp=(s in amp) if s in cna else ""))
    n_both = sum(1 for s in seq & cna)
    print(f"{st.split('_')[0].upper():9} seq={len(seq):4} cna={len(cna):4} RAS-hotspot={len(ras):4} MYC-amp={len(amp):4} both={len(set(ras) & amp):3}")
with open(f"{OUT}/driver_status.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)

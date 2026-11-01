import csv, statistics
T = "/tmp/claude-501/two"
sigs = ["OIS", "RS", "SIPS", "CICQ", "SSCQ"]
classes = ["UP-UP", "UP-DOWN", "DOWN-UP", "DOWN-DOWN"]
def load(run, sig):
    return {(r["cancer_type"], r["overlap_class"]): r for r in
            csv.DictReader(open(f"{T}/results_{run}/heatmaps/TCGA_{sig}_heatmap_data.csv"))}
print(f"{'sig':5} {'class':10} | {'full: enr/dep/ns':>16} {'med log2OR':>10} | {'no-prolif: enr/dep/ns':>21} {'med log2OR':>10}")
for s in sigs:
    m, p = load("main", s), load("sens", s)
    for c in classes:
        out = []
        for d in (m, p):
            rr = [v for k, v in d.items() if k[1] == c]
            enr = sum(r["significant"] == "TRUE" and float(r["log2_OR"]) > 0 for r in rr)
            dep = sum(r["significant"] == "TRUE" and float(r["log2_OR"]) < 0 for r in rr)
            med = statistics.median(float(r["log2_OR"]) for r in rr)
            out.append((f"{enr}/{dep}/{len(rr)-enr-dep}", med))
        print(f"{s:5} {c:10} | {out[0][0]:>16} {out[0][1]:>10.2f} | {out[1][0]:>21} {out[1][1]:>10.2f}")
# per-cancer outliers for OIS and others: cancers that break the pattern
print("\nDepleted cells (full, two-sided):")
for s in sigs:
    d = load("main", s)
    dep = sorted(f"{k[0]}:{k[1]}({float(v['log2_OR']):.2f})" for k, v in d.items()
                 if v["significant"] == "TRUE" and float(v["log2_OR"]) < 0)
    print(f"  {s}: {len(dep)} -> " + ", ".join(dep))

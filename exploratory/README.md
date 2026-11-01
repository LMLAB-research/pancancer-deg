# Exploratory: two-sided overlaps (local only, not for GitHub)

- `check_fisher.R`: fisher.test("greater") reproduces GeneOverlap p-values; two-sided detects depletion.
- `run_both.sh`: two-sided overlap pipeline on the full signatures and on signatures with
  Hallmark proliferation genes removed (E2F, G2M, mitotic spindle, MYC v1/v2).
- `summarise.py`: enriched/depleted counts per signature x overlap class, both runs.
- `explore.R`: recurrent overlap genes (>=12/17 cancers) and Hallmark content per class
  (reads the one-sided pipeline outputs; gene lists are the same either way).

Outputs: `results/exploratory/two_sided_2026-11-01/` in the main checkout (gitignored).

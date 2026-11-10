# pancancer-deg

Pan-cancer TCGA tumour-vs-normal DEGs (`analyses/TCGA_deg/01–04`, run from the repo
root) and their overlap with senescence/quiescence signatures
(`TCGA_senescence_overlapping_pipeline/`). See `README.md` for setup and layout.

## Running inside the Claude Code sandbox

- `.Rprofile` forces basic proxy auth (curl handle default + `httr::set_config`) when
  `https_proxy` is the sandbox's credentialed `127.0.0.1` proxy; otherwise it changes
  nothing. Without it every curl/httr/httr2 request (all TCGAbiolinks GDC calls) fails
  with "Proxy CONNECT aborted". Don't move this into a home-folder R startup file.
- `renv::restore()` / `setup/setup.sh` fail in the sandbox ("couldn't find open port":
  renv's installer needs a local socket). Install missing packages into the project
  library with `BiocManager::install(..., lib = renv::paths$library())` instead.

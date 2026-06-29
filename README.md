# pancancer-deg

pancancer-deg contains analysis code for the SenescenceCancer pan-cancer differential expression project. The current codebase focuses on TCGA differential expression analysis, with the repository layout prepared for additional analyses later in the project.

## R setup

This project uses `renv` for reproducible R package management.

```bash
bash setup/setup.sh
```

The setup script installs `renv` if needed, restores packages from `renv.lock`,
and runs a package availability check. To run only the check:

```bash
RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript --vanilla setup/check_setup.R
```

## Project Structure

```text
pancancer-deg/
|-- analyses/
|   `-- TCGA_deg/
|-- data/
|-- results/
|   |-- figures/
|   |   `-- TCGA_deg/
|   `-- tables/
|       `-- TCGA_deg/
|-- setup/
|   |-- check_setup.R
|   |-- install_dependencies.R
|   `-- setup.sh
|-- README.md
|-- renv.lock
`-- pancancer_deg.Rproj
```

The `analyses/` directory contains analysis scripts grouped by analysis module. The current module is `TCGA_deg/`; future modules can be added as sibling directories under `analyses/`.

The `results/` directory is organized by output type first, then by analysis module. Generated figures for the TCGA differential expression analysis should go under `results/figures/TCGA_deg/`, and generated tables should go under `results/tables/TCGA_deg/`.

The `data/` directory is intended for local input data, downloaded resources, and intermediate analysis objects. Large or generated data files should generally stay out of Git.

## TCGA_deg Analysis

Placeholder: this section will summarize the purpose, inputs, major processing steps, and outputs of the `TCGA_deg` analysis once the scripts are cleaned and finalized.

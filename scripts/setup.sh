#!/usr/bin/env bash
set -euo pipefail

echo "[1/4] Checking R installation..."
if ! command -v Rscript >/dev/null 2>&1; then
  echo "Error: Rscript not found. Please install R and try again."
  exit 1
fi

echo "[2/4] Installing renv if missing..."
RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript --vanilla -e "if (!requireNamespace('renv', quietly = TRUE)) install.packages('renv', repos = 'https://cloud.r-project.org')"

if [[ -f "renv.lock" ]]; then
  echo "[3/4] Restoring project packages from renv.lock..."
  RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript --vanilla -e "project <- normalizePath('.'); renv::restore(prompt = FALSE, project = project, library = renv::paths\$library(project = project))"
else
  echo "[3/4] renv.lock not found; installing project packages and creating lockfile..."
  RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript --vanilla scripts/install_dependencies.R
fi

echo "[4/4] Running setup checks..."
RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript --vanilla scripts/check_setup.R

echo "Setup complete."

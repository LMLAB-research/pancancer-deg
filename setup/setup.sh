#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_ROOT}"

echo "[1/5] Checking R installation..."
if ! command -v Rscript >/dev/null 2>&1; then
  echo "Error: Rscript not found. Please install R and try again."
  exit 1
fi

echo "[2/5] Installing renv if missing..."
RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript --vanilla -e "if (!requireNamespace('renv', quietly = TRUE)) install.packages('renv', repos = 'https://cloud.r-project.org')"

if [[ -f "renv.lock" ]]; then
  echo "[3/5] Restoring project packages from renv.lock..."
  RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript --vanilla -e "project <- normalizePath('.'); renv::restore(prompt = FALSE, project = project, library = renv::paths\$library(project = project))"
else
  echo "[3/5] renv.lock not found; skipping restore."
fi

echo "[4/5] Ensuring project packages are installed..."
RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript --vanilla setup/install_dependencies.R

echo "[5/5] Running setup checks..."
RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript --vanilla setup/check_setup.R

echo "Setup complete."

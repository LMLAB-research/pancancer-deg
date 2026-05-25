# pancancer-deg

## R setup

This project uses `renv` for reproducible R package management.

```bash
bash scripts/setup.sh
```

The setup script installs `renv` if needed, restores packages from `renv.lock`,
and runs a package availability check. To run only the check:

```bash
RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript --vanilla scripts/check_setup.R
```

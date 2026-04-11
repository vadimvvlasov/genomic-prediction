# AGENTS.md

## What this repo is
- R package `gp` for genomic prediction and cross-validation across/within populations.
- Main entrypoint is `gp(args)` in `R/main.R`; CLI wrapper is `inst/exec_Rscript/gp.R`.
- Core module flow is `io.R -> distances.R -> cross_validation.R -> models.R -> metrics.R` (orchestrated by `R/main.R`).

## Verified dev commands (use these)
- Enter env (Nix): `nix-shell --run bash --pure`
- Or create Conda env: `conda env create -f conda.yml`
- Run tests: `Rscript -e 'devtools::test()'`
- Run full package checks (CI-equivalent): `Rscript -e 'devtools::check()'`
- CI uses Ubuntu + R 4.3.1 and runs `devtools::check()` (`.github/workflows/r.yml`).

## Architecture and execution facts that matter
- `gp(args)` returns a path to an output `.Rds` file, not an in-memory object (`R/main.R`).
- Output list keys are consumed by the Shiny app in `inst/plot_gs_gp/app.R` (`TRAIT_NAME`, `POPULATION`, `METRICS_*`, `YPRED_*`, `GENOMIC_PREDICTIONS`, `ADDITIVE_GENETIC_EFFECTS`).
- Missing phenotypes are predicted only after model selection from within-pop CV (best mean `corr`) in `R/main.R`.
- Cross-validation uses `fn_cross_validation_preparation()` + `fn_cv_1()`; parallelism is CPU process-based via `parallel::mclapply`, not GPU (`R/cross_validation.R`).

## Data/CLI contract (easy to miss)
- Phenotype file must have: col1 `id`, col2 `pop`, col3+ traits (one trait per column). `inst/exec_Rscript/2-gp_slurm_job.sh` slices traits by column index.
- Genotype input for batch scripts is expected as `.Rds` matrix (`inst/exec_Rscript/2-gp_slurm_job.sh`); `gp.R` CLI itself supports broader formats via package loaders.
- In SLURM array mode, across-pop CV is intentionally run only for one pop index to avoid duplicate across-pop runs (`BOOL_ACROSS` logic in `inst/exec_Rscript/2-gp_slurm_job.sh`).

## Batch script quirks
- `inst/exec_Rscript/config.txt` is positional (line-number based parsing in `0-submit.sh`); reordering lines breaks substitution.
- `0-submit.sh` generates run-specific scripts (`*RAND*.sh`) and rewrites defaults via `sed`.
- `inst/exec_Rscript/` writes many temporary outputs (`bglr*`, `GENOMIC_PREDICTIONS_OUTPUT-*`, `output/`); these are intentionally gitignored.

## Model/runtime constraints
- CPU-oriented stack: `glmnet`, `BGLR`, `sommer`; no CUDA path in repo.
- Bayesian models use defaults `nIter=12000`, `burnIn=2000` unless overridden in code paths (`R/cross_validation.R`, `R/main.R`), so they are the main runtime cost.
- Parallel thread count is bounded by memory estimate (`fn_estimate_memory_footprint` in `R/io.R`), and code can auto-fallback to non-parallel.

## Editing conventions for this repo
- `NAMESPACE` is roxygen-generated; do not hand-edit (`NAMESPACE` header says generated).
- Package build excludes `inst/` via `.Rbuildignore`; demo/SLURM helpers are repo tooling, not shipped package artifacts.

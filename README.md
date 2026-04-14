# gp — Genomic Prediction

Genomic prediction cross-validation using penalised, Bayesian and mixed linear models.

|**Build Status**|**License**|
|:--------------:|:---------:|
| <a href="https://github.com/jeffersonfparil/gp/actions"><img src="https://github.com/jeffersonfparil/gp/actions/workflows/r.yml/badge.svg"></a> | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) |

A streamlined interface to calculate the breeding values of entries in breeding trials.

> **This project contains a Python port (`gp_py/`) of the original R package.**
> The Python implementation provides the same genomic prediction functionality with
> modern tooling (uv, scikit-learn, Pydantic, Typer CLI).

---

## Python Version & Requirements

| Component | Version |
|-----------|---------|
| **Python** | **>= 3.11** (developed & tested on **3.13**) |
| Package manager | **[uv](https://github.com/astral-sh/uv)** |
| Core deps | numpy >= 1.26, pandas >= 2.2, scipy >= 1.12, scikit-learn >= 1.4 |
| Optional | xgboost >= 2.0, lightgbm >= 4.0, pymc >= 5.16, rpy2 >= 3.5 |

---

## Installation

```bash
# Clone the repository
git clone https://github.com/jeffersonfparil/gp.git
cd gp

# Install dependencies with uv
uv sync

# Optional: install extra dependencies
uv sync --extra io --extra ml --extra notebook --extra bayes
```

---

## Quick Start

### CLI

```bash
# Show help
uv run gp-py --help

# Run genomic prediction
uv run gp-py \
    --geno inst/exec_Rscript/input/test_geno.Rds \
    --pheno inst/exec_Rscript/input/test_pheno.tsv \
    --trait-col 3 \
    --n-folds 3 \
    --n-reps 2
```

### Jupyter Notebooks

```bash
uv run jupyter lab notebooks/
```

---

## Architecture

### Python (`gp_py/`)

```
gp_py/
├── cli.py          # Typer CLI interface
├── io.py           # Input/output: genotype (RDS), phenotype (TSV), filtering, merge
├── cv.py           # K-fold cross-validation (within-population)
├── models.py       # Genomic prediction models (ridge, lasso, Bayes A/B/C, gBLUP, ML)
├── metrics.py      # Accuracy metrics: corr, rmse, mae, r2, heritability (h²)
├── pipeline.py     # Full pipeline: CV → best model → predict missing
├── schema.py       # Pydantic dataclasses (MergedData, CVRun, etc.)
└── plots.py        # Visualization: predicted vs actual scatter plots
```

Pipeline flow: **io → cv → models → metrics → pipeline**

### R (original)

```shell
R/
└── main.R
    ├── io.R
    ├── distances.R
    └── cross_validation.R
        └── models.R
            └── metrics.R
```

1. main.R - main function
2. io.R - input, output, filtering, and simulation
3. distances.R - genetic relationship, distance matrices and their inverses
4. cross_validation.R - k-fold cross validation within and across populations, pairwise-cross-validation, and leave-one-population-out cross-validation
5. models.R - genomic prediction models with the consistent signatures
6. metrics.R - genomic prediction accuracy metrics

---

## Models

### Python Implementation

| # | Model | Type | Backend |
|---|-------|------|---------|
| 1 | **Ridge** (RR-BLUP/GBLUP) | Linear | sklearn |
| 2 | **Lasso** | Linear | sklearn |
| 3 | **Elastic Net** | Linear | sklearn |
| 4 | **Bayes A** | Bayesian | native / BGLR |
| 5 | **Bayes B** | Bayesian | native / BGLR |
| 6 | **Bayes C** | Bayesian | native / BGLR |
| 7 | **gBLUP** | Mixed Model | statsmodels MixedLM |
| 8 | **SVR** | ML | sklearn |
| 9 | **RandomForest** | ML | sklearn |
| 10 | **XGBoost** | ML | xgboost |
| 11 | **LightGBM** | ML | lightgbm |

### R Implementation (original)

1. [Ridge](https://en.wikipedia.org/wiki/Ridge_regression) (a.k.a. GBLUP): $Cost_{ridge} = \Sigma(y - X\beta)^2 + \lambda\Sigma\beta^2$, where $\hat{\beta} = {(X^TX + \lambda I)^{-1} X^Ty}$
2. [Lasso](https://en.wikipedia.org/wiki/Lasso_(statistics)): $Cost_{lasso} = \Sigma(y - X\beta)^2 + \lambda\Sigma|\beta|$
3. [Elastic net](https://en.wikipedia.org/wiki/Elastic_net_regularization): $Cost_{ridge} = \Sigma(y - X\beta)^2 + (1-\alpha)\lambda\Sigma\beta^2 + \alpha\lambda\Sigma|\beta|$
4. [Bayes A](https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf): scaled t-distributed effects
5. [Bayes B](https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf): scaled t-distributed effects with probability $\pi$; and zero effects with probability $1-\pi$, where $\pi \sim \beta(\theta_1, \theta_2)$.
6. [Bayes C](https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf): normally distributed effects ($N(0, \sigma^2_{\beta})$) with probability $\pi$; and zero effects with probability $1-\pi$, where $\pi \sim \beta(\theta_1, \theta_2)$.
7. [gBLUP](https://link.springer.com/protocol/10.1007/978-1-62703-447-0_13): genotype best linear unbiased prediction (gBLUP) using genomic relationship matrix to predict missing breeding values via Direct-Inversion Newton-Raphson or Average Information (via the [sommer R package](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4894563/)).

---

## Development

### Python Dev Commands

```bash
# Install dependencies
uv sync

# Run CLI
uv run gp-py --help

# Run tests
uv run pytest

# Lint
uv run ruff check gp_py/

# Type check
uv run mypy gp_py/

# Run notebooks
uv run jupyter lab notebooks/
```

### R Dev Commands

With Nix: `nix-shell --run bash --pure`.
With Conda: `conda env create -f conda.yml`.

```R
devtools::test()
```

Or check the entire library:

```R
devtools::check()
```

---

## Output

- **Python:** `gp()` returns a path to output `.parquet` file
- **R:** `gp(args)` returns a path to an output `.Rds` file

Output list keys are consumed by the Shiny app in `inst/plot_gs_gp/app.R` (`TRAIT_NAME`, `POPULATION`, `METRICS_*`, `YPRED_*`, `GENOMIC_PREDICTIONS`, `ADDITIVE_GENETIC_EFFECTS`).

---

## Documentation

```R
?gp::gp
```

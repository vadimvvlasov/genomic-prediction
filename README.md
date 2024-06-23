# gp

Genomic prediction cross-validation using penalised, Bayesian and mixed linear models

|**Build Status**|**License**|
|:--------------:|:---------:|
| <a href="https://github.com/jeffersonfparil/gp/actions"><img src="https://github.com/jeffersonfparil/gp/actions/workflows/r.yml/badge.svg"></a> | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) |

A streamlined interface to calculate the breeding values of entries in breeding trials.

## Installation

```R
devtools::install_github("jeffersonfparil/gp")
```

## Architecture

```shell
R/
└── main.R
    ├── io.R
    └── cross_validation.R
        └── models.R
            └── metrics.R
```

1. main.R - main function
2. io.R - input, output, filtering, and simulation
3. cross_validation.R - k-fold cross validation within and across populations, pairwise-cross-validation, and leave-one-population-out cross-validation
4. models.R - genomic prediction models with the consistent signatures
5. metrics.R - genomic prediction accuracy metrics

## Models

1. [Ridge](https://en.wikipedia.org/wiki/Ridge_regression) (a.k.a. GBLUP): $Cost_{ridge} = \Sigma(y - X\beta)^2 + \lambda\Sigma\beta^2$, where $\hat{\beta} = {(X^TX + \lambda I)^{-1} X^Ty}$
2. [Lasso](https://en.wikipedia.org/wiki/Lasso_(statistics)): $Cost_{lasso} = \Sigma(y - X\beta)^2 + \lambda\Sigma|\beta|$
3. [Elastic net](https://en.wikipedia.org/wiki/Elastic_net_regularization): $Cost_{ridge} = \Sigma(y - X\beta)^2 + (1-\alpha)\lambda\Sigma\beta^2 + \alpha\lambda\Sigma|\beta|$
4. [Bayes A](https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf): scaled t-distributed effects
5. [Bayes B](https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf): scaled t-distributed effects with probability $\pi$; and zero effects with probability $1-\pi$, where $\pi \sim \beta(\theta_1, \theta_2)$.
6. [Bayes C](https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf): normally distributed effects ($N(0, \sigma^2_{\beta})$) with probability $\pi$; and zero effects with probability $1-\pi$, where $\pi \sim \beta(\theta_1, \theta_2)$.
7. [gBLUP](https://link.springer.com/protocol/10.1007/978-1-62703-447-0_13): genotype best linear unbiased prediction (gBLUP) using genomic relationship matrix to predict missing breeding values via Direct-Inversion Newton-Raphson or Average Information (via the [sommer R package](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4894563/)).

## Documentation and links to function-specific dodumentation

```R
?gp::gp
```

## Unit tests

With Nix: `nix-shell --run bash --pure`.
With Conda: `conda env create -f conda.yml`.

```R
devtools::test()
```

Or per module as:

```R
devtools::load_all()
library(sommer)
library(testthat)
source("tests/testthat/test-helpers.R")
source("tests/testthat/test-univariate_gx1.R")
source("tests/testthat/test-univariate_gxe.R")
```

Or check the entire library:

```R
devtools::check()
```

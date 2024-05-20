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

## Models



## Documentation

Write documentation or comments into your code a much as you can especially on top of your function definitions.

```R
devtools::document()
```

### Unit tests

Value modularity and write tests for each function definition.

See the main tests function and each Rscript for individual unit tests:

```R
devtools::check()
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

Test new models in [`tests/testthat/test-helpers.R`](tests/testthat/test-helpers.R) and adjust parsing accordingly via: [`R/helpers.R::fn_henderson_vs_newtonraphson_fit()`](R/helpers.R).
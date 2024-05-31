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

See [`ARCHITECTURE.md`](./ARCHITECTURE.md) for details.

```shell
R/
└── main.R
    ├── load.R
    └── cross_validation.R
        └── models.R
            └── metrics.R

```

## Models

1. [Ridge](https://en.wikipedia.org/wiki/Ridge_regression) (a.k.a. GBLUP): $Cost_{ridge} = \Sigma(y - X\beta)^2 + \lambda\Sigma\beta^2$, where $\hat{\beta} = {(X^TX + \lambda I)^{-1} X^Ty}$
2. [Lasso](https://en.wikipedia.org/wiki/Lasso_(statistics)): $Cost_{lasso} = \Sigma(y - X\beta)^2 + \lambda\Sigma|\beta|$
3. [Elastic net](https://en.wikipedia.org/wiki/Elastic_net_regularization): $Cost_{ridge} = \Sigma(y - X\beta)^2 + (1-\alpha)\lambda\Sigma\beta^2 + \alpha\lambda\Sigma|\beta|$
4. [Bayes A](https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf): scaled t-distributed effects
5. [Bayes B](https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf): scaled t-distributed effects with probability $\pi$; and zero effects with probability $1-\pi$, where $\pi \sim \beta(\theta_1, \theta_2)$.
6. [Bayes C](https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf): normally distributed effects ($N(0, \sigma^2_{\beta})$) with probability $\pi$; and zero effects with probability $1-\pi$, where $\pi \sim \beta(\theta_1, \theta_2)$.
7. [gBLUP](https://link.springer.com/protocol/10.1007/978-1-62703-447-0_13): genotype best linear unbiased prediction (gBLUP) using genomic relationship matrix to predict missing breeding values via Direct-Inversion Newton-Raphson or Average Information (via the [sommer R package](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4894563/)).

The best models are identified through k-fold cross-fold validation (within populations), and leave-one-population-out cross-validation (across populations) using Pearson's correlation as the primary genomic prediction accuracy metric. If there are samples with known genotype data but missing phenotype data, their phenotypes will be predicted. The genomic prediction model used depends on the population the samples are part of. If the samples are part of a population with known genotype and phenotype data, the best model for that population is used, otherwise the best overall model is used.

## Execution function

Open the [`main function's`](./src/main.R) the documentation:

```shell
Rscript src/main.R -h
```

## Input data

1. Genotype data (assumes pre-filtered data; `-g; --genotype-file`)
    - [Variant call format (`*.vcf`)](https://www.internationalgenome.org/wiki/Analysis/vcf4.0/#:~:text=VCF%20is%20a%20text%20file,for%20each%20position%20or%20not.) which requires either the genotype call (**GT**) or allele depth (**AD**) field, where the AD field has precedence if both are present (i.e., allele frequencies rather than binary genotypes are assumed by default). See [tests/test.vcf](tests/test.vcf) for a very small example.
    - [R object file format (`*.rds`)](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/readRDS) which codes for a numeric matrix where each row correspond to a sample and each column to a locus represented by a single allele. See [tests/test.rds](tests/test.rds) for a very small example but note that the column or loci names does not need to follow any format, i.e. the allele names does not need to be included.
        + The genotype data can be coded as any numeric range of values, e.g. (0,1,2), (-1,0,1), and (0.00,0.25,0.50,0.75,1.00) or as biallelic characters, e.g. for diploids: "AA", "AB", "BB", and for tetraploids: "AAAA", "AAAB", "AABB", "ABBB", and "BBBB".. It is recommended that this data should be filtered and imputed beforehand.
        + The rows are expected to have names of the samples corresponding to the names in the phenotype file.
        + The columns are expected to contain the loci names but does need to follow a specific format: chromosome name and position separated by a tab character (`\t`) and an optional allele identifier, e.g. `chr-1\t12345\tallele_A`
    - [Allele frequency table format (`*.tsv`)]

    **TODO**

    - It is highly recommended to use genotype data without missing values. Missing data will be naively imputed using mean genotype value per allele/locus.
2. Phenotype data in text-delimited format (`-p; --phenotype-file`)
    - e.g. tab-delimited (`*.txt`) and comma-separated (`*.csv`)
    - may or may not have a header line
    - at least three columns are required:
        + names of the samples
        + population or grouping ID of each sample
        + numeric phenotype data (may have missing data coded as whatever is convenient; but cannot be completely missing)
    - Samples with known genotype and unknown phenotype data will have their phenotypes predicted. The best genomic prediction model (identified after k-fold cross-validation) for the population the samples belong to will be used. However, if the samples do not belong to populations with known genotype and phenotype data, then the best overall model will be used.
    - See [tests/test_pheno.csv](tests/test_pheno.csv) for a very small example.
3. Covariate matrix (`-c; --covariate-file`)
    - *.Rds* format
    - row names overlap with the row names (sample names) of the genotype and phenotype files [Default=NULL]
    - Each column will be standard normalised per column prior to cross-validations and genomic prediction.
    - For ridge, Lasso, and elastic net regression, the covariates will be placed before the genotype data and treated the same way as the genotype matrix.
    - For Bayes A, B, and C, the covariates will be treated as fixed effects with a flat prior.
4. Remove the minor allele from the genotype matrix? (`--retain-minus-one-alleles-per-locus`) [Default=TRUE]
5. Delimiter of the phenotype file (`--phenotype-delimiter`) [Default=',']
6. Does the phenotype file have a header line? (`--phenotype-header`) [Default=TRUE]
7. Which column contains the sample names? (`--phenotype-column-id`) [Default=1]
8. Which column contains the population or group names? (`--phenotype-column-pop`) [Default=2]
9. Which column/s contain the phenotype data? If more than 1 column is requested, then separate them with commas, e.g. "3,4,5". (`--phenotype-column-data`) [Default="3"]
10. How are the missing data in the phenotype file coded? (`--phenotype-missing-strings`) [Default=c("", "-", "NA", "na", "NaN", "missing", "MISSING")]
11. Genomic prediction models to use (`--models`) [Default=c("ridge","lasso","elastic_net","Bayes_A","Bayes_B","Bayes_C")]
12. Number of folds to perform in within population k-fold cross-validation (`--k-folds`) [Default=10]
13. Number of reshuffled replications per fold of within population k-fold cross-validation (`--n-reps`) [Default=3]
14. Number of computing threads to use in within and across population cross-validations (parallel computation across replications, folds, and models) (`--n-threads`) [Default=8]
15. Prefix of the filename of the output Rds file containing a list of lists (each corresponding to a trait requested by `--phenotype-column-data`) of metrics and predicted phenotype dataframes. The base name of the input genotype file will be used by default. (`-o; --output-file-prefix`) [Default='']"
16. Do you wish to print detailed progress reports of the workflow? (`--verbose`) [Default=FALSE]

## Output data

Rds file containing a list of lists each corresponding to a trait requested by `--phenotype-column-data`. Each list or trait consists of:

1. `TRAIT_NAME`: name of the trait
2. `SUMMARY`: table of the best models per population (Field names: *pop* (population names), *model* (top model), *corr* (mean correlation), *corr_sd* (standard deviation of the correlation))
3. `METRICS_WITHIN_POP`: table of genomic prediction accuracy metrics across replications, folds, and models per population
4. `YPRED_WITHIN_POP`: table of predicted and expected phenotypes across samples, replications and models per population.
5. `METRICS_ACROSS_POP`: table of genomic prediction accuracy metrics across 1 replication, 1 fold, and all models per validation population, i.e., via leave-one-population-out cross-validation **(present if there are at least 2 populations)**.
6. `YPRED_ACROSS_POP`: table of predicted and expected phenotypes across samples, 1 replication and all models per validation population, i.e., via leave-one-population-out cross-validation **(present if there are at least 2 populations)**.
7. `GENOMIC_PREDICTIONS`: table of predicted phenotypes of samples with known genotype and missing phenotypes **(present if there are samples with known genotype and unknown phenotype data)**. Note that if you have samples missing phenotypes which belong to a population without known phenotype data and wish to use the overall best model to predict their GEBVs then please do not use [`00_gs_slurm_job_wrapper.sh`](./00_gs_slurm_job_wrapper.sh) because it runs per trait per population. Instead, run [`src/main.R`](./src/main.R) separately using `--populations-to-include="all"`.


These dataframes have the following fields:

`SUMMARY` fields:
| pop | model | corr | corr_sd |
|:---:|:-----:|:----:|:-------:|
|:---:|:-----:|:----:|:-------:|

`METRICS_*_POP` fields:
| r | k | model | mbe | mae | rmse | r2 | corr | n_non_zero | duration_mins | pop (within) or validation_pop (across) |
|:-:|:-:|:-----:|:---:|:---:|:----:|:--:|:----:|:----------:|:-------------:|:---------------------------------------:|
|:-:|:-:|:-----:|:---:|:---:|:----:|:--:|:----:|:----------:|:-------------:|:---------------------------------------:|

`YPRED_*_POP` fields:
| entry | rep | model | y_pred | y_true |
|:-----:|:---:|:-----:|:------:|:------:|
|:-----:|:---:|:-----:|:------:|:------:|

`GENOMIC_PREDICTIONS` fields:
| pop | entry | y_pred | top_model | corr_from_kfoldcv |
|:---:|:-----:|:------:|:---------:|:-----------------:|
|:---:|:-----:|:------:|:---------:|:-----------------:|

These fields are defined as:

- **pop** or **validation_pop**: population name which refer to the name of the current population for within population k-fold cross-validations (**pop**) or the names of the validation population for leave-one-population-out across population cross-validation (**validation_pop**)
- **model**: genomic prediction model, i.e.,
    + ridge
    + lasso
    + elastic_net
    + Bayes_A
    + Bayes_B
    + Bayes_C
- **corr**: Pearson's product moment correlation = $E((y_{true}-\mu_{y_{true}}) (y_{predicted}-\mu_{y_{predicted}})) / \sqrt{ E((y_{true}-\mu_{y_{true}})^2) E((y_{predicted}-\mu_{y_{predicted}})^2) }$ (averaged across folds and replications)
- **corr_sd**: Standard deviation of the correlation coefficient across folds and replications
- **r**: replication
- **k**: fold
- **mbe**: mean bias error = $E(y_{true} - y_{predicted})$
- **mae**: mean absolute error = $E(| y_{true} - y_{predicted} |)$
- **rmse**: root mean square error = $\sqrt{E((y_{true} - y_{predicted})^2)}$
- **r2**: coefficient of determination = $1.00 - (E((y_{true} - y_{predicted})^2) / E((y_{true} - E(y_{true}))^2))$
- **corr**: Pearson's product moment correlation = $Cov(y_{true}, y_{predicted}) / (\sigma(y_{true}) \sigma(y_{predicted}))$
- **n_non_zero**: number of non-zero estimates of allele/locus effects
- **duration_mins**: run time in minutes
-  **entry**: names of the samples which may be individual genotypes or pool/groups of individuals
- **rep**: same as **r**, i.e., replication
- **y_pred**: predicted phenotype values
- **y_true**: true phenotype values
- **top_model**: model used in predicting unknown phenotypes, i.e. best model within population or overall best model
- **corr_from_kfoldcv**: prediction accuracy of the model used in predicting unknown phenotypes as shown by within population k-fold cross-validation, e.g. 0.42(±0.1)|pop1 which refers to 0.42 ± 0.1 correlation between predicted and expected phenotypes as assessed within populutaion pop1

For each model and fold of cross-validation (k-fold and leave-one-population-out) across replications and traits, an intermediate output file is generated containing the replication, fold, model name, and genomic prediction metrics with 2 lines: the header and the values. The filenames of these intermediate files are structured as `${TRAIT}-pop_${POPULATION}-metrics-rep_${REP}-fold_${FOLD}-model_${MODEL}-${DATE_TIME}.${RANODM_NUMBER}.csv`. Additionally, the Bayesian models generate intermediate files containing parameter estimates, and residual variances with the following filename structures: 

- `bglr_${BAYES_MODEL}-${TRAIT}-pop_${POPULATION}--${DATE_TIME}-${RANDOM_NUMBER}varE.dat`
- `bglr_${BAYES_MODEL}-${TRAIT}-pop_${POPULATION}--${DATE_TIME}-${RANDOM_NUMBER}mu.dat`
- `bglr_Bayes_C-${TRAIT}-pop_${POPULATION}--${DATE_TIME}-${RANDOM_NUMBER}ETA_MRK_parBayesC.dat`

These intermediate files are removed iteratively with `unlink` after merging all the output into a single Rds file. However, if for some reason the run of `src/main.R` was halted before finishing, then these intermediate files may remain and will need to be cleaned up or can be used to **salvage the failed run**.




## Documentation

Write documentation or comments into your code a much as you can especially on top of your function definitions.

```R
devtools::document()
```

## Unit tests

On Nix: `nix-shell --run bash --pure`. See [shell.nix](./shell.nix).

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
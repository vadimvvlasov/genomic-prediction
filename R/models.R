# source("R/io.R")
# source("R/metrics.R")
suppressWarnings(suppressPackageStartupMessages(library(glmnet)))
suppressWarnings(suppressPackageStartupMessages(library(BGLR)))
suppressWarnings(suppressPackageStartupMessages(library(sommer)))
# suppressWarnings(suppressPackageStartupMessages(library(epistasisfeatures)))

### New models can be added by having the same input structure, i.e.:
###     1.) nxp matrix of numeric genotype data without missing values
###     2.) n-sized vector of numeric phenotype data without missing values
###     3.) vector of numeric indexes corresponding to the index of the training set in the genotype matrix and phenotype vector
###     4.) vector of numeric indexes corresponding to the index of the validation set in the genotype matrix and phenotype vector
###     5.) a list of other parameters required by the internal functions

#' Ordinary least squares model
#'
#' @param list_merged list of merged genotype matrix, and phenotype vector, as well as an optional covariate matrix
#'  - $G: numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'      Row names can be any string of characters which identify the sample or entry or pool names.
#'      Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'      the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'      subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'  - $list_pheno:
#'      + $y: named vector of numeric phenotype data
#'      + $pop: population or groupings corresponding to each element of y
#'      + $trait_name: name of the trait
#'  - $COVAR: numeric n samples x k covariates matrix with non-null row and column names.
#' @param vec_idx_training vector of numeric indexes referring to the training set
#' @param vec_idx_validation vector of numeric indexes referring to the validation set
#' @param other_params list of additional parameters, specifically $diag_inflate which refers to the value used for 
#'  diagonal inflation for singular matrices (Default=list(diag_inflate=1e-4))
#' @param verbose show ordinary least squares regression messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'      + $list_perf:
#'          - $mbe: mean bias error
#'          - $mae: mean absolute error
#'          - $rmse: root mean squared error
#'          - $r2: coefficient of determination
#'          - $corr: Pearson's product moment correlation
#'          - $power_t10: fraction of observed top 10 phenotype values correctly predicted
#'          - $power_b10: fraction of observed bottom 10 phenotype values correctly predicted
#'          - $var_additive: variance of predicted phenotype values (estimator of additive genetic variance)
#'          - $var_residual: variance of differnce between observed and predicted phenotype values (estimator of residual variance)
#'          - $h2: narrow-sense heritability estimate
#'      + $df_y_validation:
#'          - $id: names of the samples/entries/pools, 
#'          - $pop: population from which the sample/entry/pool belongs to
#'          - $y_true: observed phenotype values
#'          - $y_pred: predicted phenotype values
#'      + $vec_effects: named numeric vector of estimated effects, where the names correspond to the
#'          SNP/allele identity including chromosome/scaffold, position and optionally allele.
#'      + $n_non_zero: number of non-zero estimated effects (effects greater than machine epsilon ~2.2e-16)
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
#' n = nrow(list_merged$G)
#' vec_idx_training = sample(c(1:n), floor(n/2))
#' vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
#' list_ols = fn_ols(list_merged, vec_idx_training, vec_idx_validation, verbose=TRUE)
#' @export
fn_ols = function(list_merged, vec_idx_training, vec_idx_validation, other_params=list(diag_inflate=1e-4), verbose=FALSE) {
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # COVAR = G %*% t(G); # rownames(COVAR) = NULL
    # list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    # n = nrow(list_merged$G)
    # vec_idx_training = sample(c(1:n), floor(n/2))
    # vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
    # other_params=list(diag_inflate=1e-6)
    # verbose = TRUE
    ###################################################
    ### Input sanity check
    if (methods::is(list_merged, "gpError")) {
        error = chain(list_merged, 
            methods::new("gpError",
                code=400,
                message=paste0(
                    "Error in models::fn_ols(...). ",
                    "Input data (list_merged) is an error type."
                )))
        return(error)
    }
    if (is.logical(c(vec_idx_training, vec_idx_training)) |
        (sum(is.na(c(vec_idx_training, vec_idx_training))) > 0) |
        (min(c(vec_idx_training, vec_idx_training)) < 0) |
        (max(c(vec_idx_training, vec_idx_training)) > nrow(list_merged$G))) {
            error = methods::new("gpError",
                code=401,
                message=paste0(
                    "Error in models::fn_ols(...). ",
                    "Please make sure that the vector of indexes for the training and validation sets are not booleans. ",
                    "We require the indexes to be positive integers without missing values. ",
                    "Also, the maximum index (", max(c(vec_idx_training, vec_idx_training)), 
                    ") should be less than or equal to the total number of samples/entries/pools in the input data: (", 
                    nrow(list_merged$G), ")."))
            return(error)
    }
    X_training = list_merged$G[vec_idx_training, , drop=FALSE]
    y_training = list_merged$list_pheno$y[vec_idx_training]
    X_training = X_training[!is.na(y_training), , drop=FALSE] ### Remove missing phenotype data from the training set
    y_training = y_training[!is.na(y_training)] ### Remove missing phenotype data from the training set
    X_validation = list_merged$G[vec_idx_validation, , drop=FALSE]
    y_validation = list_merged$list_pheno$y[vec_idx_validation]
    ### Adding covariate to explanatory matrix
    if (!is.null(list_merged$COVAR)) {
        X_training = cbind(list_merged$COVAR[vec_idx_training, , drop=FALSE], X_training)
        X_validation = cbind(list_merged$COVAR[vec_idx_validation, , drop=FALSE], X_validation)
    }
    ### Solve the least squares equations
    n = nrow(X_training)
    p = ncol(X_training)
    if (n >= p) {
        ### Canonical least squares equations
        A = t(X_training) %*% X_training
        A_inv = tryCatch(solve(A), error=function(x) {
            diag(A) = diag(A) + other_params$diag_inflate
            x_inv = tryCatch(solve(A), error=function(x) {
                NA
            })
        })
        if (is.na(A_inv[1])) {
            error = methods::new("gpError",
                code=402,
                message=paste0(
                    "Error in models::fn_ols(...). ",
                    "Failed to compute the inverse of X'X. ",
                    "Consider increasing the other_params$diag_inflate: ", other_params$diag_inflate, 
                    " to force the calculation of the inverse of X'X."
                ))
            return(error)
        }
        b = A_inv %*% t(X_training) %*% y_training
    } else {
        ### Non-canonical least squares equations
        A = X_training %*% t(X_training)
        A_inv = tryCatch(solve(A), error=function(A) {
            B = A
            diag(B) = diag(B) + other_params$diag_inflate
            B_inv = tryCatch(solve(B), error=function(B) {
                return(NA)
            })
            return(B_inv)
        })
        if (is.na(A_inv[1])) {
            error = methods::new("gpError",
                code=403,
                message=paste0(
                    "Error in models::fn_ols(...). ",
                    "Failed to compute the inverse of XX'. ",
                    "Consider increasing the other_params$diag_inflate: ", other_params$diag_inflate, 
                    " to force the calculation of the inverse of XX'."
                ))
            return(error)
        }
        b_hat = t(X_training) %*% A_inv %*% y_training
    }
    ### Name the effects and count the number of non-zero effects
    b_hat = b_hat[,1]
    names(b_hat) = colnames(X_training)
    n_non_zero = sum(abs(b_hat) > .Machine$double.eps)
    if (verbose) {
        print("Allele effects distribution:")
        tryCatch(txtplot::txtdensity(b_hat[!is.na(b_hat) & !is.infinite(b_hat)]), error=function(e){print("Empty plot")})
        print(paste0("Number of non-zero effects: ", n_non_zero, " (", round(100*p/n_non_zero), "%)"))
    }
    ### Predict and assess prediction accuracy
    y_pred = X_validation %*% b_hat
    df_y_validation = merge(
        data.frame(id=names(y_validation), pop=list_merged$list_pheno$pop[vec_idx_validation], y_true=y_validation), 
        data.frame(id=rownames(y_pred), y_pred=y_pred),
        by="id")
    list_perf = fn_prediction_performance_metrics(y_true=df_y_validation$y_true, y_pred=df_y_validation$y_pred, verbose=verbose)
    ### Output
    return(list(
        list_perf=list_perf,
        df_y_validation=df_y_validation,
        vec_effects=b_hat,
        n_non_zero=n_non_zero))
}

#' Ridge regression model (a.k.a. GBLUP; regularisation via Tikhonov regularisation; alpha=0)
#'
#' @param list_merged list of merged genotype matrix, and phenotype vector, as well as an optional covariate matrix
#'  - $G: numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'      Row names can be any string of characters which identify the sample or entry or pool names.
#'      Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'      the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'      subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'  - $list_pheno:
#'      + $y: named vector of numeric phenotype data
#'      + $pop: population or groupings corresponding to each element of y
#'      + $trait_name: name of the trait
#'  - $COVAR: numeric n samples x k covariates matrix with non-null row and column names.
#' @param vec_idx_training vector of numeric indexes referring to the training set
#' @param vec_idx_validation vector of numeric indexes referring to the validation set
#' @param other_params list of additional parameters, specifically $n_folds which refers to the number of
#'  internal cross-fold validation to find the optimum lambda at which the deviance is minimised (Default=list(n_folds=10))
#' @param verbose show ridge regression messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'      + $list_perf:
#'          - $mbe: mean bias error
#'          - $mae: mean absolute error
#'          - $rmse: root mean squared error
#'          - $r2: coefficient of determination
#'          - $corr: Pearson's product moment correlation
#'          - $power_t10: fraction of observed top 10 phenotype values correctly predicted
#'          - $power_b10: fraction of observed bottom 10 phenotype values correctly predicted
#'          - $var_additive: variance of predicted phenotype values (estimator of additive genetic variance)
#'          - $var_residual: variance of differnce between observed and predicted phenotype values (estimator of residual variance)
#'          - $h2: narrow-sense heritability estimate
#'      + $df_y_validation:
#'          - $id: names of the samples/entries/pools, 
#'          - $pop: population from which the sample/entry/pool belongs to
#'          - $y_true: observed phenotype values
#'          - $y_pred: predicted phenotype values
#'      + $vec_effects: named numeric vector of estimated effects, where the names correspond to the
#'          SNP/allele identity including chromosome/scaffold, position and optionally allele.
#'      + $n_non_zero: number of non-zero estimated effects (effects greater than machine epsilon ~2.2e-16)
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
#' n = nrow(list_merged$G)
#' vec_idx_training = sample(c(1:n), floor(n/2))
#' vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
#' list_ols = fn_ols(list_merged, vec_idx_training, vec_idx_validation, verbose=TRUE)
#' @export
fn_ridge = function(list_merged, vec_idx_training, vec_idx_validation, other_params=list(n_folds=10), verbose=FALSE) {
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # COVAR = G %*% t(G); # rownames(COVAR) = NULL
    # list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    # n = nrow(list_merged$G)
    # vec_idx_training = sample(c(1:n), floor(n/2))
    # vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
    # other_params=list(n_folds=10)
    # verbose = TRUE
    ###################################################
    ### Input sanity check
    if (methods::is(list_merged, "gpError")) {
        error = chain(list_merged, 
            methods::new("gpError",
                code=404,
                message=paste0(
                    "Error in models::fn_ridge(...). ",
                    "Input data (list_merged) is an error type."
                )))
        return(error)
    }
    if (is.logical(c(vec_idx_training, vec_idx_training)) |
        (sum(is.na(c(vec_idx_training, vec_idx_training))) > 0) |
        (min(c(vec_idx_training, vec_idx_training)) < 0) |
        (max(c(vec_idx_training, vec_idx_training)) > nrow(list_merged$G))) {
            error = methods::new("gpError",
                code=405,
                message=paste0(
                    "Error in models::fn_ridge(...). ",
                    "Please make sure that the vector of indexes for the training and validation sets are not booleans. ",
                    "We require the indexes to be positive integers without missing values. ",
                    "Also, the maximum index (", max(c(vec_idx_training, vec_idx_training)), 
                    ") should be less than or equal to the total number of samples/entries/pools in the input data: (", 
                    nrow(list_merged$G), ")."))
            return(error)
    }
    X_training = list_merged$G[vec_idx_training, , drop=FALSE]
    y_training = list_merged$list_pheno$y[vec_idx_training]
    X_training = X_training[!is.na(y_training), , drop=FALSE] ### Remove missing phenotype data from the training set
    y_training = y_training[!is.na(y_training)] ### Remove missing phenotype data from the training set
    X_validation = list_merged$G[vec_idx_validation, , drop=FALSE]
    y_validation = list_merged$list_pheno$y[vec_idx_validation]
    ### Adding covariate to explanatory matrix
    if (!is.null(list_merged$COVAR)) {
        X_training = cbind(list_merged$COVAR[vec_idx_training, , drop=FALSE], X_training)
        X_validation = cbind(list_merged$COVAR[vec_idx_validation, , drop=FALSE], X_validation)
    }
    ### Solve via ridge regularisation
    sol = tryCatch(
        glmnet::cv.glmnet(x=X_training, y=y_training, alpha=0, nfolds=other_params$n_folds, parallel=FALSE),
        error = function(e) {NA})
    if (is.na(sol)) {
        return(list(
            list_perf=NA,
            df_y_validation=NA,
            vec_effects=NA,
            n_non_zero=NA))
    }
    ### Find the first lambda with the lowest squared error (deviance) while having non-zero SNP effects
    vec_idx_decreasing_deviance = order(sol$glmnet.fit$dev.ratio, decreasing=FALSE)
    idx_start = which(sol$lambda == sol$lambda.min)[1]
    idx_start_vec_idx_decreasing_deviance = which(vec_idx_decreasing_deviance == idx_start)
    vec_idx_decreasing_deviance = vec_idx_decreasing_deviance[idx_start_vec_idx_decreasing_deviance:length(vec_idx_decreasing_deviance)]
    vec_idx_decreasing_deviance = vec_idx_decreasing_deviance[!is.na(vec_idx_decreasing_deviance)]
    for (i in vec_idx_decreasing_deviance) {
        # i = idx_start
        b_hat = c(sol$glmnet.fit$a0[i], sol$glmnet.fit$beta[, i])
        names(b_hat) = c("intercept", colnames(X_training))
        y_pred = cbind(rep(1, length(vec_idx_validation)), X_validation) %*% b_hat
        n_non_zero = sum(abs(b_hat) >= .Machine$double.eps)
        ### Skip models with only intercept effects
        if (n_non_zero == 1) {
            next
        } else {
            break
        }
    }
    if (verbose) {
        p = ncol(X_training)
        print("Allele effects distribution:")
        tryCatch(txtplot::txtdensity(b_hat[!is.na(b_hat) & !is.infinite(b_hat)]), error=function(e){print("Empty plot")})
        print("Relative positions of allele effects across the genome:")
        tryCatch(txtplot::txtplot(b_hat[!is.na(b_hat) & !is.infinite(b_hat)]), error=function(e){print("Empty plot")})
        print(paste0("Number of non-zero effects: ", n_non_zero, " (", round(100*n_non_zero/p), "%)"))
    }
    ### Evaluate prediction performance
    df_y_validation = merge(
        data.frame(id=names(y_validation), pop=list_merged$list_pheno$pop[vec_idx_validation], y_true=y_validation), 
        data.frame(id=rownames(y_pred), y_pred=y_pred),
        by="id")
    list_perf = fn_prediction_performance_metrics(y_true=df_y_validation$y_true, y_pred=df_y_validation$y_pred, verbose=verbose)
    ### Output
    return(list(
        list_perf=list_perf,
        df_y_validation=df_y_validation,
        vec_effects=b_hat,
        n_non_zero=n_non_zero))
}

#' Lasso regression model (regularisation via least absolute shrinkage and selection operator; alpha=1)
#'
#' @param list_merged list of merged genotype matrix, and phenotype vector, as well as an optional covariate matrix
#'  - $G: numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'      Row names can be any string of characters which identify the sample or entry or pool names.
#'      Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'      the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'      subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'  - $list_pheno:
#'      + $y: named vector of numeric phenotype data
#'      + $pop: population or groupings corresponding to each element of y
#'      + $trait_name: name of the trait
#'  - $COVAR: numeric n samples x k covariates matrix with non-null row and column names.
#' @param vec_idx_training vector of numeric indexes referring to the training set
#' @param vec_idx_validation vector of numeric indexes referring to the validation set
#' @param other_params list of additional parameters, specifically $n_folds which refers to the number of
#'  internal cross-fold validation to find the optimum lambda at which the deviance is minimised (Default=list(n_folds=10))
#' @param verbose show Lasso regression messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'      + $list_perf:
#'          - $mbe: mean bias error
#'          - $mae: mean absolute error
#'          - $rmse: root mean squared error
#'          - $r2: coefficient of determination
#'          - $corr: Pearson's product moment correlation
#'          - $power_t10: fraction of observed top 10 phenotype values correctly predicted
#'          - $power_b10: fraction of observed bottom 10 phenotype values correctly predicted
#'          - $var_additive: variance of predicted phenotype values (estimator of additive genetic variance)
#'          - $var_residual: variance of differnce between observed and predicted phenotype values (estimator of residual variance)
#'          - $h2: narrow-sense heritability estimate
#'      + $df_y_validation:
#'          - $id: names of the samples/entries/pools, 
#'          - $pop: population from which the sample/entry/pool belongs to
#'          - $y_true: observed phenotype values
#'          - $y_pred: predicted phenotype values
#'      + $vec_effects: named numeric vector of estimated effects, where the names correspond to the
#'          SNP/allele identity including chromosome/scaffold, position and optionally allele.
#'      + $n_non_zero: number of non-zero estimated effects (effects greater than machine epsilon ~2.2e-16)
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
#' n = nrow(list_merged$G)
#' vec_idx_training = sample(c(1:n), floor(n/2))
#' vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
#' list_ols = fn_ols(list_merged, vec_idx_training, vec_idx_validation, verbose=TRUE)
#' @export
fn_lasso = function(list_merged, vec_idx_training, vec_idx_validation, other_params=list(n_folds=10), verbose=FALSE) {
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # COVAR = G %*% t(G); # rownames(COVAR) = NULL
    # list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    # n = nrow(list_merged$G)
    # vec_idx_training = sample(c(1:n), floor(n/2))
    # vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
    # other_params=list(n_folds=10)
    # verbose = TRUE
    ###################################################
    ### Input sanity check
    if (methods::is(list_merged, "gpError")) {
        error = chain(list_merged, 
            methods::new("gpError",
                code=406,
                message=paste0(
                    "Error in models::fn_lasso(...). ",
                    "Input data (list_merged) is an error type."
                )))
        return(error)
    }
    if (is.logical(c(vec_idx_training, vec_idx_training)) |
        (sum(is.na(c(vec_idx_training, vec_idx_training))) > 0) |
        (min(c(vec_idx_training, vec_idx_training)) < 0) |
        (max(c(vec_idx_training, vec_idx_training)) > nrow(list_merged$G))) {
            error = methods::new("gpError",
                code=407,
                message=paste0(
                    "Error in models::fn_lasso(...). ",
                    "Please make sure that the vector of indexes for the training and validation sets are not booleans. ",
                    "We require the indexes to be positive integers without missing values. ",
                    "Also, the maximum index (", max(c(vec_idx_training, vec_idx_training)), 
                    ") should be less than or equal to the total number of samples/entries/pools in the input data: (", 
                    nrow(list_merged$G), ")."))
            return(error)
    }
    X_training = list_merged$G[vec_idx_training, , drop=FALSE]
    y_training = list_merged$list_pheno$y[vec_idx_training]
    X_training = X_training[!is.na(y_training), , drop=FALSE] ### Remove missing phenotype data from the training set
    y_training = y_training[!is.na(y_training)] ### Remove missing phenotype data from the training set
    X_validation = list_merged$G[vec_idx_validation, , drop=FALSE]
    y_validation = list_merged$list_pheno$y[vec_idx_validation]
    ### Adding covariate to explanatory matrix
    if (!is.null(list_merged$COVAR)) {
        X_training = cbind(list_merged$COVAR[vec_idx_training, , drop=FALSE], X_training)
        X_validation = cbind(list_merged$COVAR[vec_idx_validation, , drop=FALSE], X_validation)
    }
    ### Solve via Least absolute shrinkage selection operator (Lasso) regularisation
    sol = tryCatch(
        glmnet::cv.glmnet(x=X_training, y=y_training, alpha=1, nfolds=other_params$n_folds, parallel=FALSE),
        error = function(e) {NA})
    if (is.na(sol)) {
        return(list(
            list_perf=NA,
            df_y_validation=NA,
            vec_effects=NA,
            n_non_zero=NA))
    }
    ### Find the first lambda with the lowest squared error (deviance) while having non-zero SNP effects
    vec_idx_decreasing_deviance = order(sol$glmnet.fit$dev.ratio, decreasing=FALSE)
    idx_start = which(sol$lambda == sol$lambda.min)[1]
    idx_start_vec_idx_decreasing_deviance = which(vec_idx_decreasing_deviance == idx_start)
    vec_idx_decreasing_deviance = vec_idx_decreasing_deviance[idx_start_vec_idx_decreasing_deviance:length(vec_idx_decreasing_deviance)]
    vec_idx_decreasing_deviance = vec_idx_decreasing_deviance[!is.na(vec_idx_decreasing_deviance)]
    for (i in vec_idx_decreasing_deviance) {
        # i = idx_start
        b_hat = c(sol$glmnet.fit$a0[i], sol$glmnet.fit$beta[, i])
        names(b_hat) = c("intercept", colnames(X_training))
        y_pred = cbind(rep(1, length(vec_idx_validation)), X_validation) %*% b_hat
        n_non_zero = sum(abs(b_hat) >= .Machine$double.eps)
        ### Skip models with only intercept effects
        if (n_non_zero == 1) {
            next
        } else {
            break
        }
    }
    if (verbose) {
        p = ncol(X_training)
        print("Allele effects distribution:")
        tryCatch(txtplot::txtdensity(b_hat[!is.na(b_hat) & !is.infinite(b_hat)]), error=function(e){print("Empty plot")})
        print("Relative positions of allele effects across the genome:")
        tryCatch(txtplot::txtplot(b_hat[!is.na(b_hat) & !is.infinite(b_hat)]), error=function(e){print("Empty plot")})
        print(paste0("Number of non-zero effects: ", n_non_zero, " (", round(100*n_non_zero/p), "%)"))
    }
    ### Evaluate prediction performance
    df_y_validation = merge(
        data.frame(id=names(y_validation), pop=list_merged$list_pheno$pop[vec_idx_validation], y_true=y_validation), 
        data.frame(id=rownames(y_pred), y_pred=y_pred),
        by="id")
    list_perf = fn_prediction_performance_metrics(y_true=df_y_validation$y_true, y_pred=df_y_validation$y_pred, verbose=verbose)
    ### Output
    return(list(
        list_perf=list_perf,
        df_y_validation=df_y_validation,
        vec_effects=b_hat,
        n_non_zero=n_non_zero))
}

#' Elastic-net regression model (regularisation via combination of ridge and Lasso; alpha is optimised)
#'
#' @param list_merged list of merged genotype matrix, and phenotype vector, as well as an optional covariate matrix
#'  - $G: numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'      Row names can be any string of characters which identify the sample or entry or pool names.
#'      Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'      the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'      subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'  - $list_pheno:
#'      + $y: named vector of numeric phenotype data
#'      + $pop: population or groupings corresponding to each element of y
#'      + $trait_name: name of the trait
#'  - $COVAR: numeric n samples x k covariates matrix with non-null row and column names.
#' @param vec_idx_training vector of numeric indexes referring to the training set
#' @param vec_idx_validation vector of numeric indexes referring to the validation set
#' @param other_params list of additional parameters, specifically $n_folds which refers to the number of
#'  internal cross-fold validation to find the optimum lambda at which the deviance is minimised (Default=list(n_folds=10))
#' @param verbose show elastic-net regression messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'      + $list_perf:
#'          - $mbe: mean bias error
#'          - $mae: mean absolute error
#'          - $rmse: root mean squared error
#'          - $r2: coefficient of determination
#'          - $corr: Pearson's product moment correlation
#'          - $power_t10: fraction of observed top 10 phenotype values correctly predicted
#'          - $power_b10: fraction of observed bottom 10 phenotype values correctly predicted
#'          - $var_additive: variance of predicted phenotype values (estimator of additive genetic variance)
#'          - $var_residual: variance of differnce between observed and predicted phenotype values (estimator of residual variance)
#'          - $h2: narrow-sense heritability estimate
#'      + $df_y_validation:
#'          - $id: names of the samples/entries/pools, 
#'          - $pop: population from which the sample/entry/pool belongs to
#'          - $y_true: observed phenotype values
#'          - $y_pred: predicted phenotype values
#'      + $vec_effects: named numeric vector of estimated effects, where the names correspond to the
#'          SNP/allele identity including chromosome/scaffold, position and optionally allele.
#'      + $n_non_zero: number of non-zero estimated effects (effects greater than machine epsilon ~2.2e-16)
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
#' n = nrow(list_merged$G)
#' vec_idx_training = sample(c(1:n), floor(n/2))
#' vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
#' list_ols = fn_ols(list_merged, vec_idx_training, vec_idx_validation, verbose=TRUE)
#' @export
fn_elastic_net = function(list_merged, vec_idx_training, vec_idx_validation, other_params=list(n_folds=10), verbose=FALSE) {
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # COVAR = G %*% t(G); # rownames(COVAR) = NULL
    # list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    # n = nrow(list_merged$G)
    # vec_idx_training = sample(c(1:n), floor(n/2))
    # vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
    # other_params=list(n_folds=10)
    # verbose = TRUE
    ###################################################
    ### Input sanity check
    if (methods::is(list_merged, "gpError")) {
        error = chain(list_merged, 
            methods::new("gpError",
                code=408,
                message=paste0(
                    "Error in models::fn_elastic_net(...). ",
                    "Input data (list_merged) is an error type."
                )))
        return(error)
    }
    if (is.logical(c(vec_idx_training, vec_idx_training)) |
        (sum(is.na(c(vec_idx_training, vec_idx_training))) > 0) |
        (min(c(vec_idx_training, vec_idx_training)) < 0) |
        (max(c(vec_idx_training, vec_idx_training)) > nrow(list_merged$G))) {
            error = methods::new("gpError",
                code=409,
                message=paste0(
                    "Error in models::fn_elastic_net(...). ",
                    "Please make sure that the vector of indexes for the training and validation sets are not booleans. ",
                    "We require the indexes to be positive integers without missing values. ",
                    "Also, the maximum index (", max(c(vec_idx_training, vec_idx_training)), 
                    ") should be less than or equal to the total number of samples/entries/pools in the input data: (", 
                    nrow(list_merged$G), ")."))
            return(error)
    }
    X_training = list_merged$G[vec_idx_training, , drop=FALSE]
    y_training = list_merged$list_pheno$y[vec_idx_training]
    X_training = X_training[!is.na(y_training), , drop=FALSE] ### Remove missing phenotype data from the training set
    y_training = y_training[!is.na(y_training)] ### Remove missing phenotype data from the training set
    X_validation = list_merged$G[vec_idx_validation, , drop=FALSE]
    y_validation = list_merged$list_pheno$y[vec_idx_validation]
    ### Adding covariate to explanatory matrix
    if (!is.null(list_merged$COVAR)) {
        X_training = cbind(list_merged$COVAR[vec_idx_training, , drop=FALSE], X_training)
        X_validation = cbind(list_merged$COVAR[vec_idx_validation, , drop=FALSE], X_validation)
    }
    ### Solve via Elastic-net regularisation
    sol = tryCatch(
        glmnet::cv.glmnet(x=X_training, y=y_training),
        error = function(e) {NA})
    if (is.na(sol)) {
        return(list(
            list_perf=NA,
            df_y_validation=NA,
            vec_effects=NA,
            n_non_zero=NA))
    }
    ### Find the first lambda with the lowest squared error (deviance) while having non-zero SNP effects
    vec_idx_decreasing_deviance = order(sol$glmnet.fit$dev.ratio, decreasing=FALSE)
    idx_start = which(sol$lambda == sol$lambda.min)[1]
    idx_start_vec_idx_decreasing_deviance = which(vec_idx_decreasing_deviance == idx_start)
    vec_idx_decreasing_deviance = vec_idx_decreasing_deviance[idx_start_vec_idx_decreasing_deviance:length(vec_idx_decreasing_deviance)]
    vec_idx_decreasing_deviance = vec_idx_decreasing_deviance[!is.na(vec_idx_decreasing_deviance)]
    for (i in vec_idx_decreasing_deviance) {
        # i = idx_start
        b_hat = c(sol$glmnet.fit$a0[i], sol$glmnet.fit$beta[, i])
        names(b_hat) = c("intercept", colnames(X_training))
        y_pred = cbind(rep(1, length(vec_idx_validation)), X_validation) %*% b_hat
        n_non_zero = sum(abs(b_hat) >= .Machine$double.eps)
        ### Skip models with only intercept effects
        if (n_non_zero == 1) {
            next
        } else {
            break
        }
    }
    if (verbose) {
        p = ncol(X_training)
        print("Allele effects distribution:")
        tryCatch(txtplot::txtdensity(b_hat[!is.na(b_hat) & !is.infinite(b_hat)]), error=function(e){print("Empty plot")})
        print("Relative positions of allele effects across the genome:")
        tryCatch(txtplot::txtplot(b_hat[!is.na(b_hat) & !is.infinite(b_hat)]), error=function(e){print("Empty plot")})
        print(paste0("Number of non-zero effects: ", n_non_zero, " (", round(100*n_non_zero/p), "%)"))
    }
    ### Evaluate prediction performance
    df_y_validation = merge(
        data.frame(id=names(y_validation), pop=list_merged$list_pheno$pop[vec_idx_validation], y_true=y_validation), 
        data.frame(id=rownames(y_pred), y_pred=y_pred),
        by="id")
    list_perf = fn_prediction_performance_metrics(y_true=df_y_validation$y_true, y_pred=df_y_validation$y_pred, verbose=verbose)
    ### Output
    return(list(
        list_perf=list_perf,
        df_y_validation=df_y_validation,
        vec_effects=b_hat,
        n_non_zero=n_non_zero))
}

#' Bayes A regression model (scaled t-distributed effects)
#'
#' @param list_merged list of merged genotype matrix, and phenotype vector, as well as an optional covariate matrix
#'  - $G: numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'      Row names can be any string of characters which identify the sample or entry or pool names.
#'      Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'      the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'      subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'  - $list_pheno:
#'      + $y: named vector of numeric phenotype data
#'      + $pop: population or groupings corresponding to each element of y
#'      + $trait_name: name of the trait
#'  - $COVAR: numeric n samples x k covariates matrix with non-null row and column names.
#' @param vec_idx_training vector of numeric indexes referring to the training set
#' @param vec_idx_validation vector of numeric indexes referring to the validation set
#' @param other_params list of additional parameters, 
#'  - other_params$nIter: total number of iterations (Default=12e3)
#'  - other_params$burnIn: number of burn-in iterations (Default=2e3)
#'  - other_params$out_prefix: prefix of temporary output files (Default="bglr_bayesA-")
#' @param verbose show Bayes A regression messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'      + $list_perf:
#'          - $mbe: mean bias error
#'          - $mae: mean absolute error
#'          - $rmse: root mean squared error
#'          - $r2: coefficient of determination
#'          - $corr: Pearson's product moment correlation
#'          - $power_t10: fraction of observed top 10 phenotype values correctly predicted
#'          - $power_b10: fraction of observed bottom 10 phenotype values correctly predicted
#'          - $var_additive: variance of predicted phenotype values (estimator of additive genetic variance)
#'          - $var_residual: variance of differnce between observed and predicted phenotype values (estimator of residual variance)
#'          - $h2: narrow-sense heritability estimate
#'      + $df_y_validation:
#'          - $id: names of the samples/entries/pools, 
#'          - $pop: population from which the sample/entry/pool belongs to
#'          - $y_true: observed phenotype values
#'          - $y_pred: predicted phenotype values
#'      + $vec_effects: named numeric vector of estimated effects, where the names correspond to the
#'          SNP/allele identity including chromosome/scaffold, position and optionally allele.
#'      + $n_non_zero: number of non-zero estimated effects (effects greater than machine epsilon ~2.2e-16)
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
#' n = nrow(list_merged$G)
#' vec_idx_training = sample(c(1:n), floor(n/2))
#' vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
#' list_ols = fn_ols(list_merged, vec_idx_training, vec_idx_validation, verbose=TRUE)
#' @export
fn_Bayes_A = function(list_merged, vec_idx_training, vec_idx_validation, 
    other_params=list(nIter=12e3, burnIn=2e3, out_prefix="bglr_bayesA-"), verbose=FALSE) {
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # COVAR = G %*% t(G); # rownames(COVAR) = NULL
    # list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    # n = nrow(list_merged$G)
    # vec_idx_training = sample(c(1:n), floor(n/2))
    # vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
    # other_params=list(
    #     nIter=1000,
    #     burnIn=100,
    #     out_prefix="bglr_bayesA",
    #     covariate=matrix(stats::rnorm(2*n), ncol=2))
    # verbose = TRUE
    ###################################################
    ### Input sanity check
    if (methods::is(list_merged, "gpError")) {
        error = chain(list_merged, 
            methods::new("gpError",
                code=410,
                message=paste0(
                    "Error in models::fn_elastic_net(...). ",
                    "Input data (list_merged) is an error type."
                )))
        return(error)
    }
    if (is.logical(c(vec_idx_training, vec_idx_training)) |
        (sum(is.na(c(vec_idx_training, vec_idx_training))) > 0) |
        (min(c(vec_idx_training, vec_idx_training)) < 0) |
        (max(c(vec_idx_training, vec_idx_training)) > nrow(list_merged$G))) {
            error = methods::new("gpError",
                code=411,
                message=paste0(
                    "Error in models::fn_Bayes_A(...). ",
                    "Please make sure that the vector of indexes for the training and validation sets are not booleans. ",
                    "We require the indexes to be positive integers without missing values. ",
                    "Also, the maximum index (", max(c(vec_idx_training, vec_idx_training)), 
                    ") should be less than or equal to the total number of samples/entries/pools in the input data: (", 
                    nrow(list_merged$G), ")."))
            return(error)
    }
    ### Set validation set to missing
    yNA = list_merged$list_pheno$y
    yNA[vec_idx_validation] = NA
    ### Define Bayes A regression function including covariates if they exist
    if (!is.null(list_merged$COVAR)) {
        ETA = list(MRK=list(X=list_merged$G, model="BayesA"),
                       list(X=list_merged$COVAR, model="FIXED"))
    } else {
        ETA = list(MRK=list(X=list_merged$G, model="BayesA", saveEffects=TRUE))
    }
    ### Attempt at preventing overwrites to the running Gibbs samplers in parallel
    other_params$out_prefix = gsub(":", ".", gsub(" ", "-", paste(other_params$out_prefix, Sys.time(), stats::runif(1), sep="-")))
    ### Solve via Bayes A (a priori assume heritability at 50%, i.e. R2=0.5 below)
    sol = BGLR::BGLR(y=yNA, ETA=ETA, R2=0.5, nIter=other_params$nIter, burnIn=other_params$burnIn, saveAt=other_params$out_prefix, verbose=FALSE)
    ### Extract effects including the intercept and fixed effects
    if (!is.null(list_merged$COVAR)) {
        b_hat = c(sol$mu, sol$ETA[[2]]$b, sol$ETA$MRK$b)
    } else {
        b_hat = c(sol$mu, sol$ETA$MRK$b)
    }
    names(b_hat)[1] = "intercept"
    n_non_zero = sum(abs(b_hat) >= .Machine$double.eps)
    if (verbose) {
        p = ncol(list_merged$G)
        print("Allele effects distribution:")
        tryCatch(txtplot::txtdensity(b_hat[!is.na(b_hat) & !is.infinite(b_hat)]), error=function(e){print("Empty plot")})
        print("Relative positions of allele effects across the genome:")
        tryCatch(txtplot::txtplot(b_hat[!is.na(b_hat) & !is.infinite(b_hat)]), error=function(e){print("Empty plot")})
        print(paste0("Number of non-zero effects: ", n_non_zero, " (", round(100*n_non_zero/p), "%)"))
    }
    ### Evaluate prediction performance
    df_y_validation = merge(
        data.frame(id=names(list_merged$list_pheno$y[vec_idx_validation]), pop=list_merged$list_pheno$pop[vec_idx_validation], y_true=list_merged$list_pheno$y[vec_idx_validation]), 
        data.frame(id=names(yNA)[vec_idx_validation], y_pred=sol$yHat[vec_idx_validation]),
        by="id")
    list_perf = fn_prediction_performance_metrics(y_true=df_y_validation$y_true, y_pred=df_y_validation$y_pred, verbose=verbose)
    ### Clean-up temporary files
    for (f in list.files(dirname(other_params$out_prefix), pattern=basename(other_params$out_prefix))) {
        file.remove(file.path(dirname(other_params$out_prefix), f))
    }
    ### Output
    return(list(
        list_perf=list_perf,
        df_y_validation=df_y_validation,
        vec_effects=b_hat,
        n_non_zero=n_non_zero))
}

#' Bayes B regression model 
#'  (scaled t-distributed effects with probability \eqn{\pi}; and zero effects with probability \eqn{1-\pi}, 
#'  where \eqn{\pi \sim \beta(\theta_1, \theta_2)})
#'
#' @param list_merged list of merged genotype matrix, and phenotype vector, as well as an optional covariate matrix
#'  - $G: numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'      Row names can be any string of characters which identify the sample or entry or pool names.
#'      Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'      the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'      subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'  - $list_pheno:
#'      + $y: named vector of numeric phenotype data
#'      + $pop: population or groupings corresponding to each element of y
#'      + $trait_name: name of the trait
#'  - $COVAR: numeric n samples x k covariates matrix with non-null row and column names.
#' @param vec_idx_training vector of numeric indexes referring to the training set
#' @param vec_idx_validation vector of numeric indexes referring to the validation set
#' @param other_params list of additional parameters, 
#'  - other_params$nIter: total number of iterations (Default=12e3)
#'  - other_params$burnIn: number of burn-in iterations (Default=2e3)
#'  - other_params$out_prefix: prefix of temporary output files (Default="bglr_bayesB-")
#' @param verbose show Bayes B regression messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'      + $list_perf:
#'          - $mbe: mean bias error
#'          - $mae: mean absolute error
#'          - $rmse: root mean squared error
#'          - $r2: coefficient of determination
#'          - $corr: Pearson's product moment correlation
#'          - $power_t10: fraction of observed top 10 phenotype values correctly predicted
#'          - $power_b10: fraction of observed bottom 10 phenotype values correctly predicted
#'          - $var_additive: variance of predicted phenotype values (estimator of additive genetic variance)
#'          - $var_residual: variance of differnce between observed and predicted phenotype values (estimator of residual variance)
#'          - $h2: narrow-sense heritability estimate
#'      + $df_y_validation:
#'          - $id: names of the samples/entries/pools, 
#'          - $pop: population from which the sample/entry/pool belongs to
#'          - $y_true: observed phenotype values
#'          - $y_pred: predicted phenotype values
#'      + $vec_effects: named numeric vector of estimated effects, where the names correspond to the
#'          SNP/allele identity including chromosome/scaffold, position and optionally allele.
#'      + $n_non_zero: number of non-zero estimated effects (effects greater than machine epsilon ~2.2e-16)
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
#' n = nrow(list_merged$G)
#' vec_idx_training = sample(c(1:n), floor(n/2))
#' vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
#' list_ols = fn_ols(list_merged, vec_idx_training, vec_idx_validation, verbose=TRUE)
#' @export
fn_Bayes_B = function(list_merged, vec_idx_training, vec_idx_validation, 
    other_params=list(nIter=12e3, burnIn=2e3, out_prefix="bglr_bayesB-"), verbose=FALSE) {
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # COVAR = G %*% t(G); # rownames(COVAR) = NULL
    # list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    # n = nrow(list_merged$G)
    # vec_idx_training = sample(c(1:n), floor(n/2))
    # vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
    # other_params=list(
    #     nIter=1000,
    #     burnIn=100,
    #     out_prefix="bglr_bayesB",
    #     covariate=matrix(stats::rnorm(2*n), ncol=2))
    # verbose = TRUE
    ###################################################
    ### Input sanity check
    if (methods::is(list_merged, "gpError")) {
        error = chain(list_merged, 
            methods::new("gpError",
                code=412,
                message=paste0(
                    "Error in models::fn_elastic_net(...). ",
                    "Input data (list_merged) is an error type."
                )))
        return(error)
    }
    if (is.logical(c(vec_idx_training, vec_idx_training)) |
        (sum(is.na(c(vec_idx_training, vec_idx_training))) > 0) |
        (min(c(vec_idx_training, vec_idx_training)) < 0) |
        (max(c(vec_idx_training, vec_idx_training)) > nrow(list_merged$G))) {
            error = methods::new("gpError",
                code=413,
                message=paste0(
                    "Error in models::fn_Bayes_B(...). ",
                    "Please make sure that the vector of indexes for the training and validation sets are not booleans. ",
                    "We require the indexes to be positive integers without missing values. ",
                    "Also, the maximum index (", max(c(vec_idx_training, vec_idx_training)), 
                    ") should be less than or equal to the total number of samples/entries/pools in the input data: (", 
                    nrow(list_merged$G), ")."))
            return(error)
    }
    ### Set validation set to missing
    yNA = list_merged$list_pheno$y
    yNA[vec_idx_validation] = NA
    ### Define Bayes B regression function including covariates if they exist
    if (!is.null(list_merged$COVAR)) {
        ETA = list(MRK=list(X=list_merged$G, model="BayesB"),
                       list(X=list_merged$COVAR, model="FIXED"))
    } else {
        ETA = list(MRK=list(X=list_merged$G, model="BayesB"))
    }
    ### Attempt at preventing overwrites to the running Gibbs samplers in parallel
    other_params$out_prefix = gsub(":", ".", gsub(" ", "-", paste(other_params$out_prefix, Sys.time(), stats::runif(1), sep="-")))
    ### Solve via Bayes B (a priori assume heritability at 50%, i.e. R2=0.5 below)
    sol = BGLR::BGLR(y=yNA, ETA=ETA, R2=0.5, nIter=other_params$nIter, burnIn=other_params$burnIn, saveAt=other_params$out_prefix, verbose=FALSE)
    ### Extract effects including the intercept and fixed effects
    if (!is.null(list_merged$COVAR)) {
        b_hat = c(sol$mu, sol$ETA[[2]]$b, sol$ETA$MRK$b)
    } else {
        b_hat = c(sol$mu, sol$ETA$MRK$b)
    }
    names(b_hat)[1] = "intercept"
    n_non_zero = sum(abs(b_hat) >= .Machine$double.eps)
    if (verbose) {
        p = ncol(list_merged$G)
        print("Allele effects distribution:")
        tryCatch(txtplot::txtdensity(b_hat[!is.na(b_hat) & !is.infinite(b_hat)]), error=function(e){print("Empty plot")})
        print("Relative positions of allele effects across the genome:")
        tryCatch(txtplot::txtplot(b_hat[!is.na(b_hat) & !is.infinite(b_hat)]), error=function(e){print("Empty plot")})
        print(paste0("Number of non-zero effects: ", n_non_zero, " (", round(100*n_non_zero/p), "%)"))
    }
    ### Evaluate prediction performance
    df_y_validation = merge(
        data.frame(id=names(list_merged$list_pheno$y[vec_idx_validation]), pop=list_merged$list_pheno$pop[vec_idx_validation], y_true=list_merged$list_pheno$y[vec_idx_validation]), 
        data.frame(id=names(yNA)[vec_idx_validation], y_pred=sol$yHat[vec_idx_validation]),
        by="id")
    list_perf = fn_prediction_performance_metrics(y_true=df_y_validation$y_true, y_pred=df_y_validation$y_pred, verbose=verbose)
    ### Clean-up temporary files
    for (f in list.files(dirname(other_params$out_prefix), pattern=basename(other_params$out_prefix))) {
        file.remove(file.path(dirname(other_params$out_prefix), f))
    }
    ### Output
    return(list(
        list_perf=list_perf,
        df_y_validation=df_y_validation,
        vec_effects=b_hat,
        n_non_zero=n_non_zero))
}

#' Bayes C regression model
#'  (normally distributed effects (\eqn{N(0, \sigma^2_{\beta})}) with probability \eqn{\pi}; and zero effects 
#'  with probability \eqn{1-\pi}, where \eqn{\pi \sim \beta(\theta_1, \theta_2)})
#'
#' @param list_merged list of merged genotype matrix, and phenotype vector, as well as an optional covariate matrix
#'  - $G: numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'      Row names can be any string of characters which identify the sample or entry or pool names.
#'      Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'      the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'      subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'  - $list_pheno:
#'      + $y: named vector of numeric phenotype data
#'      + $pop: population or groupings corresponding to each element of y
#'      + $trait_name: name of the trait
#'  - $COVAR: numeric n samples x k covariates matrix with non-null row and column names.
#' @param vec_idx_training vector of numeric indexes referring to the training set
#' @param vec_idx_validation vector of numeric indexes referring to the validation set
#' @param other_params list of additional parameters, 
#'  - other_params$nIter: total number of iterations (Default=12e3)
#'  - other_params$burnIn: number of burn-in iterations (Default=2e3)
#'  - other_params$out_prefix: prefix of temporary output files (Default="bglr_bayesB-")
#' @param verbose show Bayes C regression messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'      + $list_perf:
#'          - $mbe: mean bias error
#'          - $mae: mean absolute error
#'          - $rmse: root mean squared error
#'          - $r2: coefficient of determination
#'          - $corr: Pearson's product moment correlation
#'          - $power_t10: fraction of observed top 10 phenotype values correctly predicted
#'          - $power_b10: fraction of observed bottom 10 phenotype values correctly predicted
#'          - $var_additive: variance of predicted phenotype values (estimator of additive genetic variance)
#'          - $var_residual: variance of differnce between observed and predicted phenotype values (estimator of residual variance)
#'          - $h2: narrow-sense heritability estimate
#'      + $df_y_validation:
#'          - $id: names of the samples/entries/pools, 
#'          - $pop: population from which the sample/entry/pool belongs to
#'          - $y_true: observed phenotype values
#'          - $y_pred: predicted phenotype values
#'      + $vec_effects: named numeric vector of estimated effects, where the names correspond to the
#'          SNP/allele identity including chromosome/scaffold, position and optionally allele.
#'      + $n_non_zero: number of non-zero estimated effects (effects greater than machine epsilon ~2.2e-16)
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
#' n = nrow(list_merged$G)
#' vec_idx_training = sample(c(1:n), floor(n/2))
#' vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
#' list_ols = fn_ols(list_merged, vec_idx_training, vec_idx_validation, verbose=TRUE)
#' @export
fn_Bayes_C = function(list_merged, vec_idx_training, vec_idx_validation, 
    other_params=list(nIter=12e3, burnIn=2e3, out_prefix="bglr_bayesC-"), verbose=FALSE) {
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # COVAR = G %*% t(G); # rownames(COVAR) = NULL
    # list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    # n = nrow(list_merged$G)
    # vec_idx_training = sample(c(1:n), floor(n/2))
    # vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
    # other_params=list(
    #     nIter=1000,
    #     burnIn=100,
    #     out_prefix="bglr_bayesC",
    #     covariate=matrix(stats::rnorm(2*n), ncol=2))
    # verbose = TRUE
    ###################################################
    ### Input sanity check
    if (methods::is(list_merged, "gpError")) {
        error = chain(list_merged, 
            methods::new("gpError",
                code=414,
                message=paste0(
                    "Error in models::fn_elastic_net(...). ",
                    "Input data (list_merged) is an error type."
                )))
        return(error)
    }
    if (is.logical(c(vec_idx_training, vec_idx_training)) |
        (sum(is.na(c(vec_idx_training, vec_idx_training))) > 0) |
        (min(c(vec_idx_training, vec_idx_training)) < 0) |
        (max(c(vec_idx_training, vec_idx_training)) > nrow(list_merged$G))) {
            error = methods::new("gpError",
                code=415,
                message=paste0(
                    "Error in models::fn_Bayes_C(...). ",
                    "Please make sure that the vector of indexes for the training and validation sets are not booleans. ",
                    "We require the indexes to be positive integers without missing values. ",
                    "Also, the maximum index (", max(c(vec_idx_training, vec_idx_training)), 
                    ") should be less than or equal to the total number of samples/entries/pools in the input data: (", 
                    nrow(list_merged$G), ")."))
            return(error)
    }
    ### Set validation set to missing
    yNA = list_merged$list_pheno$y
    yNA[vec_idx_validation] = NA
    ### Define Bayes C regression function including covariates if they exist
    if (!is.null(list_merged$COVAR)) {
        ETA = list(MRK=list(X=list_merged$G, model="BayesC"),
                       list(X=list_merged$COVAR, model="FIXED"))
    } else {
        ETA = list(MRK=list(X=list_merged$G, model="BayesC"))
    }
    ### Attempt at preventing overwrites to the running Gibbs samplers in parallel
    other_params$out_prefix = gsub(":", ".", gsub(" ", "-", paste(other_params$out_prefix, Sys.time(), stats::runif(1), sep="-")))
    ### Solve via Bayes C (a priori assume heritability at 50%, i.e. R2=0.5 below)
    sol = BGLR::BGLR(y=yNA, ETA=ETA, R2=0.5, nIter=other_params$nIter, burnIn=other_params$burnIn, saveAt=other_params$out_prefix, verbose=FALSE)
    ### Extract effects including the intercept and fixed effects
    if (!is.null(list_merged$COVAR)) {
        b_hat = c(sol$mu, sol$ETA[[2]]$b, sol$ETA$MRK$b)
    } else {
        b_hat = c(sol$mu, sol$ETA$MRK$b)
    }
    names(b_hat)[1] = "intercept"
    n_non_zero = sum(abs(b_hat) >= .Machine$double.eps)
    if (verbose) {
        p = ncol(list_merged$G)
        print("Allele effects distribution:")
        tryCatch(txtplot::txtdensity(b_hat[!is.na(b_hat) & !is.infinite(b_hat)]), error=function(e){print("Empty plot")})
        print("Relative positions of allele effects across the genome:")
        tryCatch(txtplot::txtplot(b_hat[!is.na(b_hat) & !is.infinite(b_hat)]), error=function(e){print("Empty plot")})
        print(paste0("Number of non-zero effects: ", n_non_zero, " (", round(100*n_non_zero/p), "%)"))
    }
    ### Evaluate prediction performance
    df_y_validation = merge(
        data.frame(id=names(list_merged$list_pheno$y[vec_idx_validation]), pop=list_merged$list_pheno$pop[vec_idx_validation], y_true=list_merged$list_pheno$y[vec_idx_validation]), 
        data.frame(id=names(yNA)[vec_idx_validation], y_pred=sol$yHat[vec_idx_validation]),
        by="id")
    list_perf = fn_prediction_performance_metrics(y_true=df_y_validation$y_true, y_pred=df_y_validation$y_pred, verbose=verbose)
    ### Clean-up temporary files
    for (f in list.files(dirname(other_params$out_prefix), pattern=basename(other_params$out_prefix))) {
        file.remove(file.path(dirname(other_params$out_prefix), f))
    }
    ### Output
    return(list(
        list_perf=list_perf,
        df_y_validation=df_y_validation,
        vec_effects=b_hat,
        n_non_zero=n_non_zero))
}

#' Mixed linear genotype value model (gBLUP)
#'  (where genomic information is only used to estimate genetic relationships between samples/entries/pools)
#'  and not to estimate the additive effects of each SNP/locus/allele which differentiates it from GBLUP.
#'  The genotypic value of each sample/entry/pool is estimated as the best linear unbiased predictors or BLUPs.)
#'
#' @param list_merged list of merged genotype matrix, and phenotype vector, as well as an optional covariate matrix
#'  - $G: numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'      Row names can be any string of characters which identify the sample or entry or pool names.
#'      Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'      the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'      subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'  - $list_pheno:
#'      + $y: named vector of numeric phenotype data
#'      + $pop: population or groupings corresponding to each element of y
#'      + $trait_name: name of the trait
#'  - $COVAR: numeric n samples x k covariates matrix with non-null row and column names.
#' @param vec_idx_training vector of numeric indexes referring to the training set
#' @param vec_idx_validation vector of numeric indexes referring to the validation set
#' @param other_params list of additional parameters which is NULL
#' @param verbose show gBLUP regression messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'      + $list_perf:
#'          - $mbe: mean bias error
#'          - $mae: mean absolute error
#'          - $rmse: root mean squared error
#'          - $r2: coefficient of determination
#'          - $corr: Pearson's product moment correlation
#'          - $power_t10: fraction of observed top 10 phenotype values correctly predicted
#'          - $power_b10: fraction of observed bottom 10 phenotype values correctly predicted
#'          - $var_additive: variance of predicted phenotype values (estimator of additive genetic variance)
#'          - $var_residual: variance of differnce between observed and predicted phenotype values (estimator of residual variance)
#'          - $h2: narrow-sense heritability estimate
#'      + $df_y_validation:
#'          - $id: names of the samples/entries/pools, 
#'          - $pop: population from which the sample/entry/pool belongs to
#'          - $y_true: observed phenotype values
#'          - $y_pred: predicted phenotype values
#'      + $vec_effects: named numeric vector of estimated effects, where the names correspond to the
#'          SNP/allele identity including chromosome/scaffold, position and optionally allele.
#'      + $n_non_zero: number of non-zero estimated fixed (intercept and covariate/s) and random (genotype values)
#'           effects (effects greater than machine epsilon ~2.2e-16)
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
#' n = nrow(list_merged$G)
#' vec_idx_training = sample(c(1:n), floor(n/2))
#' vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
#' list_ols = fn_ols(list_merged, vec_idx_training, vec_idx_validation, verbose=TRUE)
#' @export
fn_gBLUP = function(list_merged, vec_idx_training, vec_idx_validation, other_params=list(), verbose=FALSE) {
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # COVAR = G %*% t(G); # rownames(COVAR) = NULL
    # list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    # n = nrow(list_merged$G)
    # vec_idx_training = sample(c(1:n), floor(n/2))
    # vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
    # other_params=list()
    # verbose = TRUE
    ###################################################
    ### Input sanity check
    if (methods::is(list_merged, "gpError")) {
        error = chain(list_merged, 
            methods::new("gpError",
                code=416,
                message=paste0(
                    "Error in models::fn_elastic_net(...). ",
                    "Input data (list_merged) is an error type."
                )))
        return(error)
    }
    if (is.logical(c(vec_idx_training, vec_idx_training)) |
        (sum(is.na(c(vec_idx_training, vec_idx_training))) > 0) |
        (min(c(vec_idx_training, vec_idx_training)) < 0) |
        (max(c(vec_idx_training, vec_idx_training)) > nrow(list_merged$G))) {
            error = methods::new("gpError",
                code=417,
                message=paste0(
                    "Error in models::fn_Bayes_C(...). ",
                    "Please make sure that the vector of indexes for the training and validation sets are not booleans. ",
                    "We require the indexes to be positive integers without missing values. ",
                    "Also, the maximum index (", max(c(vec_idx_training, vec_idx_training)), 
                    ") should be less than or equal to the total number of samples/entries/pools in the input data: (", 
                    nrow(list_merged$G), ")."))
            return(error)
    }
    ### Set validation set to missing
    yNA = list_merged$list_pheno$y
    yNA[vec_idx_validation] = NA
    ### Generate the genetic relationship matrix from the genotype data
    A = sommer::A.mat(list_merged$G)
    ### Build the phenotype and id data frame
    df_training = data.frame(y=yNA, id=as.factor(names(yNA)))
    colnames(df_training) = c("y", "id")
    ### Build the gBLUP model with/without fixed covariate/s
    ### Note that if we have more than 10 fixed covariate then the model fitting is more likely to fail.
    ### To prevent this we use the first 10 PCs of the fixed covariates instead
    if (!is.null(list_merged$COVAR)) {
        X = scale(list_merged$COVAR)
        if (ncol(X) > 10) {
            X = stats::prcomp(X)$x[, 1:10]
        }
        for (j in 1:ncol(X)) {
            eval(parse(text=paste0("df_training$covariate_", j, " = X[, j]")))
        }
        covariates_string = paste(paste0("covariate_", 1:ncol(X)), collapse="+")
        eval(parse(text=paste0("mod = sommer::mmer(y ~ 1 + ", covariates_string, ", random= ~vsr(id, Gu=A ), rcov= ~vsr(units), data=df_training, dateWarning=FALSE, verbose=FALSE)")))
    } else {
        mod = sommer::mmer(y ~ 1, random= ~vsr(id, Gu=A), rcov= ~vsr(units), data=df_training, dateWarning=FALSE, verbose=FALSE)
    }
    ### Extract effects
    b_hat = mod$Beta$Estimate; names(b_hat) = mod$Beta$Effect
    u_hat = mod$U$`u:id`$y
    vec_effects = c(b_hat, u_hat)
    n_non_zero = sum(abs(vec_effects) >= .Machine$double.eps)
    if (verbose) {
        p = length(vec_effects)
        if (length(b_hat) > 1) {
            print("Fixed effects distribution:")
            tryCatch(txtplot::txtdensity(b_hat[!is.na(b_hat) & !is.infinite(b_hat)]), error=function(e){print("Empty plot")})
        } else {
            print(paste0("Intercept = ", b_hat))
        }
        print("Random effects distribution, i.e. BLUPs of each sample:")
        tryCatch(txtplot::txtdensity(u_hat[!is.na(u_hat) & !is.infinite(u_hat)]), error=function(e){print("Empty plot")})
        print(paste0("Number of non-zero effects: ", n_non_zero, " (", round(100*n_non_zero/p), "%)"))
    }
    ### Evaluate prediction performance
    df_y_validation = merge(
        data.frame(id=names(list_merged$list_pheno$y[vec_idx_validation]), pop=list_merged$list_pheno$pop[vec_idx_validation], y_true=list_merged$list_pheno$y[vec_idx_validation]), 
        data.frame(id=names(mod$U$`u:id`$y), y_pred=unlist(mod$U)),
        by="id")
    list_perf = fn_prediction_performance_metrics(y_true=df_y_validation$y_true, y_pred=df_y_validation$y_pred, verbose=verbose)
    ### Output
    return(list(
        list_perf=list_perf,
        df_y_validation=df_y_validation,
        vec_effects=vec_effects,
        n_non_zero=n_non_zero))





    ### Merge expected and predicted GEBVs
    df_y_validation = data.frame(id=names(list_merged$list_pheno$y[vec_idx_validation]), list_merged$list_pheno$y[vec_idx_validation])
    df_y_predicted = data.frame(id=names(mod$U$`u:id`$y), pred=mod$U)
    df_y = merge(df_y_validation, df_y_predicted, by="id")
    colnames(df_y) = c("id", "true", "pred")
    perf = fn_prediction_performance_metrics(y_true=df_y$true, y_pred=df_y$pred)
    perf$y_pred = df_y$pred; names(perf$y_pred) = df_y$id
    perf$n_non_zero = NA; names(perf$n_non_zero) = NULL ### This is empty because individual alleles or markers were not used directly in the gBLUP model
    return(perf)
}

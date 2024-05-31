# source("R/load.R")
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

fn_ols = function(list_merged, vec_idx_training, vec_idx_validation, other_params=list(diag_inflate=1e-4)) {
    ###################################################
    ### TEST
    list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    COVAR = G %*% t(G); # rownames(COVAR) = NULL
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    n = nrow(list_merged$G)
    vec_idx_training = sample(c(1:n), floor(n/2))
    vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
    other_params=list(diag_inflate=1e-4)
    ###################################################
    if (methods::is(list_merged, "gpError")) {
        error = chain(list_merged, 
            methods::new("gpError",
                code=000,
                message=paste0(
                    "Error in models::fn_ols(list_merged, vec_idx_training, vec_idx_validation, other_params). ",
                    "Input data is an error type."
                )))
        return(error)
    }
    X_training = list_merged$G[vec_idx_training, ]
    y_training = list_merged$list_pheno$y[vec_idx_training]
    X_validation = list_merged$G[vec_idx_validation, ]
    y_validation = list_merged$list_pheno$y[vec_idx_validation]
    ### Adding covariate to explanatory matrix
    if (!is.null(list_merged$COVAR)) {
        X_training = cbind(list_merged$COVAR[vec_idx_training, ], X_training)
        X_validation = cbind(list_merged$COVAR[vec_idx_validation, ], X_validation)
    }

    if (nrow(X_training) > ncol(X_training)) {
        A = t(X_training) %*% X_training
        inv_A = tryCatch(solve(A), error=function(x) {
            B = A
            diag(B) = diag(B) + other_params$diag_inflate
            out = tryCatch(solve(B), error=function(e){return(NA)})
            return(out)
        })
        diag(A %*% inv_A)
        if (is.na(inv_A)) {
            error = methods::new("gpError",
                code=000,
                message=paste0(
                    "Error in models::fn_ols(list_merged, vec_idx_training, vec_idx_validation, other_params)"
                ))
        }
    }



    A = t(X_training) %*% X_training
    inv_A = solve()
    b = inv_A %*% t(X_training) %*% y_training
    y_pred = X_validation %*% b
    df_y = merge(data.frame(id=rownames(y)[vec_idx_validation], true=y_validation), data.frame(id=rownames(y_pred), pred=y_pred), by="id")
    perf = fn_prediction_performance_metrics(y_true=df_y$true, y_pred=df_y$pred)
    perf$y_pred = y_pred
    perf$n_non_zero = sum(abs(b) > .Machine$double.eps)
    return(perf)
}

fn_ridge = function(G, y, vec_idx_training, vec_idx_validation, other_params=list(covariate=NULL)) {
    # n = 100; n_alleles = 2
    # G = simquantgen::fn_simulate_genotypes(n=n, l=500, ploidy=2, n_alleles=n_alleles, min_allele_freq=0.001, n_chr=5, max_pos=135e6, dist_bp_at_50perc_r2=5e6, n_threads=2)
    # y = simquantgen::fn_simulate_phenotypes(G, n_alleles=n_alleles, dist_effects="norm", n_effects=10, purely_additive=TRUE, n_networks=1, n_effects_per_network=50, h2=0.75, pheno_reps=1)$Y
    # y = scale(y, center=TRUE, scale=TRUE)
    # vec_idx_validation = c(1:round(n*0.1))
    # vec_idx_training = !(c(1:n) %in% vec_idx_validation)
    # other_params=list(covariate=NULL)
    ### Make sure that our genotype and phenotype data are labelled, i.e. have row and column names
    if (is.null(rownames(G)) | is.null(colnames(G)) | is.null(rownames(y))) {
        print("Error: Genotype and/or phenotype data have no labels. Please include the row names (entry names) and column names (loci names).")
        return(-1)
    }
    X_training = G[vec_idx_training, ]
    y_training = y[vec_idx_training]
    X_validation = G[vec_idx_validation, ]
    y_validation = y[vec_idx_validation]
    ### Adding covariate to explanatory matrix
    C = other_params$covariate
    if (is.null(C) == FALSE) {
        if (nrow(C) < max(c(vec_idx_training, vec_idx_validation))) {
            print("Error: The covariate matrix is incompatible with the genotype data, i.e. there are less rows in the covariate matrix that the indices in vec_idx_training oridx_validation.")
            return(-1)
        }
        X_training = cbind(C[vec_idx_training, ], X_training)
        X_validation = cbind(C[vec_idx_validation, ], X_validation)
    }
    sol = glmnet::cv.glmnet(x=X_training, y=y_training, alpha=0) ### Ridge -> alpha = 0.0
    ### If the optimisation picks a lambda which does not provide estimates of the effects of the genotype data then we use the next best lambda
    n_non_zero = sol$nzero[sol$lambda==sol$lambda.min]
    if (n_non_zero == 0) {
        idx = which(sol$lambda==sol$lambda.min)
        for (i in 1:(length(sol$lambda)-1)) {
            if (sol$nzero[idx] == 0) {
                idx = idx + 1
            }
        }
        n_non_zero = sol$nzero[idx]
        b_hat = c(sol$glmnet.fit$a0[idx], sol$glmnet.fit$beta[, idx])
        y_pred = cbind(rep(1, length(vec_idx_validation)), X_validation) %*% b_hat
    } else {
        y_pred = stats::predict(sol, newx=X_validation, s="lambda.min")
    }
    colnames(y_pred) = "pred"
    df_y = merge(data.frame(id=rownames(y)[vec_idx_validation], true=y_validation), data.frame(id=rownames(y_pred), pred=y_pred), by="id")
    perf = fn_prediction_performance_metrics(y_true=df_y$true, y_pred=df_y$pred)
    perf$y_pred = y_pred
    perf$n_non_zero = n_non_zero; names(perf$n_non_zero) = NULL
    return(perf)
}

fn_lasso = function(G, y, vec_idx_training, vec_idx_validation, other_params=list(covariate=NULL)) {
    # n = 100; n_alleles = 2
    # G = simquantgen::fn_simulate_genotypes(n=n, l=500, ploidy=2, n_alleles=n_alleles, min_allele_freq=0.001, n_chr=5, max_pos=135e6, dist_bp_at_50perc_r2=5e6, n_threads=2)
    # y = simquantgen::fn_simulate_phenotypes(G, n_alleles=n_alleles, dist_effects="norm", n_effects=10, purely_additive=TRUE, n_networks=1, n_effects_per_network=50, h2=0.75, pheno_reps=1)$Y
    # y = scale(y, center=TRUE, scale=TRUE)
    # vec_idx_validation = c(1:round(n*0.1))
    # vec_idx_training = !(c(1:n) %in% vec_idx_validation)
    # other_params=list(covariate=NULL)
    ### Make sure that our genotype and phenotype data are labelled, i.e. have row and column names
    if (is.null(rownames(G)) | is.null(colnames(G)) | is.null(rownames(y))) {
        print("Error: Genotype and/or phenotype data have no labels. Please include the row names (entry names) and column names (loci names).")
        return(-1)
    }
    X_training = G[vec_idx_training, ]
    y_training = y[vec_idx_training]
    X_validation = G[vec_idx_validation, ]
    y_validation = y[vec_idx_validation]
    ### Adding covariate to explanatory matrix
    C = other_params$covariate
    if (is.null(C) == FALSE) {
        if (nrow(C) < max(c(vec_idx_training, vec_idx_validation))) {
            print("Error: The covariate matrix is incompatible with the genotype data, i.e. there are less rows in the covariate matrix that the indices in vec_idx_training oridx_validation.")
            return(-1)
        }
        X_training = cbind(C[vec_idx_training, ], X_training)
        X_validation = cbind(C[vec_idx_validation, ], X_validation)
    }
    sol = glmnet::cv.glmnet(x=X_training, y=y_training, alpha=1) ### Lasso -> alpha = 1.0
    ### If the optimisation picks a lambda which does not provide estimates of the effects of the genotype data then we use the next best lambda
    n_non_zero = sol$nzero[sol$lambda==sol$lambda.min]
    if (n_non_zero == 0) {
        idx = which(sol$lambda==sol$lambda.min)
        for (i in 1:(length(sol$lambda)-1)) {
            if (sol$nzero[idx] == 0) {
                idx = idx + 1
            }
        }
        n_non_zero = sol$nzero[idx]
        b_hat = c(sol$glmnet.fit$a0[idx], sol$glmnet.fit$beta[, idx])
        y_pred = cbind(rep(1, length(vec_idx_validation)), X_validation) %*% b_hat
    } else {
        y_pred = stats::predict(sol, newx=X_validation, s="lambda.min")
    }
    colnames(y_pred) = "pred"
    df_y = merge(data.frame(id=rownames(y)[vec_idx_validation], true=y_validation), data.frame(id=rownames(y_pred), pred=y_pred), by="id")
    perf = fn_prediction_performance_metrics(y_true=df_y$true, y_pred=df_y$pred)
    perf$y_pred = y_pred
    perf$n_non_zero = n_non_zero; names(perf$n_non_zero) = NULL
    return(perf)
}

fn_elastic_net = function(G, y, vec_idx_training, vec_idx_validation, other_params=list(covariate=NULL)) {
    # n = 100; n_alleles = 2
    # G = simquantgen::fn_simulate_genotypes(n=n, l=500, ploidy=2, n_alleles=n_alleles, min_allele_freq=0.001, n_chr=5, max_pos=135e6, dist_bp_at_50perc_r2=5e6, n_threads=2)
    # y = simquantgen::fn_simulate_phenotypes(G, n_alleles=n_alleles, dist_effects="norm", n_effects=10, purely_additive=TRUE, n_networks=1, n_effects_per_network=50, h2=0.75, pheno_reps=1)$Y
    # y = scale(y, center=TRUE, scale=TRUE)
    # vec_idx_validation = c(1:round(n*0.1))
    # vec_idx_training = !(c(1:n) %in% vec_idx_validation)
    # other_params=list(covariate=NULL)
    ### Make sure that our genotype and phenotype data are labelled, i.e. have row and column names
    if (is.null(rownames(G)) | is.null(colnames(G)) | is.null(rownames(y))) {
        print("Error: Genotype and/or phenotype data have no labels. Please include the row names (entry names) and column names (loci names).")
        return(-1)
    }
    X_training = G[vec_idx_training, ]
    y_training = y[vec_idx_training]
    X_validation = G[vec_idx_validation, ]
    y_validation = y[vec_idx_validation]
    ### Adding covariate to explanatory matrix
    C = other_params$covariate
    if (is.null(C) == FALSE) {
        if (nrow(C) < max(c(vec_idx_training, vec_idx_validation))) {
            print("Error: The covariate matrix is incompatible with the genotype data, i.e. there are less rows in the covariate matrix that the indices in vec_idx_training oridx_validation.")
            return(-1)
        }
        X_training = cbind(C[vec_idx_training, ], X_training)
        X_validation = cbind(C[vec_idx_validation, ], X_validation)
    }
    sol = glmnet::cv.glmnet(x=X_training, y=y_training) ### Elastic-net -> alpha is optimised
    ### If the optimisation picks a lambda which does not provide estimates of the effects of the genotype data then we use the next best lambda
    n_non_zero = sol$nzero[sol$lambda==sol$lambda.min]
    if (n_non_zero == 0) {
        idx = which(sol$lambda==sol$lambda.min)
        for (i in 1:(length(sol$lambda)-1)) {
            if (sol$nzero[idx] == 0) {
                idx = idx + 1
            }
        }
        n_non_zero = sol$nzero[idx]
        b_hat = c(sol$glmnet.fit$a0[idx], sol$glmnet.fit$beta[, idx])
        y_pred = cbind(rep(1, length(vec_idx_validation)), X_validation) %*% b_hat
    } else {
        y_pred = stats::predict(sol, newx=X_validation, s="lambda.min")
    }
    colnames(y_pred) = "pred"
    df_y = merge(data.frame(id=rownames(y)[vec_idx_validation], true=y_validation), data.frame(id=rownames(y_pred), pred=y_pred), by="id")
    perf = fn_prediction_performance_metrics(y_true=df_y$true, y_pred=df_y$pred)
    perf$y_pred = y_pred
    perf$n_non_zero = n_non_zero; names(perf$n_non_zero) = NULL
    return(perf)
}

fn_Bayes_A = function(G, y, vec_idx_training, vec_idx_validation, other_params=list(nIter=12e3, burnIn=2e3, h2=0.5, out_prefix="bglr_bayesA-", covariate=NULL)) {
    # n = 100; n_alleles = 2
    # G = simquantgen::fn_simulate_genotypes(n=n, l=500, ploidy=2, n_alleles=n_alleles, min_allele_freq=0.001, n_chr=5, max_pos=135e6, dist_bp_at_50perc_r2=5e6, n_threads=2)
    # y = simquantgen::fn_simulate_phenotypes(G, n_alleles=n_alleles, dist_effects="norm", n_effects=10, purely_additive=TRUE, n_networks=1, n_effects_per_network=50, h2=0.75, pheno_reps=1)$Y
    # y = scale(y, center=TRUE, scale=TRUE)
    # vec_idx_validation = c(1:round(n*0.1))
    # vec_idx_training = !(c(1:n) %in% vec_idx_validation)
    # other_params=list(nIter=1000,
    #     burnIn=100,
    #     h2=0.5,
    #     out_prefix="bglr_bayesA",
    #     covariate=matrix(stats::rnorm(2*n), ncol=2))
    ### Make sure that our genotype and phenotype data are labelled, i.e. have row and column names
    if (is.null(rownames(G)) | is.null(colnames(G)) | is.null(rownames(y))) {
        print("Error: Genotype and/or phenotype data have no labels. Please include the row names (entry names) and column names (loci names).")
        return(-1)
    }
    ### Set validation set to missing
    yNA = y
    yNA[vec_idx_validation] = NA
    if (is.null(other_params$covariate)==FALSE) {
        if (nrow(other_params$covariate) != nrow(G)) {
            print("The covariate and genotype matrices are incompatible, i.e. unequal number of rows.")
            return(-1)
        }
        ETA = list(MRK=list(X=G, model="BayesA"),
                       list(X=other_params$covariate, model="FIXED"))
    } else {
        ETA = list(MRK=list(X=G, model="BayesA"))
    }
    other_params$out_prefix = gsub(" ", "-", paste(other_params$out_prefix, Sys.time(), stats::runif(1), sep="-")) ### attempt at preventing overwrites to the running Gibbs samplers in parallel
    sol = BGLR::BGLR(y=yNA, ETA=ETA, nIter=other_params$nIter, burnIn=other_params$burnIn, R2=other_params$h2, saveAt=other_params$out_prefix, verbose=FALSE)
    y_pred = sol$yHat[vec_idx_validation]
    names(y_pred) = rownames(yNA)[vec_idx_validation]
    df_y = merge(data.frame(id=rownames(y)[vec_idx_validation], true=y[vec_idx_validation]), data.frame(id=names(y_pred), pred=y_pred), by="id")
    perf = fn_prediction_performance_metrics(y_true=df_y$true, y_pred=df_y$pred)
    perf$y_pred = matrix(y_pred, ncol=1)
    colnames(perf$y_pred) = "pred"
    rownames(perf$y_pred) = names(y_pred)
    perf$n_non_zero = sum(abs(sol$ETA$MRK$b) > 0.0)
    for (f in list.files(dirname(other_params$out_prefix), pattern=basename(other_params$out_prefix))) {
        file.remove(f)    
    }
    return(perf)
}

fn_Bayes_B = function(G, y, vec_idx_training, vec_idx_validation, other_params=list(nIter=12e3, burnIn=2e3, h2=0.5, out_prefix="bglr_bayesB-", covariate=NULL)) {
    # n = 100; n_alleles = 2
    # G = simquantgen::fn_simulate_genotypes(n=n, l=500, ploidy=2, n_alleles=n_alleles, min_allele_freq=0.001, n_chr=5, max_pos=135e6, dist_bp_at_50perc_r2=5e6, n_threads=2)
    # y = simquantgen::fn_simulate_phenotypes(G, n_alleles=n_alleles, dist_effects="norm", n_effects=10, purely_additive=TRUE, n_networks=1, n_effects_per_network=50, h2=0.75, pheno_reps=1)$Y
    # y = scale(y, center=TRUE, scale=TRUE)
    # vec_idx_validation = c(1:round(n*0.1))
    # vec_idx_training = !(c(1:n) %in% vec_idx_validation)
    # other_params=list(nIter=1000,
    #     burnIn=100,
    #     h2=0.5,
    #     out_prefix="bglr_bayesB",
    #     covariate=matrix(stats::rnorm(2*n), ncol=2))
    ### Make sure that our genotype and phenotype data are labelled, i.e. have row and column names
    if (is.null(rownames(G)) | is.null(colnames(G)) | is.null(rownames(y))) {
        print("Error: Genotype and/or phenotype data have no labels. Please include the row names (entry names) and column names (loci names).")
        return(-1)
    }
    ### Set validation set to missing
    yNA = y
    yNA[vec_idx_validation] = NA
    if (is.null(other_params$covariate)==FALSE) {
        if (nrow(other_params$covariate) != nrow(G)) {
            print("The covariate and genotype matrices are incompatible, i.e. unequal number of rows.")
            return(-1)
        }
        ETA = list(MRK=list(X=G, model="BayesB"),
                       list(X=other_params$covariate, model="FIXED"))
    } else {
        ETA = list(MRK=list(X=G, model="BayesB"))
    }
    other_params$out_prefix = gsub(" ", "-", paste(other_params$out_prefix, Sys.time(), stats::runif(1), sep="-")) ### attempt at preventing overwrites to the running Gibbs samplers in parallel
    sol = BGLR::BGLR(y=yNA, ETA=ETA, nIter=other_params$nIter, burnIn=other_params$burnIn, R2=other_params$h2, saveAt=other_params$out_prefix, verbose=FALSE)
    y_pred = sol$yHat[vec_idx_validation]
    names(y_pred) = rownames(yNA)[vec_idx_validation]
    df_y = merge(data.frame(id=rownames(y)[vec_idx_validation], true=y[vec_idx_validation]), data.frame(id=names(y_pred), pred=y_pred), by="id")
    perf = fn_prediction_performance_metrics(y_true=df_y$true, y_pred=df_y$pred)
    perf$y_pred = matrix(y_pred, ncol=1)
    colnames(perf$y_pred) = "pred"
    rownames(perf$y_pred) = names(y_pred)
    perf$n_non_zero = sum(abs(sol$ETA$MRK$b) > .Machine$double.eps)
    for (f in list.files(dirname(other_params$out_prefix), pattern=basename(other_params$out_prefix))) {
        file.remove(f)    
    }
    return(perf)
}

fn_Bayes_C = function(G, y, vec_idx_training, vec_idx_validation, other_params=list(nIter=12e3, burnIn=2e3, h2=0.5, out_prefix="bglr_bayesC-", covariate=NULL)) {
    # n = 100; n_alleles = 2
    # G = simquantgen::fn_simulate_genotypes(n=n, l=500, ploidy=2, n_alleles=n_alleles, min_allele_freq=0.001, n_chr=5, max_pos=135e6, dist_bp_at_50perc_r2=5e6, n_threads=2)
    # y = simquantgen::fn_simulate_phenotypes(G, n_alleles=n_alleles, dist_effects="norm", n_effects=10, purely_additive=TRUE, n_networks=1, n_effects_per_network=50, h2=0.75, pheno_reps=1)$Y
    # y = scale(y, center=TRUE, scale=TRUE)
    # vec_idx_validation = c(1:round(n*0.1))
    # vec_idx_training = !(c(1:n) %in% vec_idx_validation)
    # other_params=list(nIter=1000,
    #     burnIn=100,
    #     h2=0.5,
    #     out_prefix="bglr_bayesC",
    #     covariate=matrix(stats::rnorm(2*n), ncol=2))
    ### Make sure that our genotype and phenotype data are labelled, i.e. have row and column names
    if (is.null(rownames(G)) | is.null(colnames(G)) | is.null(rownames(y))) {
        print("Error: Genotype and/or phenotype data have no labels. Please include the row names (entry names) and column names (loci names).")
        return(-1)
    }
    ### Set validation set to missing
    yNA = y
    yNA[vec_idx_validation] = NA
    if (is.null(other_params$covariate)==FALSE) {
        if (nrow(other_params$covariate) != nrow(G)) {
            print("The covariate and genotype matrices are incompatible, i.e. unequal number of rows.")
            return(-1)
        }
        ETA = list(MRK=list(X=G, model="BayesC"),
                       list(X=other_params$covariate, model="FIXED"))
    } else {
        ETA = list(MRK=list(X=G, model="BayesC"))
    }
    other_params$out_prefix = gsub(" ", "-", paste(other_params$out_prefix, Sys.time(), stats::runif(1), sep="-")) ### attempt at preventing overwrites to the running Gibbs samplers in parallel
    sol = BGLR::BGLR(y=yNA, ETA=ETA, nIter=other_params$nIter, burnIn=other_params$burnIn, R2=other_params$h2, saveAt=other_params$out_prefix, verbose=FALSE)
    y_pred = sol$yHat[vec_idx_validation]
    names(y_pred) = rownames(yNA)[vec_idx_validation]
    df_y = merge(data.frame(id=rownames(y)[vec_idx_validation], true=y[vec_idx_validation]), data.frame(id=names(y_pred), pred=y_pred), by="id")
    perf = fn_prediction_performance_metrics(y_true=df_y$true, y_pred=df_y$pred)
    perf$y_pred = matrix(y_pred, ncol=1)
    colnames(perf$y_pred) = "pred"
    rownames(perf$y_pred) = names(y_pred)
    perf$n_non_zero = sum(abs(sol$ETA$MRK$b) > .Machine$double.eps)
    for (f in list.files(dirname(other_params$out_prefix), pattern=basename(other_params$out_prefix))) {
        file.remove(f)    
    }
    return(perf)
}

fn_gBLUP = function(G, y, vec_idx_training, vec_idx_validation, other_params=list(covariate=NULL, mem_mb=1000)) {
    # n = 100; n_alleles = 2
    # G = simquantgen::fn_simulate_genotypes(n=n, l=500, ploidy=2, n_alleles=n_alleles, min_allele_freq=0.001, n_chr=5, max_pos=135e6, dist_bp_at_50perc_r2=5e6, n_threads=2)
    # y = simquantgen::fn_simulate_phenotypes(G, n_alleles=n_alleles, dist_effects="norm", n_effects=10, purely_additive=TRUE, n_networks=1, n_effects_per_network=50, h2=0.75, pheno_reps=1)$Y
    # y = scale(y, center=TRUE, scale=TRUE)
    # vec_idx_validation = c(1:round(n*0.1))
    # vec_idx_training = !(c(1:n) %in% vec_idx_validation)
    # other_params=list(covariate=NULL, mem_mb=1000)
    ### Make sure that our genotype and phenotype data are labelled, i.e. have row and column names
    if (is.null(rownames(G)) | is.null(colnames(G)) | is.null(rownames(y))) {
        print("Error: Genotype and/or phenotype data have no labels. Please include the row names (entry names) and column names (loci names).")
        return(-1)
    }
    if (length(rownames(y)) > length(unique(rownames(y)))) {
        print("Error: Redundant genotype IDs in the phenotype vector, i.e. redundant row names.")
        return(-2)
    }
    if (length(rownames(G)) > length(unique(rownames(G)))) {
        print("Error: Redundant genotype IDs in the genotype matrix, i.e. redundant row names.")
        return(-3)
    }
    ### Define training and validation datasets
    y_training = y
    y_training[vec_idx_validation] = NA
    y_validation = y[vec_idx_validation, 1, drop=FALSE]
    ### Generate the genetic relationship matrix from the genotype data
    A = sommer::A.mat(G)
    ### Build the phenotype and id data frame
    df_training = data.frame(y=y_training, id=as.factor(rownames(y_training))); colnames(df_training) = c("y", "id")
    ### Build the gBLUP model with/without fixed covariate/s
    ### Note that if we have more than 10 fixed covariate then the modeli fitting is more likely to fail. To prevent this we use the first 10 PCs of the fixed covariates instead
    if (is.null(other_params$covariate) == FALSE) {
        X = scale(other_params$covariate)
        if (ncol(X) > 10) {
            X = stats::prcomp(X)$x[, 1:10]
        }
        for (j in 1:ncol(X)) {
            eval(parse(text=paste0("df_training$x_", j, " = X[, j]")))
        }
        covariates_string = paste(paste0("x_", 1:ncol(X)), collapse="+")
        eval(parse(text=paste0("mod = sommer::mmer(y ~ 1 + ", covariates_string, ", random= ~vsr(id, Gu=A ), rcov= ~vsr(units), data=df_training, dateWarning=FALSE, verbose=FALSE)")))
    } else {
        mod = sommer::mmer(y ~ 1, random= ~vsr(id, Gu=A ), rcov= ~vsr(units), data=df_training, dateWarning=FALSE, verbose=FALSE)
    }
    ### Merge expected and predicted GEBVs
    df_y_validation = data.frame(id=rownames(y_validation), y_validation)
    df_y_predicted = data.frame(id=names(mod$U$`u:id`$y), pred=mod$U)
    df_y = merge(df_y_validation, df_y_predicted, by="id")
    colnames(df_y) = c("id", "true", "pred")
    perf = fn_prediction_performance_metrics(y_true=df_y$true, y_pred=df_y$pred)
    perf$y_pred = df_y$pred; names(perf$y_pred) = df_y$id
    perf$n_non_zero = NA; names(perf$n_non_zero) = NULL ### This is empty because individual alleles or markers were not used directly in the gBLUP model
    return(perf)
}

#####################################################################
############################### Tests ###############################
#####################################################################
tests_models = function() {
    source(paste0(dirname_functions, "/simquantgen/R/sim.R"))
    set.seed(123)
    n = 100; n_alleles = 2
    G = simquantgen::fn_simulate_genotypes(n=n, l=500, ploidy=2, n_alleles=n_alleles, min_allele_freq=0.001, n_chr=5, max_pos=135e6, dist_bp_at_50perc_r2=5e6, n_threads=2, verbose=TRUE)
    y = simquantgen::fn_simulate_phenotypes(G, n_alleles=n_alleles, dist_effects="norm", n_effects=10, purely_additive=TRUE, n_networks=1, n_effects_per_network=50, h2=0.5, pheno_reps=1, verbose=TRUE)$Y
    z = fn_simulate_phenotypes(G, n_alleles=n_alleles, dist_effects="norm", n_effects=10, purely_additive=FALSE, n_networks=50, n_effects_per_network=500, h2=0.5, pheno_reps=1, verbose=TRUE)$Y
    C = matrix(stats::rnorm(3*n), ncol=3)
    vec_idx_validation = c(1:round(n*0.1))
    vec_idx_training = !(c(1:n) %in% vec_idx_validation)
    test_that(
        "fn_ols", {
            print("fn_ols:")
            perf = fn_ols(G, y, vec_idx_training, vec_idx_validation)
            expect_equal(perf$corr > 0.0, TRUE)
            expect_equal(perf$n_non_zero == 500, TRUE)
            perf = fn_ols(G, y, vec_idx_training, vec_idx_validation, other_params=list(tol=.Machine$double.eps, covariate=C))
            expect_equal(perf$corr > 0.0, TRUE)
            expect_equal(perf$n_non_zero == 503, TRUE)
        }
    )

    test_that(
        "fn_ridge", {
            print("fn_ridge:")
            perf = fn_ridge(G, y, vec_idx_training, vec_idx_validation)
            expect_equal(perf$corr > 0.0, TRUE)
            expect_equal(perf$n_non_zero == 500, TRUE)
            perf = fn_ridge(G, y, vec_idx_training, vec_idx_validation, other_params=list(covariate=C))
            expect_equal(perf$corr > 0.0, TRUE)
            expect_equal(perf$n_non_zero == 503, TRUE)
        }
    )

    test_that(
        "fn_lasso", {
            print("fn_lasso:")
            perf = fn_lasso(G, y, vec_idx_training, vec_idx_validation)
            expect_equal(perf$corr > 0.1, TRUE)
            expect_equal(perf$n_non_zero < 500, TRUE)
            perf = fn_lasso(G, y, vec_idx_training, vec_idx_validation, other_params=list(covariate=C))
            expect_equal(perf$corr > 0.1, TRUE)
            expect_equal(perf$n_non_zero < 503, TRUE)
        }
    )

    test_that(
        "fn_elastic_net", {
            print("fn_elastic_net:")
            perf = fn_elastic_net(G, y, vec_idx_training, vec_idx_validation)
            expect_equal(perf$corr > 0.1, TRUE)
            expect_equal(perf$n_non_zero < 500, TRUE)
            perf = fn_elastic_net(G, y, vec_idx_training, vec_idx_validation, other_params=list(covariate=C))
            expect_equal(perf$corr > 0.1, TRUE)
            expect_equal(perf$n_non_zero < 503, TRUE)
        }
    )

    test_that(
        "fn_Bayes_A", {
            print("fn_Bayes_A:")
            perf = fn_Bayes_A(G, y, vec_idx_training, vec_idx_validation, other_params=list(nIter=1000, burnIn=100, h2=0.5, out_prefix="bglr_bayesA"))
            expect_equal(perf$corr > 0.1, TRUE)
            expect_equal(perf$n_non_zero == 500, TRUE)
            perf = fn_Bayes_A(G, y, vec_idx_training, vec_idx_validation, other_params=list(nIter=1000, burnIn=100, h2=0.5, out_prefix="bglr_bayesA", covariate=C))
            expect_equal(perf$corr > 0.1, TRUE)
            expect_equal(perf$n_non_zero == 500, TRUE) ### Excludes the fixed covariate effects
        }
    )

    test_that(
        "fn_Bayes_B", {
            print("fn_Bayes_B:")
            perf = fn_Bayes_B(G, y, vec_idx_training, vec_idx_validation, other_params=list(nIter=1000, burnIn=100, h2=0.5, out_prefix="bglr_bayesB"))
            expect_equal(perf$corr > 0.1, TRUE)
            expect_equal(perf$n_non_zero == 500, TRUE)
            perf = fn_Bayes_B(G, y, vec_idx_training, vec_idx_validation, other_params=list(nIter=1000, burnIn=100, h2=0.5, out_prefix="bglr_bayesB", covariate=C))
            expect_equal(perf$corr > 0.1, TRUE)
            expect_equal(perf$n_non_zero == 500, TRUE) ### Excludes the fixed covariate effects
        }
    )

    test_that(
        "fn_Bayes_C", {
            print("fn_Bayes_C:")
            perf = fn_Bayes_C(G, y, vec_idx_training, vec_idx_validation, other_params=list(nIter=1000, burnIn=100, h2=0.5, out_prefix="bglr_bayesC"))
            expect_equal(perf$corr > 0.1, TRUE)
            expect_equal(perf$n_non_zero == 500, TRUE)
            perf = fn_Bayes_C(G, y, vec_idx_training, vec_idx_validation, other_params=list(nIter=1000, burnIn=100, h2=0.5, out_prefix="bglr_bayesC", covariate=C))
            expect_equal(perf$corr > 0.1, TRUE)
            expect_equal(perf$n_non_zero == 500, TRUE) ### Excludes the fixed covariate effects
        }
    )

    test_that(
        "fn_gBLUP", {
            print("fn_gBLUP:")
            perf = fn_gBLUP(G, y, vec_idx_training, vec_idx_validation)
            expect_equal(perf$corr > 0.1, TRUE)
            expect_equal(is.na(perf$n_non_zero), TRUE)
            perf = fn_gBLUP(G, y, vec_idx_training, vec_idx_validation, other_params=list(covariate=C))
            expect_equal(perf$corr > 0.1, TRUE)
            expect_equal(is.na(perf$n_non_zero), TRUE) ### Excludes the fixed covariate effects
        }
    )
}

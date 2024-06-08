# source("R/io.R")
# source("R/metrics.R")
# source("R/models.R")

#' Cross-validate on a single fold, replicate, and model
#'
#' @param i index referring to the row in df_params (below) from which 
#'  the replicate id, fold number and model will be sourced from
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
#' @param df_params data frame containing all possible or just a defined non-exhaustive combinations of:
#'  - $rep: replication number
#'  - $fold: fold number
#'  - $model: model name
#' @param mat_idx_shuffle numeric n sample x r replications matrix of sample/entry/pool index shuffling where
#'  each column refer to a random shuffling of samples/entry/pool from which the identities of the
#'  training and validation sets will be sourced from
#' @param vec_set_partition_groupings vector of numeric partitioning indexes where each index refer to the
#'  fold which will serve as the validation population
#' @param prefix_tmp string referring to the prefix of the temporary files, 
#'  i.e. prefix (which can include an existing directory) of Bayesian (BGLR) model temporary files
#' @param verbose show cross-validation on a single fold, replicate, and model messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'      + $df_metrics:
#'          - $rep: replication number
#'          - $fold: fold number
#'          - $model: genomic prediction model name
#'          - $pop_training: population/s used in the training set (separated by commas if more than 1)
#'          - $pop_validation: population/s used in the validation set (separated by commas if more than 1)
#'          - $n_training: number of samples/entries/pools in the training set
#'          - $n_validation: number of samples/entries/pools in the validation set
#'          - $duration_mins: time taken in minutes to fit the genomic prediction model and assess the prediction accuracies
#'          - $n_non_zero: number of non-zero estimated effects (effects greater than machine epsilon ~2.2e-16)
#'          - $n_total_features: total number of zero and non-zero effects, i.e. total number of SNPs/alleles/features
#'          - $mbe: mean bias error
#'          - $mae: mean absolute error
#'          - $rmse: root mean squared error
#'          - $r2: coefficient of determination
#'          - $corr: Pearson's product moment correlation
#'          - $power_t10: fraction of observed top 10 phenotype values correctly predicted
#'          - $power_b10: fraction of observed bottom 10 phenotype values correctly predicted
#'          - $var_pred: variance of predicted phenotype values (estimator of additive genetic variance)
#'          - $var_true: variance of observed phenotype values (estimator of total phenotypic variance)
#'          - $h2: narrow-sense heritability estimate
#'      + $df_y_validation: 
#'          - $rep: replication number
#'          - $fold: fold number
#'          - $model: genomic prediction model name
#'          - $pop_training: population/s used in the training set (separated by commas if more than 1)
#'          - $id: names of the samples/entries/pools, 
#'          - $pop_validation: population from which the sample/entry/pool belongs to
#'          - $y_true: observed phenotype values
#'          - $y_pred: predicted phenotype values
#'      + $fname_metrics_out: filename of df_metrics saved as a the tab-delimited file with 2 rows
#'      + $fname_y_validation_out: filename of df_y_validation saved as a the tab-delimited file
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' COVAR = G %*% t(G)
#' n = nrow(G)
#' n_reps = 3
#' n_folds = 5
#' set_size = floor(n / n_folds)
#' vec_models_to_test = c("ridge", "lasso", "elastic_net", "Bayes_A", "Bayes_B", "Bayes_C", "gBLUP")
#' list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
#' df_params = expand.grid(rep=c(1:n_reps), fold=c(1:n_folds), model=vec_models_to_test)
#' mat_idx_shuffle = matrix(sample(1:n, size=n, replace=FALSE), ncol=1)
#' if (n_reps > 1) {
#'     for (r in 2:n_reps) {
#'         mat_idx_shuffle = cbind(mat_idx_shuffle, sample(1:n, size=n, replace=FALSE))
#'     }
#' }
#' vec_set_partition_groupings = rep(1:n_folds, each=set_size)
#' if (length(vec_set_partition_groupings) < n) {
#'     vec_set_partition_groupings = c(
#'          vec_set_partition_groupings, 
#'          rep(n_folds, times=(n-length(vec_set_partition_groupings)))
#'      )
#' }
#' list_cv_1 = fn_cv_1(
#'     i=2, 
#'     list_merged=list_merged, 
#'     df_params=df_params, 
#'     mat_idx_shuffle=mat_idx_shuffle, 
#'     vec_set_partition_groupings=vec_set_partition_groupings,
#'     prefix_tmp=tempfile(),
#'     verbose=TRUE)
#' @export
fn_cv_1 = function(i, list_merged, df_params, mat_idx_shuffle, vec_set_partition_groupings, prefix_tmp="gsTmp", verbose=FALSE) {
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # COVAR = G %*% t(G)
    # n = nrow(G)
    # n_reps = 3
    # n_folds = 5
    # set_size = floor(n / n_folds)
    # vec_models_to_test = c("ridge", "lasso", "elastic_net", "Bayes_A", "Bayes_B", "Bayes_C", "gBLUP")
    # list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    # i = 1
    # df_params = expand.grid(rep=c(1:n_reps), fold=c(1:n_folds), model=vec_models_to_test)
    # mat_idx_shuffle = matrix(sample(1:n, size=n, replace=FALSE), ncol=1)
    # if (n_reps > 1) {
    #     for (r in 2:n_reps) {
    #         mat_idx_shuffle = cbind(mat_idx_shuffle, sample(1:n, size=n, replace=FALSE))
    #     }
    # }
    # vec_set_partition_groupings = rep(1:n_folds, each=set_size)
    # if (length(vec_set_partition_groupings) < n) {
    #     vec_set_partition_groupings = c(vec_set_partition_groupings, rep(n_folds, times=(n-length(vec_set_partition_groupings))))
    # }
    # prefix_tmp="gsTmp"
    # verbose = TRUE
    ###################################################
    ### Input sanity check
    if (methods::is(list_merged, "gpError")) {
        error = chain(list_merged, 
            methods::new("gpError",
                code=000,
                message=paste0(
                    "Error in cross_validation::fn_cv_1(...). ",
                    "Input data (list_merged) is an error type."
                )))
        return(error)
    }
    if ((i < 1) | (i > nrow(df_params))) {
        error = methods::new("gpError",
            code=000,
            message=paste0(
                "Error in cross_validation::fn_cv_1(...). ",
                "The index (i) of df_params is beyond the number of rows in df_params (may also be less than 1)."
            ))
        return(error)
    }
    if (sum((colnames(df_params) == c("rep", "fold", "model"))) != 3) {
        error = methods::new("gpError",
            code=000,
            message=paste0(
                "Error in cross_validation::fn_cv_1(...). ",
                "The data frame of parameters is incorrect. We are expecting the following columns in order: 'rep', 'fold', and 'model'.",
                "The supplied data frame has the following columns or fields: ", paste(colnames(df_params), collapse=", ")
            ))
        return(error)
    }
    if (nrow(mat_idx_shuffle) != nrow(list_merged$G)) {
        error = methods::new("gpError",
            code=000,
            message=paste0(
                "Error in cross_validation::fn_cv_1(...). ",
                "The number of rows in the shuffling matrix (mat_idx_shuffle; ", nrow(mat_idx_shuffle), " rows) ",
                "does not match the number of samples in the input genotype and phenotype (and covariate) data (",
                nrow(list_merged$G) , " rows)."
            ))
        return(error)
    }
    if (ncol(mat_idx_shuffle) != max(df_params$rep)) {
        error = methods::new("gpError",
            code=000,
            message=paste0(
                "Error in cross_validation::fn_cv_1(...). ",
                "The number of columns in the shuffling matrix (mat_idx_shuffle; ", ncol(mat_idx_shuffle), " columns) ",
                "does not match the replications requested (", max(df_params$rep) , " replications)."
            ))
        return(error)
    }
    if (length(vec_set_partition_groupings) != nrow(list_merged$G)) {
        error = methods::new("gpError",
            code=000,
            message=paste0(
                "Error in cross_validation::fn_cv_1(...). ",
                "The number of elements in the vector of set partitioning (vec_set_partition_groupings; ", 
                length(vec_set_partition_groupings), " elements) does not match the number of samples in ",
                "the input genotype and phenotype (and covariate) data (", nrow(list_merged$G) , " rows)."
            ))
        return(error)
    }
    if (sum(range(vec_set_partition_groupings) == range(df_params$fold)) != 2) {
        error = methods::new("gpError",
            code=000,
            message=paste0(
                "Error in cross_validation::fn_cv_1(...). ",
                "The number of folds present in the vector of set partitioning (vec_set_partition_groupings; fold ", 
                min(vec_set_partition_groupings), " to fold ", max(vec_set_partition_groupings), ") ",
                "does not match the number of folds requested (fold ", min(df_params$fold), " to fold ",
                max(df_params$fold), ")."
            ))
        return(error)
    }
    ### Define prefix of intermediate output files
    if ((prefix_tmp == "") | is.na(prefix_tmp) | is.null(prefix_tmp)) {
        prefix_tmp = "gsTmp"
    }
    ### Define rep, fold, and model
    r = df_params$rep[i]
    k = df_params$fold[i]
    model = df_params$model[i]
    ### Shuffle the samples and divide into validation and training sets
    ###     as determined by the current replication and fold
    vec_idx_shuffle = mat_idx_shuffle[, r]
    vec_idx_validation = vec_idx_shuffle[vec_set_partition_groupings==k]
    vec_idx_training = vec_idx_shuffle[vec_set_partition_groupings!=k]
    ### Remove samples with missing grouping/population information
    vec_idx_validation = vec_idx_validation[!is.na(vec_idx_validation)]
    vec_idx_training = vec_idx_training[!is.na(vec_idx_training)]
    ### Define additional model input/s
    if (grepl("Bayes", model)==TRUE) {
        ### Append the input prefix into the temporary file prefix generated by Bayesian models so we don't overwite these when performing parallel computations
        other_params = list(nIter=12e3, burnIn=2e3, out_prefix=paste0(dirname(prefix_tmp), "/bglr_", model, "-", basename(prefix_tmp), "-"))
    } else {
        other_params = list(n_folds=10)
    }
    time_ini = Sys.time()
    perf = eval(parse(text=paste0("fn_", model, "(list_merged=list_merged, vec_idx_training=vec_idx_training, vec_idx_validation=vec_idx_validation, other_params=other_params, verbose=verbose)")))
    duration_mins = difftime(Sys.time(), time_ini, units="min")
    if (methods::is(perf, "gpError")) {
        error = chain(perf, methods::new("gpError",
            code=000,
            message=paste0(
                "Error in cross_validation::fn_cv_1(...). ",
                "Unable to fit the model, ", model, " and/or assess genomic prediction accuracy."
            )))
        return(error)
    }
    ### One-liner data frame of the prediction performance metrics
    df_metrics = data.frame(
        rep=r, 
        fold=k, 
        model=model,
        pop_training=paste(sort(unique(list_merged$list_pheno$pop[vec_idx_training])), collapse=","),
        pop_validation=paste(sort(unique(list_merged$list_pheno$pop[vec_idx_validation])), collapse=","),
        n_training=sum(!is.na(list_merged$list_pheno$y[vec_idx_training])),
        n_validation=sum(!is.na(list_merged$list_pheno$y[vec_idx_validation])),
        duration_mins=as.numeric(duration_mins),
        n_non_zero_effects=perf$n_non_zero,
        n_total_features=length(perf$vec_effects),
        perf$list_perf
    )
    ### Data frame of the validation phenotype values
    df_y_validation = data.frame(
        rep=r, 
        fold=k, 
        model=model,
        pop_training=paste(sort(unique(list_merged$list_pheno$pop[vec_idx_training])), collapse=","),
        perf$df_y_validation
    )
    colnames(df_y_validation)[colnames(df_y_validation)=="pop"] = "pop_validation"
    ### Note: We are not returning the allelic/loci/SNP effects during cross-validation, 
    ###     This is because this is the assessment phase,
    ###     and not the modelling phase per se, in which case we will use the entire dataset available to maximise expected accuracy.
    ### Temporary output filenames
    time_rand_id = paste0(round(as.numeric(Sys.time())), sample.int(1e6, size=1))
    fname_metrics_out = paste0(prefix_tmp, "-metrics-rep_", r, "-fold_", k, "-model_", model, "-", time_rand_id, ".tsv")
    fname_y_validation_out = paste0(prefix_tmp, "-y_pred-rep_", r, "-fold_", k, "-model_", model, "-", time_rand_id, ".tsv")
    utils::write.table(df_metrics, file=fname_metrics_out, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    utils::write.table(df_y_validation, file=fname_y_validation_out, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    gc()
    if (verbose) {
        print("Performance:")
        print(df_metrics)
        print("Output file names:")
        print(paste0("    - ", fname_metrics_out))
        print(paste0("    - ", fname_y_validation_out))
    }
    return(list(
        df_metrics=df_metrics,
        df_y_validation=df_y_validation,
        fname_metrics_out=fname_metrics_out,
        fname_y_validation_out=fname_y_validation_out
    ))
}

#' Define the type of cross-validation scheme including the models, partitioning, and shuffling,
#'  as well as estimate the maximum number of threads which can be used in parallel without memory errors.
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
#' @param cv_type (Default=1)
#'  - 1: k-fold cross-validation across all samples/entries/pools regardless of their groupings
#'  - 2: pairwise-population cross-validation, e.g. training on population A and validation on population B
#'  - 3: leave-one-population-out cross-validation, e.g. training on populations 1 to 9 and validation on population 10
#' @param n_folds if cv_type=1, i.e. k-fold cross-validation is desired, then this defines the number of partitions or folds.
#'  However, if the resulting size per fold is less than 2, then this return an error.
#'  Additionally, this parameter is ignored is cv_type=2 or cv_type=3 (Default=10).
#' @param n_reps if cv_type=1, i.e. k-fold cross-validation is desired, then this defines the number of random shuffling
#'  for repeated k-fold cross-validation.
#'  Similar to n_folds, this parameter is ignored is cv_type=2 or cv_type=3 (Default=10).
#' @param vec_models_to_test genomic prediction models to use (uses all the models below by default). Please choose one or more from:
#'  - "ridge": \eqn{Cost_{ridge} = \Sigma(y - X\beta)^2 + \lambda\Sigma\beta^2$, where $\hat{\beta} = {(X^TX + \lambda I)^{-1} X^Ty}}
#'      https://en.wikipedia.org/wiki/Ridge_regression
#'  - "lasso": \eqn{Cost_{lasso} = \Sigma(y - X\beta)^2 + \lambda\Sigma|\beta|}
#'      https://en.wikipedia.org/wiki/Lasso_(statistics)
#'  - "elastic_net": \eqn{Cost_{ridge} = \Sigma(y - X\beta)^2 + (1-\alpha)\lambda\Sigma\beta^2 + \alpha\lambda\Sigma|\beta|}
#'      https://en.wikipedia.org/wiki/Elastic_net_regularization
#'  - "Bayes_A": scaled t-distributed effects
#'      https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf
#'  - "Bayes_B": scaled t-distributed effects with probability \eqn{\pi}; and zero effects with probability \eqn{1-\pi}, 
#'      where \eqn{\pi \sim \beta(\theta_1, \theta_2)}.
#'      https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf
#'  - "Bayes_C": normally distributed effects (\eqn{N(0, \sigma^2_{\beta})}) with probability \eqn{\pi}; 
#'      and zero effects with probability \eqn{1-\pi}, where \eqn{\pi \sim \beta(\theta_1, \theta_2)}.
#'      https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf
#'  - "gBLUP": genotype best linear unbiased prediction (gBLUP) using genomic relationship matrix to predict missing breeding values
#'      via Direct-Inversion Newton-Raphson or Average Information 
#'      using the sommer R package (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4894563/).
#'      https://link.springer.com/protocol/10.1007/978-1-62703-447-0_13
#' @param max_mem_Gb maximum memory in gigabytes available for computation (Default=15)
#' @param verbose show cross-validation parameter preparation messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'  - $df_params: data frame containing all possible or just a defined non-exhaustive combinations of:
#'      + $rep: replication number
#'      + $fold: fold number
#'      + $model: model name
#'  - $mat_idx_shuffle numeric: n sample x r replications matrix of sample/entry/pool index shuffling where
#'      each column refer to a random shuffling of samples/entry/pool from which the identities of the
#'      training and validation sets will be sourced from
#'  - $vec_set_partition_groupings vector of numeric partitioning indexes where each index refer to the
#'      fold which will serve as the validation population
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' COVAR = G %*% t(G)
#' list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
#' list_cv_params = fn_cross_validation_preparation(list_merged, cv_type=1)
#' list_cv_params = fn_cross_validation_preparation(list_merged, cv_type=3)
#' list_cv_params = fn_cross_validation_preparation(list_merged, cv_type=2)
#' list_merged$list_pheno$pop = rep(c("popA", "popB"), times=length(list_merged$list_pheno$pop)/2)
#' list_cv_params = fn_cross_validation_preparation(list_merged, cv_type=2)
#' @export
fn_cross_validation_preparation = function(list_merged, cv_type=1, n_folds=10, n_reps=10, 
    vec_models_to_test=c("ridge","lasso","elastic_net","Bayes_A","Bayes_B","Bayes_C","gBLUP"),
    max_mem_Gb=15, verbose=FALSE)
{
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # COVAR = G %*% t(G)
    # list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    # cv_type=1
    # n_folds = 10
    # n_reps = 10
    # vec_models_to_test = c("ridge", "lasso", "elastic_net", "Bayes_A", "Bayes_B", "Bayes_C", "gBLUP")
    # max_mem_Gb = 15
    # verbose = TRUE
    ###################################################
    if (verbose) {print("Preparing cross-validation parameters.")}
    ### Input sanity check
    if (methods::is(list_merged, "gpError")) {
        error = chain(list_merged, 
            methods::new("gpError",
                code=000,
                message=paste0(
                    "Error in cross_validation::fn_cross_validation_preparation(...). ",
                    "Input data (list_merged) is an error type."
                )))
        return(error)
    }
    if (cv_type == 1) {
        ###############################
        ### k-fold cross-validation ###
        ###############################
        n = nrow(list_merged$G)
        ### Adjust the number of folds so that we have at least 2 samples/entries/pools in the validation set
        set_size = floor(n / n_folds)
        while ((set_size < 2) & (n_folds > 1)) {
            if (verbose) {print(paste0("Reducing the number of folds from ", n_folds, " to ", n_folds-1, "."))}
            n_folds = n_folds - 1
        }
        if (set_size < 2) {
            error = methods::new("gpError",
                code=000,
                messages=paste0(
                    "Error in cross_validation::fn_cross_validation_preparation(...). ",
                    "The size of the data set is too small, n= ", n, "."
                ))
            return(error)
        }
        vec_set_partition_groupings = rep(1:n_folds, each=set_size)
        if (length(vec_set_partition_groupings) < n) {
            vec_set_partition_groupings = c(vec_set_partition_groupings, rep(n_folds, times=(n-length(vec_set_partition_groupings))))
        }
        ### Prepare shuffling across replications
        mat_idx_shuffle = matrix(sample(1:n, size=n, replace=FALSE), ncol=1)
        if (n_reps > 1) {
            for (r in 2:n_reps) {
                mat_idx_shuffle = cbind(mat_idx_shuffle, sample(1:n, size=n, replace=FALSE))
            }
        }
        ### Prepare matrix of rep x fold x model combinations and sort by fold so that we can start comparing as soon as results become available
        df_params = expand.grid(rep=c(1:n_reps), fold=c(1:n_folds), model=vec_models_to_test)
        df_params = df_params[order(df_params$rep), ]
        df_params = df_params[order(df_params$fold), ]
        ### Estimate the maximum number of threads which can be used without running out of memory
        list_mem = fn_estimate_memory_footprint(
            X=list_merged,
            n_models=length(vec_models_to_test),
            n_folds=n_folds,
            n_reps=n_reps, 
            memory_requested_Gb=max_mem_Gb,
            verbose=verbose)
    } else if (cv_type == 2) {
        ############################################
        ### Pairwise-population cross-validation ###
        ############################################
        n = nrow(list_merged$G)
        ### Define the partitioning as the population grouping
        vec_set_partition_groupings = as.numeric(as.factor(list_merged$list_pheno$pop))
        n_folds = length(unique(list_merged$list_pheno$pop))
        if (n_folds != 2) {
            error = methods::new("gpError",
                code=000,
                message=paste0(
                    "Error in cross_validation::fn_cross_validation_preparation(...). ",
                    "Cannot perform pairwise-population cross-validation (cv_type=2) ",
                    "because the number of populations (", n_folds, " populations) in the data set is not equal to 2."
                ))
            return(error)
        }
        ### No shuffling needed as cross-validation is not replicated
        mat_idx_shuffle = matrix(1:n, ncol=1)
        ### Prepare matrix of rep x fold x model combinations
        df_params = expand.grid(rep=c(1), fold=c(1:n_folds), model=vec_models_to_test)
        ### Estimate the maximum number of threads which can be used without running out of memory
        list_mem = fn_estimate_memory_footprint(
            X=list_merged,
            n_models=length(vec_models_to_test),
            n_folds=n_folds,
            n_reps=1, 
            memory_requested_Gb=max_mem_Gb,
            verbose=verbose)
    } else if (cv_type == 3) {
        #################################################
        ### Leave-one-population-out cross-validation ###
        #################################################
        n = nrow(list_merged$G)
        ### Define the partitioning as the population grouping
        vec_set_partition_groupings = as.numeric(as.factor(list_merged$list_pheno$pop))
        n_folds = length(unique(list_merged$list_pheno$pop))
        if (n_folds == 1) {
            error = methods::new("gpError",
                code=000,
                message=paste0(
                    "Error in cross_validation::fn_cross_validation_preparation(...). ",
                    "Cannot perform leave-one-population-out cross-validation (cv_type=3) ",
                    "because there is only one population in the data set."
                ))
            return(error)
        }
        ### No shuffling needed as cross-validation is not replicated
        mat_idx_shuffle = matrix(1:n, ncol=1)
        ### Prepare matrix of rep x fold x model combinations
        df_params = expand.grid(rep=c(1), fold=c(1:n_folds), model=vec_models_to_test)
        ### Estimate the maximum number of threads which can be used without running out of memory
        list_mem = fn_estimate_memory_footprint(
            X=list_merged,
            n_models=length(vec_models_to_test),
            n_folds=n_folds,
            n_reps=1, 
            memory_requested_Gb=max_mem_Gb,
            verbose=verbose)
    } else {
        error = methods::new("gpError",
            code=000,
            messages=paste0(
                "Error in cross_validation::fn_cross_validation_preparation(...). ",
                "The cross-validation type, cv_type=", cv_type, " is invalid. ",
                "Please choose: ",
                "  --> '1' for k-fold cross-validation across all samples/entries/pools regardless of their groupings. ",
                "  --> '2' for pairwise-population cross-validation, e.g. training on population A and validation on population B. ",
                "  --> '3' for leave-one-population-out cross-validation, e.g. training on populations 1 to 9 and validation on population 10."
            ))
        return(error)
    }
    ### Memory allocation error handling
    if (methods::is(list_mem, "gpError")) {
        error = chain(list_mem, methods::new("gpError",
            code=000,
            message=paste0(
                "Error in cross_validation::fn_cross_validation_preparation(...). ",
                "Failed to estimate memory allocation requirements for parallel computations ",
                "and the maximum number of threads which can be used to avoid out-of-memory (OOM) error."
            )))
        return(error)
    }
    ### Print the full list of cross-validation sets, replications and models combinations
    if (verbose) {
        print("Full list of cross-validation sets, replications and models combinations:")
        print(df_params)
    }
    ### Output
    return(list(
        df_params=df_params,
        mat_idx_shuffle=mat_idx_shuffle,
        vec_set_partition_groupings=vec_set_partition_groupings,
        list_mem=list_mem
    ))
}

#' K-fold cross-validation per population
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
#' @param n_folds if cv_type=1, i.e. k-fold cross-validation is desired, then this defines the number of partitions or folds.
#'  However, if the resulting size per fold is less than 2, then this return an error.
#'  Additionally, this parameter is ignored is cv_type=2 or cv_type=3 (Default=10).
#' @param n_reps if cv_type=1, i.e. k-fold cross-validation is desired, then this defines the number of random shuffling
#'  for repeated k-fold cross-validation.
#'  Similar to n_folds, this parameter is ignored is cv_type=2 or cv_type=3 (Default=10).
#' @param vec_models_to_test genomic prediction models to use (uses all the models below by default). Please choose one or more from:
#'  - "ridge": \eqn{Cost_{ridge} = \Sigma(y - X\beta)^2 + \lambda\Sigma\beta^2$, where $\hat{\beta} = {(X^TX + \lambda I)^{-1} X^Ty}}
#'      https://en.wikipedia.org/wiki/Ridge_regression
#'  - "lasso": \eqn{Cost_{lasso} = \Sigma(y - X\beta)^2 + \lambda\Sigma|\beta|}
#'      https://en.wikipedia.org/wiki/Lasso_(statistics)
#'  - "elastic_net": \eqn{Cost_{ridge} = \Sigma(y - X\beta)^2 + (1-\alpha)\lambda\Sigma\beta^2 + \alpha\lambda\Sigma|\beta|}
#'      https://en.wikipedia.org/wiki/Elastic_net_regularization
#'  - "Bayes_A": scaled t-distributed effects
#'      https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf
#'  - "Bayes_B": scaled t-distributed effects with probability \eqn{\pi}; and zero effects with probability \eqn{1-\pi}, 
#'      where \eqn{\pi \sim \beta(\theta_1, \theta_2)}.
#'      https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf
#'  - "Bayes_C": normally distributed effects (\eqn{N(0, \sigma^2_{\beta})}) with probability \eqn{\pi}; 
#'      and zero effects with probability \eqn{1-\pi}, where \eqn{\pi \sim \beta(\theta_1, \theta_2)}.
#'      https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf
#'  - "gBLUP": genotype best linear unbiased prediction (gBLUP) using genomic relationship matrix to predict missing breeding values
#'      via Direct-Inversion Newton-Raphson or Average Information 
#'      using the sommer R package (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4894563/).
#'      https://link.springer.com/protocol/10.1007/978-1-62703-447-0_13
#' @param bool_parallel perform multi-threaded cross-validation (Default=TRUE)
#' @param max_mem_Gb maximum memory in gigabytes available for computation (Default=15)
#' @param n_threads total number of computing threads available (Default=2)
#' @param dir_output output directory where temporary text and Rds files will be saved into (Default=NULL)
#' @param verbose show within population cross-validation messages? (Default=FALSE)
#' @returns
#'  - Ok: filename of temporary Rds file containing a 2-element list:
#'      + $METRICS_WITHIN_POP (row-binded df_metrics across populations, reps, folds, and models):
#'          - $rep: replication number
#'          - $fold: fold number
#'          - $model: genomic prediction model name
#'          - $pop_training: population/s used in the training set (separated by commas if more than 1)
#'          - $pop_validation: population/s used in the validation set (separated by commas if more than 1)
#'          - $n_training: number of samples/entries/pools in the training set
#'          - $n_validation: number of samples/entries/pools in the validation set
#'          - $duration_mins: time taken in minutes to fit the genomic prediction model and assess the prediction accuracies
#'          - $n_non_zero: number of non-zero estimated effects (effects greater than machine epsilon ~2.2e-16)
#'          - $n_total_features: total number of zero and non-zero effects, i.e. total number of SNPs/alleles/features
#'          - $mbe: mean bias error
#'          - $mae: mean absolute error
#'          - $rmse: root mean squared error
#'          - $r2: coefficient of determination
#'          - $corr: Pearson's product moment correlation
#'          - $power_t10: fraction of observed top 10 phenotype values correctly predicted
#'          - $power_b10: fraction of observed bottom 10 phenotype values correctly predicted
#'          - $var_pred: variance of predicted phenotype values (estimator of additive genetic variance)
#'          - $var_true: variance of observed phenotype values (estimator of total phenotypic variance)
#'          - $h2: narrow-sense heritability estimate
#'      + $YPRED_WITHIN_POP (row-binded df_y_validation across populations, reps, folds, and models):
#'          - $rep: replication number
#'          - $fold: fold number
#'          - $model: genomic prediction model name
#'          - $pop_training: population/s used in the training set (separated by commas if more than 1)
#'          - $id: names of the samples/entries/pools, 
#'          - $pop_validation: population from which the sample/entry/pool belongs to
#'          - $y_true: observed phenotype values
#'          - $y_pred: predicted phenotype values
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' COVAR = G %*% t(G)
#' list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
#' fname_within_Rds = fn_cross_validation_within_population(list_merged, n_folds=2, n_reps=1, 
#'  vec_models_to_test=c("ridge","lasso"), verbose=TRUE)
#' list_within = readRDS(fname_within_Rds)
#' head(list_within$METRICS_WITHIN_POP)
#' head(list_within$YPRED_WITHIN_POP)
#' @export
fn_cross_validation_within_population = function(list_merged, n_folds=10, n_reps=10, 
    vec_models_to_test=c("ridge","lasso","elastic_net","Bayes_A","Bayes_B","Bayes_C","gBLUP"),
    bool_parallel=TRUE, max_mem_Gb=15, n_threads=2, dir_output=NULL, verbose=FALSE)
{
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # COVAR = G %*% t(G)
    # list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    # n_folds = 2
    # n_reps = 1
    # vec_models_to_test = c("ridge", "lasso")
    # max_mem_Gb = 15
    # n_threads = 2
    # dir_output = NULL
    # bool_parallel = TRUE
    # verbose = TRUE
    ###################################################
    if (verbose) {
        print("##########################################")
        print("### Within population cross-validation ###")
        print("##########################################")
    }
    ### Input sanity check
    if (methods::is(list_merged, "gpError")) {
        error = chain(list_merged, 
            methods::new("gpError",
                code=000,
                message=paste0(
                    "Error in cross_validation::fn_cross_validation_within_population(...). ",
                    "Input data (list_merged) is an error type."
                )))
        return(error)
    }
    ### Define the output directory
    if (!is.null(dir_output)) {
        if (!dir.exists(dir_output)) {
            dir.create(dir_output, showWarnings=FALSE)
        }
        
    } else {
        dir_output = tempfile()
        dir.create(dir_output, showWarnings=FALSE)
    }
    if (!dir.exists(dir_output)) {
        error = methods::new("gpError",
            code=000,
            message=paste0(
                "Error in cross_validation::fn_cross_validation_within_population(...). ",
                "Unable to create the output directory: ", dir_output, ". ",
                "Please check your permissions to write into that directory."
            ))
        return(error)
    }
    ### Determine the number of populations
    vec_populations = sort(unique(list_merged$list_pheno$pop))
    ### Instantiate the vector of Rds filenames containing the temporary output data per population
    vec_fname_within_Rds = c()
    for (population in vec_populations) {
        # population = vec_populations[1]
        ### Subset the data per population
        vec_idx = which(list_merged$list_pheno$pop == population)
        list_merged_sub = fn_subset_merged_genotype_and_phenotype(list_merged=list_merged, vec_idx=vec_idx, verbose=verbose)
        if (methods::is(list_merged_sub, "gpError")) {
            error = chain(list_merged_sub, methods::new("gpError",
                code=000,
                message=paste0(
                    "Error in cross_validation::fn_cross_validation_within_population(...). ",
                    "Failed to subset the data set."
                )))
            return(error)
        }
        ### Define the cross-validation parameters as well as the maximum number of threads we can safely use in parallel
        list_cv_params = fn_cross_validation_preparation(
            list_merged=list_merged_sub,
            cv_type=1,
            n_folds=n_folds,
            n_reps=n_reps, 
            vec_models_to_test=vec_models_to_test,
            max_mem_Gb=max_mem_Gb,
            verbose=verbose)
        if (methods::is(list_cv_params, "gpError")) {
            error = chain(list_cv_params, methods::new("gpError",
                code=000,
                message=paste0(
                    "Error in cross_validation::fn_cross_validation_within_population(...). ",
                    "Failed to define the cross-validation parameters."
                )))
            return(error)
        }
        if (list_cv_params$list_mem$n_threads <= 1) {
            if (verbose) {
                print("Data set is too large for parallel computations.")
                print("Performing iterative cross-validation.")
            }
            bool_parallel = FALSE
        }
        ### Cross-validation
        if (bool_parallel) {
            ###########################################
            ### Multi-threaded cross-validation
            list_list_perf = parallel::mclapply(c(1:nrow(list_cv_params$df_params)), 
                FUN=fn_cv_1, 
                    list_merged=list_merged_sub,
                    df_params=list_cv_params$df_params,
                    mat_idx_shuffle=list_cv_params$mat_idx_shuffle,
                    vec_set_partition_groupings=list_cv_params$vec_set_partition_groupings,
                    prefix_tmp=file.path(dir_output, paste0("within-", population)),
                    verbose=verbose, 
                mc.cores=min(c(n_threads, list_cv_params$n_threads)))
            for (idx in 1:length(list_list_perf)) {
                if (methods::is(list_list_perf[[idx]], "gpError")) {
                    error = chain(list_list_perf[[idx]], methods::new("gpError",
                        code=000,
                        message=paste0(
                            "Error in cross_validation::fn_cross_validation_within_population(...). ",
                            "Something went wrong in the execution of multi-threaded within population k-fold cross-validation. ",
                            "Please check re-run cross_validation::fn_cross_validation_within_population(...) with ",
                            "bool_parallel=FALSE to identify the error."
                        )))
                    return(error)
                }
            }
        } else {
            ############################################
            ### Single-threaded cross-validation
            list_list_perf = list()
            for (i in 1:nrow(list_cv_params$df_params)) {
                list_perf = fn_cv_1(
                    i=i,
                    list_merged=list_merged_sub,
                    df_params=list_cv_params$df_params,
                    mat_idx_shuffle=list_cv_params$mat_idx_shuffle,
                    vec_set_partition_groupings=list_cv_params$vec_set_partition_groupings,
                    prefix_tmp=file.path(dir_output, paste0("within-", population)),
                    verbose=verbose
                )
                if (methods::is(list_perf, "gpError")) {
                    error = chain(list_perf, methods::new("gpError",
                        code=000,
                        message=paste0(
                            "Error in cross_validation::fn_cross_validation_within_population(...). ",
                            "Error running cross-validation for population: ", population, " at ",
                            "rep: ", list_cv_params$df_params$rep[i], ", ",
                            "fold: ", list_cv_params$df_params$fold[i], ", and ",
                            "model: ", list_cv_params$df_params$model[i], "."
                        )))
                    return(error)
                }
                eval(parse(text=paste0("list_list_perf$`", i, "` = list_perf")))
            }
        }
        ### Concatenate performances
        df_metrics = NULL
        df_y_validation = NULL
        for (list_perf in list_list_perf) {
            # list_perf = list_list_perf[[1]]
            if (is.null(df_metrics) & is.null(df_y_validation)) {
                df_metrics = list_perf$df_metrics
                df_y_validation = list_perf$df_y_validation
            } else {
                df_metrics = rbind(df_metrics, list_perf$df_metrics)
                df_y_validation = rbind(df_y_validation, list_perf$df_y_validation)
            }
            ### Clean-up to reduce memory footprint
            list_perf$df_metrics = NULL
            list_perf$df_y_validation = NULL
            gc()
        }
        ### Save temporary Rds output per population
        time_rand_id = paste0(round(as.numeric(Sys.time())), sample.int(1e6, size=1))
        fname_within_Rds = file.path(dir_output, paste0("within-tempout-", population, "-", time_rand_id, ".Rds"))
        saveRDS(list(
            df_metrics=df_metrics,
            df_y_validation=df_y_validation), 
            file=fname_within_Rds)
        ### Take note of the names of these temporary Rds output files
        vec_fname_within_Rds = c(vec_fname_within_Rds, fname_within_Rds)
        ### Clean-up temporary files generated by fn_cv_1 in parallel
        for (list_perf in list_list_perf) {
            unlink(list_perf$fname_metrics_out)
            unlink(list_perf$fname_y_validation_out)
        }
        ### Clean-up the memory used by subsetting the data set
        list_merged_sub = NULL 
        list_list_perf = NULL 
        list_perf = NULL 
        gc()
    }
    ### Concatenate the temporary output Rds files
    METRICS_WITHIN_POP = NULL
    YPRED_WITHIN_POP = NULL
    for (fname_within_Rds in vec_fname_within_Rds) {
        # fname_within_Rds = vec_fname_within_Rds[1]
        list_tmp = readRDS(fname_within_Rds)
        if (is.null(METRICS_WITHIN_POP) & is.null(YPRED_WITHIN_POP)) {
            METRICS_WITHIN_POP = list_tmp$df_metrics
            YPRED_WITHIN_POP = list_tmp$df_y_validation
        } else {
            METRICS_WITHIN_POP = rbind(METRICS_WITHIN_POP, list_tmp$df_metrics)
            YPRED_WITHIN_POP = rbind(YPRED_WITHIN_POP, list_tmp$df_y_validation)
        }
    }
    ### Save the concatenated output data frames across populations and clean-up
    time_rand_id = paste0(round(as.numeric(Sys.time())), sample.int(1e6, size=1))
    fname_within_Rds = file.path(dir_output, paste0("within-tempout-", time_rand_id, ".Rds"))
    saveRDS(list(
        METRICS_WITHIN_POP=METRICS_WITHIN_POP,
        YPRED_WITHIN_POP=YPRED_WITHIN_POP),
        file=fname_within_Rds)
    ### Clean-up
    for (fname in vec_fname_within_Rds) {
        unlink(fname)
    }
    rm("METRICS_WITHIN_POP")
    rm("YPRED_WITHIN_POP")
    gc()
    ### Output
    return(fname_within_Rds)
}

#' K-fold cross-validation across populations without accounting for population grouping
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
#' @param n_folds if cv_type=1, i.e. k-fold cross-validation is desired, then this defines the number of partitions or folds.
#'  However, if the resulting size per fold is less than 2, then this return an error.
#'  Additionally, this parameter is ignored is cv_type=2 or cv_type=3 (Default=10).
#' @param n_reps if cv_type=1, i.e. k-fold cross-validation is desired, then this defines the number of random shuffling
#'  for repeated k-fold cross-validation.
#'  Similar to n_folds, this parameter is ignored is cv_type=2 or cv_type=3 (Default=10).
#' @param vec_models_to_test genomic prediction models to use (uses all the models below by default). Please choose one or more from:
#'  - "ridge": \eqn{Cost_{ridge} = \Sigma(y - X\beta)^2 + \lambda\Sigma\beta^2$, where $\hat{\beta} = {(X^TX + \lambda I)^{-1} X^Ty}}
#'      https://en.wikipedia.org/wiki/Ridge_regression
#'  - "lasso": \eqn{Cost_{lasso} = \Sigma(y - X\beta)^2 + \lambda\Sigma|\beta|}
#'      https://en.wikipedia.org/wiki/Lasso_(statistics)
#'  - "elastic_net": \eqn{Cost_{ridge} = \Sigma(y - X\beta)^2 + (1-\alpha)\lambda\Sigma\beta^2 + \alpha\lambda\Sigma|\beta|}
#'      https://en.wikipedia.org/wiki/Elastic_net_regularization
#'  - "Bayes_A": scaled t-distributed effects
#'      https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf
#'  - "Bayes_B": scaled t-distributed effects with probability \eqn{\pi}; and zero effects with probability \eqn{1-\pi}, 
#'      where \eqn{\pi \sim \beta(\theta_1, \theta_2)}.
#'      https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf
#'  - "Bayes_C": normally distributed effects (\eqn{N(0, \sigma^2_{\beta})}) with probability \eqn{\pi}; 
#'      and zero effects with probability \eqn{1-\pi}, where \eqn{\pi \sim \beta(\theta_1, \theta_2)}.
#'      https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf
#'  - "gBLUP": genotype best linear unbiased prediction (gBLUP) using genomic relationship matrix to predict missing breeding values
#'      via Direct-Inversion Newton-Raphson or Average Information 
#'      using the sommer R package (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4894563/).
#'      https://link.springer.com/protocol/10.1007/978-1-62703-447-0_13
#' @param bool_parallel perform multi-threaded cross-validation (Default=TRUE)
#' @param max_mem_Gb maximum memory in gigabytes available for computation (Default=15)
#' @param n_threads total number of computing threads available (Default=2)
#' @param dir_output output directory where temporary text and Rds files will be saved into (Default=NULL)
#' @param verbose show bulk across population cross-validation messages? (Default=FALSE)
#' @returns
#'  - Ok: filename of temporary Rds file containing a 2-element list:
#'      + $METRICS_ACROSS_POP_BULK:
#'          - $rep: replication number
#'          - $fold: fold number
#'          - $model: genomic prediction model name
#'          - $pop_training: population/s used in the training set (separated by commas if more than 1)
#'          - $pop_validation: population/s used in the validation set (separated by commas if more than 1)
#'          - $n_training: number of samples/entries/pools in the training set
#'          - $n_validation: number of samples/entries/pools in the validation set
#'          - $duration_mins: time taken in minutes to fit the genomic prediction model and assess the prediction accuracies
#'          - $n_non_zero: number of non-zero estimated effects (effects greater than machine epsilon ~2.2e-16)
#'          - $n_total_features: total number of zero and non-zero effects, i.e. total number of SNPs/alleles/features
#'          - $mbe: mean bias error
#'          - $mae: mean absolute error
#'          - $rmse: root mean squared error
#'          - $r2: coefficient of determination
#'          - $corr: Pearson's product moment correlation
#'          - $power_t10: fraction of observed top 10 phenotype values correctly predicted
#'          - $power_b10: fraction of observed bottom 10 phenotype values correctly predicted
#'          - $var_pred: variance of predicted phenotype values (estimator of additive genetic variance)
#'          - $var_true: variance of observed phenotype values (estimator of total phenotypic variance)
#'          - $h2: narrow-sense heritability estimate
#'      + $YPRED_ACROSS_POP_BULK:
#'          - $rep: replication number
#'          - $fold: fold number
#'          - $model: genomic prediction model name
#'          - $pop_training: population/s used in the training set (separated by commas if more than 1)
#'          - $id: names of the samples/entries/pools, 
#'          - $pop_validation: population from which the sample/entry/pool belongs to
#'          - $y_true: observed phenotype values
#'          - $y_pred: predicted phenotype values
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' COVAR = G %*% t(G)
#' list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
#' fname_across_bulk_Rds = fn_cross_validation_across_populations_bulk(list_merged, n_folds=2, 
#'  n_reps=1, vec_models_to_test=c("ridge","lasso"), verbose=TRUE)
#' list_across_bulk = readRDS(fname_across_bulk_Rds)
#' head(list_across_bulk$METRICS_ACROSS_POP_BULK)
#' head(list_across_bulk$YPRED_ACROSS_POP_BULK)
#' @export
fn_cross_validation_across_populations_bulk = function(list_merged, n_folds=10, n_reps=10, 
    vec_models_to_test=c("ridge","lasso","elastic_net","Bayes_A","Bayes_B","Bayes_C","gBLUP"),
    bool_parallel=TRUE, max_mem_Gb=15, n_threads=2, dir_output=NULL, verbose=FALSE)
{
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # COVAR = G %*% t(G)
    # list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    # n_folds = 2
    # n_reps = 1
    # vec_models_to_test = c("ridge", "lasso")
    # max_mem_Gb = 15
    # n_threads = 2
    # dir_output = NULL
    # bool_parallel = TRUE
    # verbose = TRUE
    ###################################################
    if (verbose) {
        print("###################################################")
        print("### Across populations cross-validation in bulk ###")
        print("###################################################")
    }
    ### Input sanity check
    if (methods::is(list_merged, "gpError")) {
        error = chain(list_merged, 
            methods::new("gpError",
                code=000,
                message=paste0(
                    "Error in cross_validation::fn_cross_validation_across_populations_bulk(...). ",
                    "Input data (list_merged) is an error type."
                )))
        return(error)
    }
    ### Define the output directory
    if (!is.null(dir_output)) {
        if (!dir.exists(dir_output)) {
            dir.create(dir_output, showWarnings=FALSE)
        }
        
    } else {
        dir_output = tempfile()
        dir.create(dir_output, showWarnings=FALSE)
    }
    if (!dir.exists(dir_output)) {
        error = methods::new("gpError",
            code=000,
            message=paste0(
                "Error in cross_validation::fn_cross_validation_across_populations_bulk(...). ",
                "Unable to create the output directory: ", dir_output, ". ",
                "Please check your permissions to write into that directory."
            ))
        return(error)
    }
    ### Check if we have more than 1 population
    vec_populations = sort(unique(list_merged$list_pheno$pop))
    if (length(vec_populations) == 1) {
        error = methods::new("gpError",
            code=000,
            message=paste0(
                "Error in cross_validation::fn_cross_validation_across_populations_bulk(...). ",
                "Cannot perform bulked across populations cross-validation ",
                "because there is only 1 population in the data set."
            ))
        return(error)
    }
    ### Define the cross-validation parameters as well as the maximum number of threads we can safely use in parallel
    list_cv_params = fn_cross_validation_preparation(
        list_merged=list_merged,
        cv_type=1,
        n_folds=n_folds,
        n_reps=n_reps, 
        vec_models_to_test=vec_models_to_test,
        max_mem_Gb=max_mem_Gb,
        verbose=verbose)
    if (methods::is(list_cv_params, "gpError")) {
        error = chain(list_cv_params, methods::new("gpError",
            code=000,
            message=paste0(
                "Error in cross_validation::fn_cross_validation_within_population(...). ",
                "Failed to define the cross-validation parameters."
            )))
        return(error)
    }
    if (list_cv_params$list_mem$n_threads <= 1) {
        if (verbose) {
            print("Data set is too large for parallel computations.")
            print("Performing iterative cross-validation.")
        }
        bool_parallel = FALSE
    }
    ### Cross-validation
    if (bool_parallel) {
        ###########################################
        ### Multi-threaded cross-validation
        list_list_perf = parallel::mclapply(c(1:nrow(list_cv_params$df_params)), 
            FUN=fn_cv_1, 
                list_merged=list_merged,
                df_params=list_cv_params$df_params,
                mat_idx_shuffle=list_cv_params$mat_idx_shuffle,
                vec_set_partition_groupings=list_cv_params$vec_set_partition_groupings,
                prefix_tmp=file.path(dir_output, "across_bulk"),
                verbose=verbose, 
            mc.cores=min(c(n_threads, list_cv_params$n_threads)))
        for (idx in 1:length(list_list_perf)) {
            if (methods::is(list_list_perf[[idx]], "gpError")) {
                error = chain(list_list_perf[[idx]], methods::new("gpError",
                    code=000,
                    message=paste0(
                        "Error in cross_validation::fn_cross_validation_across_populations_bulk(...). ",
                        "Something went wrong in the execution of multi-threaded across population cross-validation, ",
                        "without regard for population groupings. ",
                        "Please check re-run cross_validation::fn_cross_validation_across_populations_bulk(...) with ",
                        "bool_parallel=FALSE to identify the error."
                    )))
                return(error)
            }
        }
    } else {
        ############################################
        ### Single-threaded cross-validation
        list_list_perf = list()
        for (i in 1:nrow(list_cv_params$df_params)) {
            # i = 1
            list_perf = fn_cv_1(
                i=i,
                list_merged=list_merged,
                df_params=list_cv_params$df_params,
                mat_idx_shuffle=list_cv_params$mat_idx_shuffle,
                vec_set_partition_groupings=list_cv_params$vec_set_partition_groupings,
                prefix_tmp=file.path(dir_output, "across_bulk"),
                verbose=verbose
            )
            if (methods::is(list_perf, "gpError")) {
                error = chain(list_perf, methods::new("gpError",
                    code=000,
                    message=paste0(
                        "Error in cross_validation::fn_cross_validation_across_populations_bulk(...). ",
                        "Error running bulked across population cross-validation at ",
                        "rep: ", list_cv_params$df_params$rep[i], ", ",
                        "fold: ", list_cv_params$df_params$fold[i], ", and ",
                        "model: ", list_cv_params$df_params$model[i], "."
                    )))
                return(error)
            }
            eval(parse(text=paste0("list_list_perf$`", i, "` = list_perf")))
        }
    }
    ### Concatenate performances
    METRICS_ACROSS_POP_BULK = NULL
    YPRED_ACROSS_POP_BULK = NULL
    for (list_perf in list_list_perf) {
        # list_perf = list_list_perf[[1]]
        if (is.null(METRICS_ACROSS_POP_BULK) & is.null(YPRED_ACROSS_POP_BULK)) {
            METRICS_ACROSS_POP_BULK = list_perf$df_metrics
            YPRED_ACROSS_POP_BULK = list_perf$df_y_validation
        } else {
            METRICS_ACROSS_POP_BULK = rbind(METRICS_ACROSS_POP_BULK, list_perf$df_metrics)
            YPRED_ACROSS_POP_BULK = rbind(YPRED_ACROSS_POP_BULK, list_perf$df_y_validation)
        }
    }
    ### Save the concatenated output data frames across populations and clean-up
    time_rand_id = paste0(round(as.numeric(Sys.time())), sample.int(1e6, size=1))
    fname_across_bulk_Rds = file.path(dir_output, paste0("across_bulk-tempout-", time_rand_id, ".Rds"))
    saveRDS(list(
        METRICS_ACROSS_POP_BULK=METRICS_ACROSS_POP_BULK,
        YPRED_ACROSS_POP_BULK=YPRED_ACROSS_POP_BULK),
        file=fname_across_bulk_Rds)
    ### Clean-up temporary files generated by fn_cv_1 in parallel
    for (list_perf in list_list_perf) {
        unlink(list_perf$fname_metrics_out)
        unlink(list_perf$fname_y_validation_out)
    }
    ### Clean-up the memory used by subsetting the data set
    rm("list_list_perf")
    rm("list_perf")
    rm("METRICS_ACROSS_POP_BULK")
    rm("YPRED_ACROSS_POP_BULK")
    gc()
    ### Output
    return(fname_across_bulk_Rds)
}

#' Pairwise-population cross-validation
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
#' @param vec_models_to_test genomic prediction models to use (uses all the models below by default). Please choose one or more from:
#'  - "ridge": \eqn{Cost_{ridge} = \Sigma(y - X\beta)^2 + \lambda\Sigma\beta^2$, where $\hat{\beta} = {(X^TX + \lambda I)^{-1} X^Ty}}
#'      https://en.wikipedia.org/wiki/Ridge_regression
#'  - "lasso": \eqn{Cost_{lasso} = \Sigma(y - X\beta)^2 + \lambda\Sigma|\beta|}
#'      https://en.wikipedia.org/wiki/Lasso_(statistics)
#'  - "elastic_net": \eqn{Cost_{ridge} = \Sigma(y - X\beta)^2 + (1-\alpha)\lambda\Sigma\beta^2 + \alpha\lambda\Sigma|\beta|}
#'      https://en.wikipedia.org/wiki/Elastic_net_regularization
#'  - "Bayes_A": scaled t-distributed effects
#'      https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf
#'  - "Bayes_B": scaled t-distributed effects with probability \eqn{\pi}; and zero effects with probability \eqn{1-\pi}, 
#'      where \eqn{\pi \sim \beta(\theta_1, \theta_2)}.
#'      https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf
#'  - "Bayes_C": normally distributed effects (\eqn{N(0, \sigma^2_{\beta})}) with probability \eqn{\pi}; 
#'      and zero effects with probability \eqn{1-\pi}, where \eqn{\pi \sim \beta(\theta_1, \theta_2)}.
#'      https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf
#'  - "gBLUP": genotype best linear unbiased prediction (gBLUP) using genomic relationship matrix to predict missing breeding values
#'      via Direct-Inversion Newton-Raphson or Average Information 
#'      using the sommer R package (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4894563/).
#'      https://link.springer.com/protocol/10.1007/978-1-62703-447-0_13
#' @param bool_parallel perform multi-threaded cross-validation (Default=TRUE)
#' @param max_mem_Gb maximum memory in gigabytes available for computation (Default=15)
#' @param n_threads total number of computing threads available (Default=2)
#' @param dir_output output directory where temporary text and Rds files will be saved into (Default=NULL)
#' @param verbose show bulk across population cross-validation messages? (Default=FALSE)
#' @returns
#'  - Ok: filename of temporary Rds file containing a 2-element list:
#'      + $METRICS_ACROSS_POP_PAIRWISE:
#'          - $rep: replication number
#'          - $fold: fold number
#'          - $model: genomic prediction model name
#'          - $pop_training: population/s used in the training set (separated by commas if more than 1)
#'          - $pop_validation: population/s used in the validation set (separated by commas if more than 1)
#'          - $n_training: number of samples/entries/pools in the training set
#'          - $n_validation: number of samples/entries/pools in the validation set
#'          - $duration_mins: time taken in minutes to fit the genomic prediction model and assess the prediction accuracies
#'          - $n_non_zero: number of non-zero estimated effects (effects greater than machine epsilon ~2.2e-16)
#'          - $n_total_features: total number of zero and non-zero effects, i.e. total number of SNPs/alleles/features
#'          - $mbe: mean bias error
#'          - $mae: mean absolute error
#'          - $rmse: root mean squared error
#'          - $r2: coefficient of determination
#'          - $corr: Pearson's product moment correlation
#'          - $power_t10: fraction of observed top 10 phenotype values correctly predicted
#'          - $power_b10: fraction of observed bottom 10 phenotype values correctly predicted
#'          - $var_pred: variance of predicted phenotype values (estimator of additive genetic variance)
#'          - $var_true: variance of observed phenotype values (estimator of total phenotypic variance)
#'          - $h2: narrow-sense heritability estimate
#'      + $YPRED_ACROSS_POP_PAIRWISE:
#'          - $rep: replication number
#'          - $fold: fold number
#'          - $model: genomic prediction model name
#'          - $pop_training: population/s used in the training set (separated by commas if more than 1)
#'          - $id: names of the samples/entries/pools, 
#'          - $pop_validation: population from which the sample/entry/pool belongs to
#'          - $y_true: observed phenotype values
#'          - $y_pred: predicted phenotype values
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' COVAR = G %*% t(G)
#' list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
#' fname_across_pairwise_Rds = fn_cross_validation_across_populations_pairwise(list_merged, 
#'  vec_models_to_test=c("ridge","lasso"), verbose=TRUE)
#' list_across_pairwise = readRDS(fname_across_pairwise_Rds)
#' head(list_across_pairwise$METRICS_ACROSS_POP_PAIRWISE)
#' head(list_across_pairwise$YPRED_ACROSS_POP_PAIRWISE)
#' @export
fn_cross_validation_across_populations_pairwise = function(list_merged, 
    vec_models_to_test=c("ridge","lasso","elastic_net","Bayes_A","Bayes_B","Bayes_C","gBLUP"),
    bool_parallel=TRUE, max_mem_Gb=15, n_threads=2, dir_output=NULL, verbose=FALSE) 
{
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # COVAR = G %*% t(G)
    # list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    # vec_models_to_test = c("ridge", "lasso")
    # max_mem_Gb = 15
    # n_threads = 2
    # dir_output = NULL
    # bool_parallel = TRUE
    # verbose = TRUE
    ###################################################
    if (verbose) {
        print("####################################################")
        print("### Across populations cross-validation pairwise ###")
        print("####################################################")
    }
    ### Input sanity check
    if (methods::is(list_merged, "gpError")) {
        error = chain(list_merged, 
            methods::new("gpError",
                code=000,
                message=paste0(
                    "Error in cross_validation::fn_cross_validation_across_populations_pairwise(...). ",
                    "Input data (list_merged) is an error type."
                )))
        return(error)
    }
    ### Define the output directory
    if (!is.null(dir_output)) {
        if (!dir.exists(dir_output)) {
            dir.create(dir_output, showWarnings=FALSE)
        }
        
    } else {
        dir_output = tempfile()
        dir.create(dir_output, showWarnings=FALSE)
    }
    if (!dir.exists(dir_output)) {
        error = methods::new("gpError",
            code=000,
            message=paste0(
                "Error in cross_validation::fn_cross_validation_across_populations_pairwise(...). ",
                "Unable to create the output directory: ", dir_output, ". ",
                "Please check your permissions to write into that directory."
            ))
        return(error)
    }
    ### Determine the number of populations
    vec_populations = sort(unique(list_merged$list_pheno$pop))
    if (length(vec_populations) == 1) {
        error = methods::new("gpError",
            code=000,
            message=paste0(
                "Error in cross_validation::fn_cross_validation_across_populations_pairwise(...). ",
                "Cannot perform pairwise-population cross-validation ",
                "because there is only 1 population in the data set."
            ))
        return(error)
    }
    ### Instantiate the vector of Rds filenames containing the temporary output data per population
    vec_fname_across_pairwise_Rds = c()
    for (idx_pop1 in 1:(length(vec_populations)-1)) {
        for (idx_pop2 in (idx_pop1+1):length(vec_populations)) {
            # idx_pop1 = 1; idx_pop2 = 2
            population1 = vec_populations[idx_pop1]
            population2 = vec_populations[idx_pop2]
            ### Skip if both populations are the same which should not happen given the iterators above
            if (population1 == population2) {
                next
            }
            ### Subset the data per population
            vec_idx = which((list_merged$list_pheno$pop == population1) | (list_merged$list_pheno$pop == population2))
            list_merged_sub = fn_subset_merged_genotype_and_phenotype(list_merged=list_merged, vec_idx=vec_idx, verbose=verbose)
            if (methods::is(list_merged_sub, "gpError")) {
                error = chain(list_merged_sub, methods::new("gpError",
                    code=000,
                    message=paste0(
                        "Error in cross_validation::fn_cross_validation_across_populations_pairwise(...). ",
                        "Failed to subset the data set."
                    )))
                return(error)
            }
            ### Define the cross-validation parameters as well as the maximum number of threads we can safely use in parallel
            list_cv_params = fn_cross_validation_preparation(
                list_merged=list_merged_sub,
                cv_type=2,
                n_folds=1,
                n_reps=1, 
                vec_models_to_test=vec_models_to_test,
                max_mem_Gb=max_mem_Gb,
                verbose=verbose)
            if (methods::is(list_cv_params, "gpError")) {
                error = chain(list_cv_params, methods::new("gpError",
                    code=000,
                    message=paste0(
                        "Error in cross_validation::fn_cross_validation_within_population(...). ",
                        "Failed to define the cross-validation parameters."
                    )))
                return(error)
            }
            if (list_cv_params$list_mem$n_threads <= 1) {
                if (verbose) {
                    print("Data set is too large for parallel computations.")
                    print("Performing iterative cross-validation.")
                }
                bool_parallel = FALSE
            }
            ### Cross-validation
            if (bool_parallel) {
                ###########################################
                ### Multi-threaded cross-validation
                list_list_perf = parallel::mclapply(c(1:nrow(list_cv_params$df_params)), 
                    FUN=fn_cv_1, 
                        list_merged=list_merged_sub,
                        df_params=list_cv_params$df_params,
                        mat_idx_shuffle=list_cv_params$mat_idx_shuffle,
                        vec_set_partition_groupings=list_cv_params$vec_set_partition_groupings,
                        prefix_tmp=file.path(dir_output, paste0("across_pairwise-", population1, "-x-", population2)),
                        verbose=verbose, 
                    mc.cores=min(c(n_threads, list_cv_params$n_threads)))
                for (idx in 1:length(list_list_perf)) {
                    if (methods::is(list_list_perf[[idx]], "gpError")) {
                        error = chain(list_list_perf[[idx]], methods::new("gpError",
                            code=000,
                            message=paste0(
                                "Error in cross_validation::fn_cross_validation_across_populations_pairwise(...). ",
                                "Something went wrong in the execution of multi-threaded pairwise-population cross-validation. ",
                                "Please check re-run cross_validation::fn_cross_validation_across_populations_pairwise(...) with ",
                                "bool_parallel=FALSE to identify the error."
                            )))
                        return(error)
                    }
                }
            } else {
                ############################################
                ### Single-threaded cross-validation
                list_list_perf = list()
                for (i in 1:nrow(list_cv_params$df_params)) {
                    list_perf = fn_cv_1(
                        i=i,
                        list_merged=list_merged_sub,
                        df_params=list_cv_params$df_params,
                        mat_idx_shuffle=list_cv_params$mat_idx_shuffle,
                        vec_set_partition_groupings=list_cv_params$vec_set_partition_groupings,
                        prefix_tmp=file.path(dir_output, paste0("across_pairwise-", population1, "-x-", population2)),
                        verbose=verbose
                    )
                    if (methods::is(list_perf, "gpError")) {
                        error = chain(list_perf, methods::new("gpError",
                            code=000,
                            message=paste0(
                                "Error in cross_validation::fn_cross_validation_across_populations_pairwise(...). ",
                                "Error running pairwise cross-validation for populations: ", 
                                population1, " and ", population2, " at ",
                                "rep: ", list_cv_params$df_params$rep[i], ", ",
                                "fold: ", list_cv_params$df_params$fold[i], ", and ",
                                "model: ", list_cv_params$df_params$model[i], "."
                            )))
                        return(error)
                    }
                    eval(parse(text=paste0("list_list_perf$`", i, "` = list_perf")))
                }
            }
            ### Concatenate performances
            df_metrics = NULL
            df_y_validation = NULL
            for (list_perf in list_list_perf) {
                # list_perf = list_list_perf[[1]]
                if (is.null(df_metrics) & is.null(df_y_validation)) {
                    df_metrics = list_perf$df_metrics
                    df_y_validation = list_perf$df_y_validation
                } else {
                    df_metrics = rbind(df_metrics, list_perf$df_metrics)
                    df_y_validation = rbind(df_y_validation, list_perf$df_y_validation)
                }
            }
            ### Save temporary Rds output per population
            time_rand_id = paste0(round(as.numeric(Sys.time())), sample.int(1e6, size=1))
            fname_across_pairwise_Rds = file.path(dir_output, paste0("across_pairwise-tempout-", population1, "_x_", population2, "-", time_rand_id, ".Rds"))
            saveRDS(list(
                df_metrics=df_metrics,
                df_y_validation=df_y_validation), 
                file=fname_across_pairwise_Rds)
            ### Take note of the names of these temporary Rds output files
            vec_fname_across_pairwise_Rds = c(vec_fname_across_pairwise_Rds, fname_across_pairwise_Rds)
            ### Clean-up temporary files generated by fn_cv_1 in parallel
            for (list_perf in list_list_perf) {
                unlink(list_perf$fname_metrics_out)
                unlink(list_perf$fname_y_validation_out)
            }
            ### Clean-up the memory used by subsetting the data set
            list_merged_sub = NULL 
            list_list_perf = NULL 
            list_perf = NULL 
            gc()
        }
    }
    ### Concatenate the temporary output Rds files
    METRICS_ACROSS_POP_PAIRWISE = NULL
    YPRED_ACROSS_POP_PAIRWISE = NULL
    for (fname_across_pairwise_Rds in vec_fname_across_pairwise_Rds) {
        # fname_across_pairwise_Rds = vec_fname_across_pairwise_Rds[1]
        list_tmp = readRDS(fname_across_pairwise_Rds)
        if (is.null(METRICS_ACROSS_POP_PAIRWISE) & is.null(YPRED_ACROSS_POP_PAIRWISE)) {
            METRICS_ACROSS_POP_PAIRWISE = list_tmp$df_metrics
            YPRED_ACROSS_POP_PAIRWISE = list_tmp$df_y_validation
        } else {
            METRICS_ACROSS_POP_PAIRWISE = rbind(METRICS_ACROSS_POP_PAIRWISE, list_tmp$df_metrics)
            YPRED_ACROSS_POP_PAIRWISE = rbind(YPRED_ACROSS_POP_PAIRWISE, list_tmp$df_y_validation)
        }
    }
    ### Save the concatenated output data frames across populations and clean-up
    time_rand_id = paste0(round(as.numeric(Sys.time())), sample.int(1e6, size=1))
    fname_across_pairwise_Rds = file.path(dir_output, paste0("across_pairwise-tempout-", time_rand_id, ".Rds"))
    saveRDS(list(
        METRICS_ACROSS_POP_PAIRWISE=METRICS_ACROSS_POP_PAIRWISE,
        YPRED_ACROSS_POP_PAIRWISE=YPRED_ACROSS_POP_PAIRWISE),
        file=fname_across_pairwise_Rds)
    for (fname in vec_fname_across_pairwise_Rds) {
        unlink(fname)
    }
    rm("METRICS_ACROSS_POP_PAIRWISE")
    rm("YPRED_ACROSS_POP_PAIRWISE")
    gc()
    ### Output
    return(fname_across_pairwise_Rds)
}

#' Leave-one-population-out (LOPO) cross-validation
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
#' @param vec_models_to_test genomic prediction models to use (uses all the models below by default). Please choose one or more from:
#'  - "ridge": \eqn{Cost_{ridge} = \Sigma(y - X\beta)^2 + \lambda\Sigma\beta^2$, where $\hat{\beta} = {(X^TX + \lambda I)^{-1} X^Ty}}
#'      https://en.wikipedia.org/wiki/Ridge_regression
#'  - "lasso": \eqn{Cost_{lasso} = \Sigma(y - X\beta)^2 + \lambda\Sigma|\beta|}
#'      https://en.wikipedia.org/wiki/Lasso_(statistics)
#'  - "elastic_net": \eqn{Cost_{ridge} = \Sigma(y - X\beta)^2 + (1-\alpha)\lambda\Sigma\beta^2 + \alpha\lambda\Sigma|\beta|}
#'      https://en.wikipedia.org/wiki/Elastic_net_regularization
#'  - "Bayes_A": scaled t-distributed effects
#'      https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf
#'  - "Bayes_B": scaled t-distributed effects with probability \eqn{\pi}; and zero effects with probability \eqn{1-\pi}, 
#'      where \eqn{\pi \sim \beta(\theta_1, \theta_2)}.
#'      https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf
#'  - "Bayes_C": normally distributed effects (\eqn{N(0, \sigma^2_{\beta})}) with probability \eqn{\pi}; 
#'      and zero effects with probability \eqn{1-\pi}, where \eqn{\pi \sim \beta(\theta_1, \theta_2)}.
#'      https://cran.r-hub.io/web/packages/BGLR/vignettes/BGLR-extdoc.pdf
#'  - "gBLUP": genotype best linear unbiased prediction (gBLUP) using genomic relationship matrix to predict missing breeding values
#'      via Direct-Inversion Newton-Raphson or Average Information 
#'      using the sommer R package (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4894563/).
#'      https://link.springer.com/protocol/10.1007/978-1-62703-447-0_13
#' @param bool_parallel perform multi-threaded cross-validation (Default=TRUE)
#' @param max_mem_Gb maximum memory in gigabytes available for computation (Default=15)
#' @param n_threads total number of computing threads available (Default=2)
#' @param dir_output output directory where temporary text and Rds files will be saved into (Default=NULL)
#' @param verbose show bulk across population cross-validation messages? (Default=FALSE)
#' @returns
#'  - Ok: filename of temporary Rds file containing a 2-element list:
#'      + $METRICS_ACROSS_POP_LOPO:
#'          - $rep: replication number
#'          - $fold: fold number
#'          - $model: genomic prediction model name
#'          - $pop_training: population/s used in the training set (separated by commas if more than 1)
#'          - $pop_validation: population/s used in the validation set (separated by commas if more than 1)
#'          - $n_training: number of samples/entries/pools in the training set
#'          - $n_validation: number of samples/entries/pools in the validation set
#'          - $duration_mins: time taken in minutes to fit the genomic prediction model and assess the prediction accuracies
#'          - $n_non_zero: number of non-zero estimated effects (effects greater than machine epsilon ~2.2e-16)
#'          - $n_total_features: total number of zero and non-zero effects, i.e. total number of SNPs/alleles/features
#'          - $mbe: mean bias error
#'          - $mae: mean absolute error
#'          - $rmse: root mean squared error
#'          - $r2: coefficient of determination
#'          - $corr: Pearson's product moment correlation
#'          - $power_t10: fraction of observed top 10 phenotype values correctly predicted
#'          - $power_b10: fraction of observed bottom 10 phenotype values correctly predicted
#'          - $var_pred: variance of predicted phenotype values (estimator of additive genetic variance)
#'          - $var_true: variance of observed phenotype values (estimator of total phenotypic variance)
#'          - $h2: narrow-sense heritability estimate
#'      + $YPRED_ACROSS_POP_LOPO:
#'          - $rep: replication number
#'          - $fold: fold number
#'          - $model: genomic prediction model name
#'          - $pop_training: population/s used in the training set (separated by commas if more than 1)
#'          - $id: names of the samples/entries/pools, 
#'          - $pop_validation: population from which the sample/entry/pool belongs to
#'          - $y_true: observed phenotype values
#'          - $y_pred: predicted phenotype values
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' COVAR = G %*% t(G)
#' list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
#' fname_across_lopo_Rds = fn_cross_validation_across_populations_lopo(list_merged, 
#'  vec_models_to_test=c("ridge","lasso"), verbose=TRUE)
#' list_across_lopo = readRDS(fname_across_lopo_Rds)
#' head(list_across_lopo$METRICS_ACROSS_POP_LOPO)
#' head(list_across_lopo$YPRED_ACROSS_POP_LOPO)
#' @export
fn_cross_validation_across_populations_lopo = function(list_merged, 
    vec_models_to_test=c("ridge","lasso","elastic_net","Bayes_A","Bayes_B","Bayes_C","gBLUP"),
    bool_parallel=TRUE, max_mem_Gb=15, n_threads=2, dir_output=NULL, verbose=FALSE) 
{
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # COVAR = G %*% t(G)
    # list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    # vec_models_to_test = c("ridge", "lasso")
    # max_mem_Gb = 15
    # n_threads = 2
    # dir_output = NULL
    # bool_parallel = TRUE
    # verbose = TRUE
    ###################################################
    if (verbose) {
        print("####################################################################")
        print("### Across populations cross-validation leave-one-population-out ###")
        print("####################################################################")
    }
    ### Input sanity check
    if (methods::is(list_merged, "gpError")) {
        error = chain(list_merged, 
            methods::new("gpError",
                code=000,
                message=paste0(
                    "Error in cross_validation::fn_cross_validation_across_populations_lopo(...). ",
                    "Input data (list_merged) is an error type."
                )))
        return(error)
    }
    ### Define the output directory
    if (!is.null(dir_output)) {
        if (!dir.exists(dir_output)) {
            dir.create(dir_output, showWarnings=FALSE)
        }
        
    } else {
        dir_output = tempfile()
        dir.create(dir_output, showWarnings=FALSE)
    }
    if (!dir.exists(dir_output)) {
        error = methods::new("gpError",
            code=000,
            message=paste0(
                "Error in cross_validation::fn_cross_validation_across_populations_lopo(...). ",
                "Unable to create the output directory: ", dir_output, ". ",
                "Please check your permissions to write into that directory."
            ))
        return(error)
    }
    ### Define the cross-validation parameters as well as the maximum number of threads we can safely use in parallel
    list_cv_params = fn_cross_validation_preparation(
        list_merged=list_merged,
        cv_type=3,
        n_folds=1,
        n_reps=1, 
        vec_models_to_test=vec_models_to_test,
        max_mem_Gb=max_mem_Gb,
        verbose=verbose)
    if (methods::is(list_cv_params, "gpError")) {
        error = chain(list_cv_params, methods::new("gpError",
            code=000,
            message=paste0(
                "Error in cross_validation::fn_cross_validation_across_populations_lopo(...). ",
                "Failed to instantiate the cross-validation parameters."
            )))
        return(error)
    }
    if (list_cv_params$list_mem$n_threads <= 1) {
        if (verbose) {
            print("Data set is too large for parallel computations.")
            print("Performing iterative cross-validation.")
        }
        bool_parallel = FALSE
    }
    ### Cross-validation
    if (bool_parallel) {
        ###########################################
        ### Multi-threaded cross-validation
        list_list_perf = parallel::mclapply(c(1:nrow(list_cv_params$df_params)), 
            FUN=fn_cv_1, 
                list_merged=list_merged,
                df_params=list_cv_params$df_params,
                mat_idx_shuffle=list_cv_params$mat_idx_shuffle,
                vec_set_partition_groupings=list_cv_params$vec_set_partition_groupings,
                prefix_tmp=file.path(dir_output, "across_lopo"),
                verbose=verbose, 
            mc.cores=min(c(n_threads, list_cv_params$n_threads)))
        for (idx in 1:length(list_list_perf)) {
            if (methods::is(list_list_perf[[idx]], "gpError")) {
                error = chain(list_list_perf[[idx]], methods::new("gpError",
                    code=000,
                    message=paste0(
                        "Error in cross_validation::fn_cross_validation_across_populations_lopo(...). ",
                        "Something went wrong in the execution of multi-threaded pairwise-population cross-validation. ",
                        "Please check re-run cross_validation::fn_cross_validation_across_populations_lopo(...) with ",
                        "bool_parallel=FALSE to identify the error."
                    )))
                return(error)
            }
        }
    } else {
        ############################################
        ### Single-threaded cross-validation
        list_list_perf = list()
        for (i in 1:nrow(list_cv_params$df_params)) {
            list_perf = fn_cv_1(
                i=i,
                list_merged=list_merged,
                df_params=list_cv_params$df_params,
                mat_idx_shuffle=list_cv_params$mat_idx_shuffle,
                vec_set_partition_groupings=list_cv_params$vec_set_partition_groupings,
                prefix_tmp=file.path(dir_output, "across_lopo"),
                verbose=verbose
            )
            if (methods::is(list_perf, "gpError")) {
                error = chain(list_perf, methods::new("gpError",
                    code=000,
                    message=paste0(
                        "Error in cross_validation::fn_cross_validation_across_populations_lopo(...). ",
                        "Error running leave-one-population-out cross-validation at ",
                        "rep: ", list_cv_params$df_params$rep[i], ", ",
                        "fold: ", list_cv_params$df_params$fold[i], ", and ",
                        "model: ", list_cv_params$df_params$model[i], "."
                    )))
                return(error)
            }
            eval(parse(text=paste0("list_list_perf$`", i, "` = list_perf")))
        }
    }
    ### Concatenate performances
    df_metrics = NULL
    df_y_validation = NULL
    for (list_perf in list_list_perf) {
        # list_perf = list_list_perf[[1]]
        if (is.null(df_metrics) & is.null(df_y_validation)) {
            df_metrics = list_perf$df_metrics
            df_y_validation = list_perf$df_y_validation
        } else {
            df_metrics = rbind(df_metrics, list_perf$df_metrics)
            df_y_validation = rbind(df_y_validation, list_perf$df_y_validation)
        }
    }
    ### Save temporary Rds output per population
    time_rand_id = paste0(round(as.numeric(Sys.time())), sample.int(1e6, size=1))
    fname_across_lopo_Rds = file.path(dir_output, paste0("across_lopo-tempout-", time_rand_id, ".Rds"))
    saveRDS(list(
        METRICS_ACROSS_POP_LOPO=df_metrics,
        YPRED_ACROSS_POP_LOPO=df_y_validation), 
        file=fname_across_lopo_Rds)
    ### Clean-up temporary files generated by fn_cv_1 in parallel
    for (list_perf in list_list_perf) {
        unlink(list_perf$fname_metrics_out)
        unlink(list_perf$fname_y_validation_out)
    }
    ### Clean-up the memory used by subsetting the data set
    rm("list_list_perf")
    rm("list_perf")
    gc()
    ### Output
    return(fname_across_lopo_Rds)
}

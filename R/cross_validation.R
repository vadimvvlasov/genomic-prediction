# source("R/io.R")
# source("R/metrics.R")
# source("R/models.R")

#' Cross-validate on a single fold, replicate, and model
#'
#' @param list_merged list of merged genotype matrix, and phenotype vector, as well as an optional covariate matrix
#'  $G: numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'      Row names can be any string of characters which identify the sample or entry or pool names.
#'      Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'      the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'      subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'  $list_pheno:
#'      $y: named vector of numeric phenotype data
#'      $pop: population or groupings corresponding to each element of y
#'      $trait_name: name of the trait
#'  $COVAR: numeric n samples x k covariates matrix with non-null row and column names.
#' @param i index referring to the row in df_params (below) from which 
#'  the replicate id, fold number and model will be sourced from
#' @param df_params data frame containing all possible or just a defined non-exhaustive set of 
#'  combinations of replication, fold and model to be used in cross-validation
#' @param mat_idx_shuffle numeric n sample x r replications matrix of sample/entry/pool index shuffling where
#'  each column refer to a random shuffling of samples/entry/pool from which the identities of the
#'  training and validation sets will be sourced from
#' @param vec_set_partition_groupings vector of numeric partitioning indexes where each index refer to the
#'  fold which will serve as the validation population
#' @param prefix_tmp string referring to the prefix of the temporary files, 
#'  i.e. prefix (which can include an existing directory) of Bayesian (BGLR) model temporary files
#' @param verbose show cross-validation on a single fold, replicate, and model messages? (Default=FALSE)
#' @returns
#'  Ok:
#'      $df_metrics:
#'          $rep: replication number
#'          $fold: fold number
#'          $model: genomic prediction model name
#'          $pop_training: population/s used in the training set (separated by commas if more than 1)
#'          $pop_validation: population/s used in the validation set (separated by commas if more than 1)
#'          $duration_mins: time taken in minutes to fit the genomic prediction model and assess the prediction accuracies
#'          $n_non_zero: number of non-zero estimated effects (effects greater than machine epsilon ~2.2e-16)
#'          $mbe: mean bias error
#'          $mae: mean absolute error
#'          $rmse: root mean squared error
#'          $r2: coefficient of determination
#'          $corr: Pearson's product moment correlation
#'          $power_t10: fraction of observed top 10 phenotype values correctly predicted
#'          $power_b10: fraction of observed bottom 10 phenotype values correctly predicted
#'          $var_pred: variance of predicted phenotype values (estimator of additive genetic variance)
#'          $var_true: variance of observed phenotype values (estimator of total phenotypic variance)
#'          $h2: narrow-sense heritability estimate
#'      $df_y_validation: 
#'          $rep: replication number
#'          $fold: fold number
#'          $model: genomic prediction model name
#'          $pop_training: population/s used in the training set (separated by commas if more than 1)
#'          $id: names of the samples/entries/pools, 
#'          $pop_validation: population from which the sample/entry/pool belongs to
#'          $y_true: observed phenotype values
#'          $y_pred: predicted phenotype values
#'      $fname_metrics_out: filename of df_metrics saved as a the tab-delimited file with 2 rows
#'      $fname_y_validation_out: filename of df_y_validation saved as a the tab-delimited file
#'  Err: gpError
#' @examples
#' list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' COVAR = G %*% t(G)
#' n = nrow(G)
#' n_reps = 3
#' k_folds = 5
#' set_size = floor(n / k_folds)
#' vec_models_to_test = c("ridge", "lasso", "elastic_net", "Bayes_A", "Bayes_B", "Bayes_C", "gBLUP")
#' list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
#' df_params = expand.grid(rep=c(1:n_reps), fold=c(1:k_folds), model=vec_models_to_test)
#' mat_idx_shuffle = matrix(sample(1:n, size=n, replace=FALSE), ncol=1)
#' if (n_reps > 1) {
#'     for (r in 2:n_reps) {
#'         mat_idx_shuffle = cbind(mat_idx_shuffle, sample(1:n, size=n, replace=FALSE))
#'     }
#' }
#' vec_set_partition_groupings = rep(1:k_folds, each=set_size)
#' if (length(vec_set_partition_groupings) < n) {
#'     vec_set_partition_groupings = c(vec_set_partition_groupings, rep(k_folds, times=(n-length(vec_set_partition_groupings))))
#' }
#' list_cv_1 = fn_cv_1(
#'     list_merged=list_merged, 
#'     i=2, 
#'     df_params=df_params, 
#'     mat_idx_shuffle=mat_idx_shuffle, 
#'     vec_set_partition_groupings=vec_set_partition_groupings,
#'     prefix_tmp="gsTmp",
#'     verbose=TRUE)
#' @export
fn_cv_1 = function(list_merged, 
    i, df_params, mat_idx_shuffle, vec_set_partition_groupings, 
    prefix_tmp="gsTmp", verbose=FALSE){
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # COVAR = G %*% t(G)
    # n = nrow(G)
    # n_reps = 3
    # k_folds = 5
    # set_size = floor(n / k_folds)
    # vec_models_to_test = c("ridge", "lasso", "elastic_net", "Bayes_A", "Bayes_B", "Bayes_C", "gBLUP")
    # list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    # i = 1
    # df_params = expand.grid(rep=c(1:n_reps), fold=c(1:k_folds), model=vec_models_to_test)
    # mat_idx_shuffle = matrix(sample(1:n, size=n, replace=FALSE), ncol=1)
    # if (n_reps > 1) {
    #     for (r in 2:n_reps) {
    #         mat_idx_shuffle = cbind(mat_idx_shuffle, sample(1:n, size=n, replace=FALSE))
    #     }
    # }
    # vec_set_partition_groupings = rep(1:k_folds, each=set_size)
    # if (length(vec_set_partition_groupings) < n) {
    #     vec_set_partition_groupings = c(vec_set_partition_groupings, rep(k_folds, times=(n-length(vec_set_partition_groupings))))
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
    ### Define additional model input/s
    if (grepl("Bayes", model)==TRUE) {
        ### Append the input prefix into the temporary file prefix generated by Bayesian models so we don't overwite these when performing parallel computations
        other_params = list(nIter=12e3, burnIn=2e3, h2=0.5, out_prefix=paste0(dirname(prefix_tmp), "/bglr_", model, "-", basename(prefix_tmp), "-"))
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
        duration_mins=as.numeric(duration_mins),
        n_non_zero_effects=perf$n_non_zero,
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




















### K-fold cross validation across GP models with parallelisation across folds, replications and models
### Generates csv files with one row each for each fold x rep x model iteration. Clean these up if you like. These were just made to make sure we have output in case the entire function fails and we still have output to look at.
### Outputs a list containing a data frame of prediction accuracy metrics and a data frame of predicted phenotypes
fn_cross_validation = function(G, y, COVAR=NULL, vec_models_to_test=c("ridge","lasso","elastic_net","Bayes_A","Bayes_B","Bayes_C","gBLUP"), k_folds=10, n_reps=3, n_threads=32, mem_mb=64000, prefix_tmp="", verbose=FALSE) {
    ### Input data dimensions
    n = nrow(G)
    l = ncol(G)
    ### Standard normalise the phenotypes
    y_names = rownames(y)
    if (is.null(y_names)) {
        y_names = names(y)
    }
    y = scale(y, center=TRUE, scale=TRUE)
    rownames(y) = y_names
    ### Set partitioning
    set_size = floor(n / k_folds)
    if (set_size < 2) {
        print("Error: set size too small (less than 2 per fold). Please consider reducing the number of folds.")
        return(list(df_metrics=NA, df_y_pred=NA, fnames_metrics=NA, fnames_y_pred=NA))
    }
    vec_set_partition_groupings = rep(1:k_folds, each=set_size)
    if (length(vec_set_partition_groupings) < n) {
        vec_set_partition_groupings = c(vec_set_partition_groupings, rep(k_folds, times=(n-length(vec_set_partition_groupings))))
    }
    ### Prepare shuffling across replications
    mat_idx_shuffle = matrix(sample(1:n, size=n, replace=FALSE), ncol=1)
    if (n_reps > 1) {
        for (r in 2:n_reps) {
            mat_idx_shuffle = cbind(mat_idx_shuffle, sample(1:n, size=n, replace=FALSE))
        }
    }
    ### Prepare matrix of rep x fold x model combinations and sort by fold so that we can start comparing as soon as results become available
    df_params = expand.grid(rep=c(1:n_reps), fold=c(1:k_folds), model=vec_models_to_test)
    df_params = df_params[order(df_params$rep), ]
    df_params = df_params[order(df_params$fold), ]
    ### Limit the number of forks by the memory available and the size of G and y
    total_memory = mem_mb / 1e3 ### in gigabytes
    data_size = 2 * max(1, c(round(as.numeric(gsub(" Gb", "", format(utils::object.size(G), units="Gb"))) + as.numeric(gsub(" Gb", "", format(utils::object.size(y), units="Gb"))))))
    n_forks = max(c(1, min(c(nrow(df_params), n_threads, floor((total_memory-data_size) / data_size))))) ### force the minimum number of threads to 1
    ### Report the folds, replication, models being tested, memory allocation and number of parallel threads
    if (verbose==TRUE) {
        print("Genomic prediction cross-validation:")
        print(paste0("     - ", k_folds, " folds x ", n_reps, " replications with ", set_size, " entries per set."))
        print(paste0("     - ", n, " entries x ", l, " loci x (n_alleles - 1)."))
        print("Models being tested:")
        for (mod in vec_models_to_test) {
            print(paste0("     - ", mod))
        }
        print(paste0("Total memory: ", total_memory, " Gb"))
        print(paste0("Memory allocated per thread: ", data_size, " Gb"))
        print(paste0("Maximum number of threads to be used: ", n_forks))
    }
    ### Multi-threaded cross-fold validation
    time_ini = Sys.time()
    list_perf = tryCatch(
        parallel::mclapply(c(1:nrow(df_params)), 
            FUN=fn_cv_1, vec_set_partition_groupings=vec_set_partition_groupings, mat_idx_shuffle=mat_idx_shuffle, df_params=df_params, G=G, y=y, COVAR=COVAR, prefix_tmp=prefix_tmp, mem_mb=(floor(total_memory/n_forks)*1e3), 
            mc.cores=n_forks, mc.preschedule=FALSE, mc.set.seed=TRUE, mc.silent=FALSE, mc.cleanup=TRUE),
        error = function(e){c()}
    )
    # ### Recover from out-of-memory (OOM) errors when too many threads are being used.
    # ### Note that this will not work if we get a segmentation error, in which case simply reduce the number of CPUs (threads or cores, see`--cpus-per-task` in [`01_gs_slurm_job.sh`](./01_gs_slurm_job.sh).
    # while ((length(list_perf)==0) | (n_forks>1)) {
    #     n_forks = max(c(1, round(n_forks / 2)))
    #     print(paste0("Reducing the maximum number of threads to be used: ", n_forks))
    #     print(gc())
    #     list_perf = tryCatch(
    #         parallel::mclapply(c(1:nrow(df_params)), 
    #             FUN=fn_cv_1, vec_set_partition_groupings=vec_set_partition_groupings, mat_idx_shuffle=mat_idx_shuffle, df_params=df_params, G=G, y=y, COVAR=COVAR, prefix_tmp=prefix_tmp, mem_mb=(floor(total_memory/n_forks)*1e3), 
    #             mc.cores=n_forks, mc.preschedule=FALSE, mc.set.seed=TRUE, mc.silent=FALSE, mc.cleanup=TRUE),
    #         error = function(e){c()}
    #     )
    # }
    # list_perf = list()
    # for (i in c(1:nrow(df_params))) {
    #     print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    #     print(i)
    #     print(df_params[i, ])
    #     print(dim(COVAR))
    #     list_perf = c(list_perf, fn_cv_1(
    #         i=i,
    #         vec_set_partition_groupings=vec_set_partition_groupings,
    #         mat_idx_shuffle=mat_idx_shuffle,
    #         df_params=df_params,
    #         G=G,
    #         y=y,
    #         COVAR=COVAR,
    #         prefix_tmp=prefix_tmp,
    #         mem_mb=(floor(total_memory/n_forks)*1e3)
    #     ))
    # }
    duration_mins = difftime(Sys.time(), time_ini, units="min")
    print(paste0("Multi-threaded ", k_folds, "-fold cross-fold validation finished after ", duration_mins, " minutes."))
    ### Collect the output
    for (i in 1:length(list_perf)) {
        perf = list_perf[[i]]
        colnames(perf$df_y_pred) = c("rep", "entry", "model", "y_pred")
        if (i==1) {
            df_metrics = perf$df_metrics
            df_y_pred = perf$df_y_pred
            fnames_metrics = perf$fname_metrics_out
            fnames_y_pred = perf$fname_y_pred_out
        } else {
            df_metrics = rbind(df_metrics, perf$df_metrics)
            df_y_pred = rbind(df_y_pred, perf$df_y_pred)
            fnames_metrics = c(fnames_metrics, perf$fname_metrics_out)
            fnames_y_pred = c(fnames_y_pred, perf$fname_y_pred_out)
        }
    }
    df_y_true = data.frame(entry=rownames(y), y_true=y); colnames(df_y_true) = c("entry", "y_true"); rownames(df_y_true) = NULL
    df_y_pred = merge(df_y_pred, df_y_true, by="entry")
    ### Scatterplot of the observed and predicted phenotype data
    if (verbose==TRUE) {
        for (model in unique(df_y_pred$model)) {
            # model = unique(df_y_pred$model)[1]
            print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
            print(paste0(prefix_tmp, " - ", model))
            idx = which(df_y_pred$model == model)
            txtplot::txtplot(x=df_y_pred$y_true[idx], y=df_y_pred$y_pred[idx], xlab="Observed", ylab="Predicted")
        }
    }
    ### Output
    return(list(df_metrics=df_metrics, df_y_pred=df_y_pred, fnames_metrics=fnames_metrics, fnames_y_pred=fnames_y_pred))
}

### Pairwise population cross-validation, i.e. 1 population for training and the other for validation
### Generates csv files with one row each for each fold x rep x model iteration. Clean these up if you like. These were just made to make sure we have output in case the entire function fails and we still have output to look at.
### Outputs a list containing a data frame of prediction accuracy metrics and a data frame of predicted phenotypes
fn_pairwise_cross_validation = function(G_training, G_validation, y_training, y_validation, COVAR=NULL, vec_models_to_test=c("ridge","lasso","elastic_net","Bayes_A","Bayes_B","Bayes_C"), n_threads=32, mem_mb=64000, prefix_tmp="", verbose=FALSE) {
    ### Data dimensions
    for (suffix in c("training", "validation")) {
        eval(parse(text=paste0("n_", suffix, " = nrow(G_", suffix, ")")))
        eval(parse(text=paste0("l_", suffix, " = ncol(G_", suffix, ")")))
    }
    ### Standard normalise the phenotypes
    for (suffix in c("training", "validation")) {
        # suffix = "training"
        eval(parse(text=paste0("y_names_", suffix, " = rownames(y_", suffix, ")")))
        eval(parse(text=paste0("if (is.null(y_names_", suffix, ")){ y_names_", suffix, " = names(y_", suffix, ") }")))
        eval(parse(text=paste0("y_", suffix, " = scale(y_", suffix, ", center=TRUE, scale=TRUE)")))
        eval(parse(text=paste0("rownames(y_", suffix, ") = y_names_", suffix)))
    }
    ### Merge the training and validation sets, where the validation set is group 1
    G = rbind(G_validation, G_training)
    y = rbind(y_validation, y_training)
    ### Set partitioning where group 1 is the validation set, i.e. pop 1 and group 2 is the training set, i.e. pop 2
    vec_set_partition_groupings = c(rep(1, each=n_validation), rep(2, each=n_training))
    ### No shuffling needed as cross-validation is not replicated
    mat_idx_shuffle = matrix(1:(n_training+n_validation), ncol=1)
    ### Prepare matrix of rep x fold x model combinations
    df_params = expand.grid(rep=c(1), fold=c(1), model=vec_models_to_test)
    ### Report the folds, replication and models being tested
    if (verbose==TRUE) {
        print("Genomic prediction cross-validation:")
        print("Models being tested:")
        for (mod in vec_models_to_test) {
            print(paste0("     - ", mod))
        }
    }
    ### Limit the number of forks by the memory available and the size of G and y
    total_memory = mem_mb / 1e3 ### in gigabytes
    data_size = 50 * as.numeric(gsub(" Gb", "", format(utils::object.size(c(G, y)), units="Gb")))
    n_forks = min(c(n_threads, floor((total_memory-10) / data_size)))
    ### Multi-threaded cross-fold validation
    time_ini = Sys.time()
    list_perf = parallel::mclapply(c(1:nrow(df_params)), FUN=fn_cv_1, vec_set_partition_groupings=vec_set_partition_groupings, mat_idx_shuffle=mat_idx_shuffle, df_params=df_params, G=G, y=y, COVAR=COVAR, prefix_tmp=prefix_tmp, mem_mb=(floor(total_memory/n_forks)*1e3), mc.cores=n_forks)
    for (i in 1:nrow(df_params)) {
        i = 7
        x = fn_cv_1(i=i, vec_set_partition_groupings=vec_set_partition_groupings, mat_idx_shuffle=mat_idx_shuffle, df_params=df_params, G=G, y=y, COVAR=COVAR, prefix_tmp=prefix_tmp, mem_mb=(floor(total_memory/n_forks)*1e3))
    }
    duration_mins = difftime(Sys.time(), time_ini, units="min")
    print(paste0("Multi-threaded pairwise cross-fold validation finished after ", duration_mins, " minutes."))
    ### Collect the output
    for (i in 1:length(list_perf)) {
        perf = list_perf[[i]]
        colnames(perf$df_y_pred) = c("rep", "entry", "model", "y_pred")
        if (i==1) {
            df_metrics = perf$df_metrics
            df_y_pred = perf$df_y_pred
            fnames_metrics = perf$fname_metrics_out
            fnames_y_pred = perf$fname_y_pred_out
        } else {
            df_metrics = rbind(df_metrics, perf$df_metrics)
            df_y_pred = rbind(df_y_pred, perf$df_y_pred)
            fnames_metrics = c(fnames_metrics, perf$fname_metrics_out)
            fnames_y_pred = c(fnames_y_pred, perf$fname_y_pred_out)
        }
    }
    df_y_true = data.frame(entry=rownames(y), y_true=y); colnames(df_y_true) = c("entry", "y_true"); rownames(df_y_true) = NULL
    df_y_pred = merge(df_y_pred, df_y_true, by="entry")
    ### Scatterplot of the observed and predicted phenotype data
    if (verbose==TRUE) {
        for (model in unique(df_y_pred$model)) {
            # model = unique(df_y_pred$model)[1]
            print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
            print(paste0(prefix_tmp, " - ", model))
            idx = which(df_y_pred$model == model)
            txtplot::txtplot(x=df_y_pred$y_true[idx], y=df_y_pred$y_pred[idx], xlab="Observed", ylab="Predicted")
        }
    }
    ### Output
    return(list(df_metrics=df_metrics, df_y_pred=df_y_pred, fnames_metrics=fnames_metrics, fnames_y_pred=fnames_y_pred))
}

## Leave-one-population-out, i.e., using 1 population as the validation set and the rest as the training set
### Generates csv files with one row each for each fold x rep x model iteration. Clean these up if you like. These were just made to make sure we have output in case the entire function fails and we still have output to look at.
### Outputs a list containing a data frame of prediction accuracy metrics (which includes the validation population field corresponding to each fold-validation) and a data frame of predicted phenotypes
fn_leave_one_population_out_cross_validation = function(G, y, pop, COVAR=NULL, vec_models_to_test=c("ridge","lasso","elastic_net","Bayes_A","Bayes_B","Bayes_C"), n_threads=32, mem_mb=64000, prefix_tmp="", verbose=FALSE) {
    # n=100; l=1000; ploidy=100; n_alleles=2; min_allele_freq=0.01; n_chr=7; max_pos=1e6; dist_bp_at_50perc_r2=1e3; n_threads=8
    # G = simquantgen::fn_simulate_genotypes(n=n, l=l, ploidy=ploidy, n_alleles=n_alleles, min_allele_freq=min_allele_freq, n_chr=n_chr, max_pos=max_pos, dist_bp_at_50perc_r2=dist_bp_at_50perc_r2, n_threads=n_threads)
    # dist_effects=c("norm", "chi2")[1]; n_effects=100; purely_additive=FALSE; n_networks=10; n_effects_per_network=50; h2=0.5; pheno_reps=1
    # y = simquantgen::fn_simulate_phenotypes(G=G, n_alleles=n_alleles, dist_effects=dist_effects, n_effects=n_effects, purely_additive=purely_additive, n_networks=n_networks, n_effects_per_network=n_effects_per_network, h2=h2, pheno_reps=pheno_reps)$Y
    # # y = simquantgen::fn_simulate_phenotypes(G=G, n_alleles=n_alleles, dist_effects=dist_effects, n_effects=10, purely_additive=TRUE, n_networks=0, n_effects_per_network=0, h2=0.99, pheno_reps=1)$Y
    # # y = simquantgen::fn_simulate_phenotypes(G=G, n_alleles=n_alleles, dist_effects=dist_effects, n_effects=10, purely_additive=FALSE, n_networks=2, n_effects_per_network=10, h2=0.99, pheno_reps=1)$Y
    # pop = sample(c("popA", "popB", "popC", "popD"), size=n, replace=TRUE)
    # COVAR=NULL; vec_models_to_test=c("ridge","lasso","elastic_net","Bayes_A","Bayes_B","Bayes_C"); n_threads=32; mem_mb=64000; prefix_tmp=""
    # verbose = TRUE
    ### Data dimensions
    expect_equal(rownames(G), rownames(y))
    n = nrow(G)
    p = ncol(G)
    # ### Estimate kinship matrix which will be used as covariate
    # K = stats::cor(t(G))
    # ### Sort the samples by population
    # idx = order(pop)
    # G = G[idx, , drop=FALSE]
    # y = y[idx, , drop=FALSE]
    # expect_equal(rownames(G), rownames(y))
    # pop = pop[idx]
    ### Define the cross-validation partitioning where 1 population is set as validation set
    uniq_pop = sort(unique(pop))
    k = length(uniq_pop)
    vec_set_partition_groupings = rep(NA, n)
    for (i in 1:k) {
        idx = pop==uniq_pop[i]
        vec_set_partition_groupings[idx] = i
    }
    ### No shuffling needed as cross-validation is not replicated
    mat_idx_shuffle = matrix(1:n, ncol=1)
    ### Prepare matrix of rep x fold x model combinations
    df_params = expand.grid(rep=c(1), fold=c(1:k), model=vec_models_to_test)
    ### Report the folds, replication and models being tested
    if (verbose==TRUE) {
        print("Genomic prediction cross-validation:")
        print("Models being tested:")
        for (mod in vec_models_to_test) {
            print(paste0("     - ", mod))
        }
    }
    ### Limit the number of forks by the memory available and the size of G and y
    total_memory = mem_mb * 1e3
    data_size = 50 * as.numeric(gsub(" Gb", "", format(utils::object.size(c(G, y)), units="Gb")))
    n_forks = min(c(n_threads, floor((total_memory-10) / data_size), length(vec_models_to_test)))
    ### Multi-threaded cross-fold validation (Note: y-predictions are standard normalised)
    time_ini = Sys.time()
    list_perf = parallel::mclapply(c(1:nrow(df_params)), FUN=fn_cv_1, vec_set_partition_groupings=vec_set_partition_groupings, mat_idx_shuffle=mat_idx_shuffle, df_params=df_params, G=G, y=y, COVAR=COVAR, prefix_tmp=prefix_tmp, mem_mb=(floor(total_memory/n_forks)*1e3), mc.cores=n_forks)
    duration_mins = difftime(Sys.time(), time_ini, units="min") # 2.8 hours!
    print(paste0("Multi-threaded leave-one-population-out cross-fold validation finished after ", duration_mins, " minutes."))
    ### Collect the output
    for (i in 1:length(list_perf)) {
        perf = list_perf[[i]]
        colnames(perf$df_y_pred) = c("rep", "entry", "model", "y_pred")
        tmp_df_metrics = perf$df_metrics
        tmp_df_y_pred = perf$df_y_pred
        tmp_df_y_pred$pop = uniq_pop[df_params$fold[i]]
        tmp_fnames_metrics = perf$fname_metrics_out
        tmp_fnames_y_pred = perf$fname_y_pred_out
        if (i==1) {
            df_metrics = tmp_df_metrics
            df_y_pred = tmp_df_y_pred
            fnames_metrics = tmp_fnames_metrics
            fnames_y_pred = tmp_fnames_y_pred
        } else {
            df_metrics = rbind(df_metrics, tmp_df_metrics)
            df_y_pred = rbind(df_y_pred, tmp_df_y_pred)
            fnames_metrics = c(fnames_metrics, tmp_fnames_metrics)
            fnames_y_pred = c(fnames_y_pred, tmp_fnames_y_pred)
        }
    }
    df_y_true = data.frame(entry=rownames(y), y_true=y); colnames(df_y_true) = c("entry", "y_true"); rownames(df_y_true) = NULL
    df_y_pred = merge(df_y_pred, df_y_true, by="entry")
    ### Scatterplot of the observed and predicted phenotype data
    if (verbose==TRUE) {
        for (model in unique(df_y_pred$model)) {
            # model = unique(df_y_pred$model)[1]
            print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
            print(paste0(prefix_tmp, " - ", model))
            idx = which(df_y_pred$model == model)
            txtplot::txtplot(x=df_y_pred$y_true[idx], y=df_y_pred$y_pred[idx], xlab="Observed", ylab="Predicted")
        }
    }
    ### Output
    df_metrics$validation_pop = uniq_pop[df_metrics$k]
    return(list(df_metrics=df_metrics, df_y_pred=df_y_pred, fnames_metrics=fnames_metrics, fnames_y_pred=fnames_y_pred))
}

### Within population, across populations, and per se genomic predictions
fn_within_across_perse_genomic_prediction = function(G, idx_col_y, args, dir_tmp) {
    # G=X; idx_col_y=vec_idx_col_y[1]; args=args
    ### Load phenotype data
    print("##################################################")
    print("Load phenotype & exclude outliers")
    list_y_pop = fn_load_phenotype(fname_csv_txt=args$fname_pheno,
                                sep=args$sep,
                                header=args$header,
                                idx_col_id=args$idx_col_id,
                                idx_col_pop=args$idx_col_pop,
                                idx_col_y=idx_col_y,
                                na.strings=args$na_strings)
    print(list_y_pop$trait_name)
    ### Remove phenotype outlier/s
    list_y_pop = fn_filter_outlying_phenotypes(list_y_pop=list_y_pop,
                                            verbose=args$verbose)
    ### Load covariate
    if (is.null(args$fname_covar_rds) == TRUE) {
        COVAR = NULL
    } else {
        print("##################################################")
        print("Load covariate and standard normalise per column")
        COVAR = readRDS(args$fname_covar_rds)
        COVAR = scale(COVAR, scale=TRUE, center=TRUE)
    }
    ### Merge genotype, phenotype, and covariate matrices, where they all need to be matrices with overlapping rownames i.e., sample names (requirement of fn_cross_validation(...))
    print("##################################################")
    print("Merge data")
    list_G_list_y_pop_COVAR = fn_merge_genotype_and_phenotype(G=G,
                                                            list_y_pop=list_y_pop,
                                                            COVAR=COVAR,
                                                            verbose=args$verbose)
    idx_no_missing = which(!is.na(list_G_list_y_pop_COVAR$list_y_pop$y))
    y = matrix(list_G_list_y_pop_COVAR$list_y_pop$y[idx_no_missing], ncol=1)
    G = list_G_list_y_pop_COVAR$G[idx_no_missing, ]
    rownames(y) = names(list_G_list_y_pop_COVAR$list_y_pop$y[idx_no_missing])
    pop = list_G_list_y_pop_COVAR$list_y_pop$pop[idx_no_missing]
    COVAR = list_G_list_y_pop_COVAR$COVAR[idx_no_missing, ]
    vec_pop = sort(unique(pop))
    if (args$vec_pop != "all") {
        vec_pop_requested = unlist(strsplit(args$vec_pop, ","))
        vec_pop = intersect(vec_pop, vec_pop_requested)
        if (length(vec_pop) == 0) {
            print("Requested populations does not exist in the input data: ")
            print(args$vec_pop)
            print("The data has the following populations: ")
            print(paste(sort(unique(pop)), collapse=","))
            quit()
        }
        print("The following populations will be included in the analysis: ")
        print(paste(vec_pop, collapse=","))
    }
    idx_entries_for_genomic_prediction = which(is.na(list_G_list_y_pop_COVAR$list_y_pop$y) & (list_G_list_y_pop_COVAR$list_y_pop$pop %in% vec_pop))
    print(paste0("Total number of samples: ", nrow(G) + length(idx_entries_for_genomic_prediction)))
    print(paste0("Number of samples with phenotype data: ", nrow(G)))
    print(paste0("Number of samples whose phenotypes will be predicted: ", length(idx_entries_for_genomic_prediction)))
    print(paste0("Total number of markers: ", ncol(G)))
    print(paste0("Number of populations: ", length(vec_pop)))
    if (is.null(COVAR)) {
        print("Number of covariates: 0")
    } else {
        print(paste0("Number of covariates: ", ncol(COVAR)))
    }
    ##################################################
    ### WITHIN POPULATION: k-fold cross-validation ###
    ##################################################
    print("##################################################")
    print("Within population cross-validation:")
    vec_fnames_metrics = c()
    vec_fnames_y_pred = c()
    for (pop_id in vec_pop) {
        # pop_id = vec_pop[1]
        print("-----------------------------------")
        idx = which(pop == pop_id)
        print(paste0("Population: ", pop_id, " (n=", length(idx), " x p=", ncol(G), ")"))
        prefix = file.path(dir_tmp, gsub(".vcf.gz$", "", ignore.case=TRUE, gsub(".vcf$", "", ignore.case=TRUE, gsub(".rds$", "", ignore.case=TRUE, basename(args$fname_rds_or_vcf)))))
        list_out = fn_cross_validation(G=G[idx, , drop=FALSE],
                                    y=y[ idx, , drop=FALSE],
                                    COVAR=COVAR[idx, , drop=FALSE],
                                    vec_models_to_test=args$vec_models_to_test,
                                    k_folds=args$k_folds,
                                    n_reps=args$n_reps,
                                    n_threads=args$n_threads,
                                    mem_mb=args$mem_mb,
                                    prefix_tmp=paste0(prefix, "-pop_", pop_id),
                                    verbose=args$verbose)
        if (length(list_out$df_metrics) == 1) {
            ### Quit, if cross-validation failed because the number size of each set is less than 2.
            ### The error message from fn_cross_validation(...) will recommend decreasing k_folds.
            if (is.na(list_out$df_metrics)) {
                quit()
            }
        }
        list_out$df_metrics$pop = pop_id
        vec_fnames_metrics = c(vec_fnames_metrics, list_out$fnames_metrics)
        vec_fnames_y_pred = c(vec_fnames_y_pred, list_out$fnames_y_pred)
    }
    METRICS_WITHIN_POP = list_out$df_metrics
    YPRED_WITHIN_POP = list_out$df_y_pred
    ### Cleanup
    for (i in 1:length(vec_fnames_metrics)) {
        unlink(vec_fnames_metrics[i])
        unlink(vec_fnames_y_pred[i])
    }
    #####################################################################
    ### ACROSS POPULATIONS: leave-one-population-out cross-validation ###
    #####################################################################
    if ((args$skip_lopo_cv==TRUE) | (length(unique(pop))==1)) {
        args$skip_lopo_cv = TRUE
        METRICS_ACROSS_POP_LOPO = NULL
        YPRED_ACROSS_POP_LOPO = NULL
        print("Skipping leave-one-population-out cross-validation as only one population was supplied.")
    } else {
        print("##################################################")
        print("Across population cross-validation:")
        print("(Leave-one-population-out CV)")



        ### REDUCE NUMBER OF THREADS HERE ARE DATA SIZE differs from WITHIN POP CV!
        ### USE MAXIMUM NUMBER OF THREADS ESTIMATION



        list_out = fn_leave_one_population_out_cross_validation(G=G, 
                                                                y=y,
                                                                pop=pop,
                                                                COVAR=COVAR,
                                                                vec_models_to_test=args$vec_models_to_test,
                                                                n_threads=args$n_threads,
                                                                mem_mb=args$mem_mb,
                                                                prefix_tmp=paste0(prefix, "-LOPO"),
                                                                verbose=args$verbose)
        METRICS_ACROSS_POP_LOPO = list_out$df_metrics
        YPRED_ACROSS_POP_LOPO = list_out$df_y_pred
        ### Cleanup
        for (i in 1:length(list_out$fnames_metrics)) {
            unlink(list_out$fnames_metrics[i])
            unlink(list_out$fnames_y_pred[i])
        }
    }
    ################################################################
    ### ACROSS POPULATIONS: pairwise population cross-validation ###
    ################################################################
    vec_pop = sort(unique(pop))
    METRICS_ACROSS_POP_PAIRWISE = NULL
    YPRED_ACROSS_POP_PAIRWISE = NULL
    if ((args$skip_pairwise_cv==TRUE) | (length(unique(pop))==1)) {
        args$skip_pairwise_cv = TRUE
        print("Skipping leave-one-population-out cross-validation as only one population was supplied.")
    } else {
        print("##################################################")
        print("Across population cross-validation:")
        print("(Pairwise population CV)")
        vec_fnames_metrics = c()
        vec_fnames_y_pred = c()
        for (pop_id_1 in vec_pop) {
            for (pop_id_2 in vec_pop) {
                # pop_id_1 = vec_pop[1]; pop_id_2 = vec_pop[2]
                idx_1 = which(pop == pop_id_1)
                idx_2 = which(pop == pop_id_2)
                if (pop_id_1==pop_id_2) {
                    next
                }
                print("-----------------------------------")
                print(paste0("Training population: ", pop_id_1, " (n=", length(idx_1), " x p=", ncol(G), ")"))
                print(paste0("Validation population: ", pop_id_2, " (n=", length(idx_2), " x p=", ncol(G), ")"))
                prefix = file.path(dir_tmp, gsub(".vcf.gz$", "", ignore.case=TRUE, gsub(".vcf$", "", ignore.case=TRUE, gsub(".rds$", "", ignore.case=TRUE, basename(args$fname_rds_or_vcf)))))
                list_out = fn_pairwise_cross_validation(G_training=G[idx_1, , drop=FALSE], 
                                                        G_validation=G[idx_2, , drop=FALSE],
                                                        y_training=y[ idx_1, , drop=FALSE],
                                                        y_validation=y[ idx_2, , drop=FALSE],
                                                        COVAR=COVAR[idx, , drop=FALSE],
                                                        vec_models_to_test=args$vec_models_to_test,
                                                        n_threads=args$n_threads,
                                                        mem_mb=args$mem_mb,
                                                        prefix_tmp=paste0(prefix, "-PAIRWISE-pop1_", pop_id_1, "-pop2_", pop_id_2),
                                                        verbose=args$verbose)
                if (length(list_out$df_metrics) == 1) {
                    ### Quit, if cross-validation failed because the number size of each set is less than 2.
                    ### The error message from fn_cross_validation(...) will recommend decreasing k_folds.
                    if (is.na(list_out$df_metrics)) {
                        quit()
                    }
                }
                list_out$df_metrics$training_pop = pop_id_1
                list_out$df_metrics$validation_pop = pop_id_2
                list_out$df_y_pred$training_pop = pop_id_1
                list_out$df_y_pred$validation_pop = pop_id_2
                if (is.null(METRICS_ACROSS_POP_PAIRWISE) & is.null(YPRED_ACROSS_POP_PAIRWISE)) {
                    METRICS_ACROSS_POP_PAIRWISE = list_out$df_metrics
                    YPRED_ACROSS_POP_PAIRWISE = list_out$df_y_pred
                } else {
                    METRICS_ACROSS_POP_PAIRWISE = rbind(METRICS_ACROSS_POP_PAIRWISE, list_out$df_metrics)
                    YPRED_ACROSS_POP_PAIRWISE = rbind(YPRED_ACROSS_POP_PAIRWISE, list_out$df_y_pred)
                }
                vec_fnames_metrics = c(vec_fnames_metrics, list_out$fnames_metrics)
                vec_fnames_y_pred = c(vec_fnames_y_pred, list_out$fnames_y_pred)
            }
        }
        ### Cleanup
        for (i in 1:length(vec_fnames_metrics)) {
            unlink(vec_fnames_metrics[i])
            unlink(vec_fnames_y_pred[i])
        }
    }
    ######################################################################
    ### PER SE GENOMIC PREDICTION: using the best model per population ###
    ######################################################################
    ### Best models per population
    print("##################################################")
    print("Identifying the best models per population")
    print("(using Pearson's correlation)")
    ### Only considering within population k-fold cross-validation to identify the best model per population
    ### Using Pearson's correlation as the genomic prediction accuracy metric
    agg_corr = stats::aggregate(corr ~ model + pop, FUN=mean, data=METRICS_WITHIN_POP)
    agg_corr_sd = stats::aggregate(corr ~ model + pop, FUN=stats::sd, data=METRICS_WITHIN_POP)
    agg_corr = merge(agg_corr, agg_corr_sd, by=c("model", "pop"))
    colnames(agg_corr) = c("model", "pop", "corr", "corr_sd")
    vec_popns = c()
    vec_model = c()
    vec_corr = c()
    vec_corr_sd = c()
    for (pop in unique(agg_corr$pop)) {
        agg_sub = agg_corr[agg_corr$pop==pop, ]
        idx = which(agg_sub$corr == max(agg_sub$corr, na.rm=TRUE))[1]
        model = as.character(agg_sub$model[idx])
        corr = agg_sub$corr[idx]
        corr_sd = agg_sub$corr_sd[idx]
        vec_popns = c(vec_popns, pop)
        vec_model = c(vec_model, model)
        vec_corr = c(vec_corr, corr)
        vec_corr_sd = c(vec_corr_sd, corr_sd)
    }
    ### Append overall best model across all populations
    agg_overall_mean = stats::aggregate(corr ~ model, FUN=mean, data=METRICS_WITHIN_POP)
    agg_overall_sdev = stats::aggregate(corr ~ model, FUN=stats::sd, data=METRICS_WITHIN_POP)
    agg_overall = merge(agg_overall_mean, agg_overall_sdev, by="model")
    colnames(agg_overall) = c("model", "corr", "corr_sd")
    agg_overall = agg_overall[which(agg_overall[,2]==max(agg_overall[,2], na.rm=TRUE))[1], ]
    vec_popns = c(vec_popns, "overall")
    vec_model = c(vec_model, as.character(agg_overall$model))
    vec_corr = c(vec_corr, agg_overall$corr)
    vec_corr_sd = c(vec_corr_sd, agg_overall$corr_sd)
    SUMMARY = data.frame(pop=vec_popns, model=vec_model, corr=vec_corr, corr_sd=vec_corr_sd)
    print(SUMMARY)
    ### Genomic prediction per se
    if (length(idx_entries_for_genomic_prediction) == 0) {
        GENOMIC_PREDICTIONS = NULL
        print("No entries with missing phenotype and known genotypes detected in the input dataset.")
        print("Entries with missing phenotype do not belong to the included populations.")
    } else {
        print("##################################################")
        print("Genomic prediction per se using the best models")
        print("(predict missing phenotypes using known genotypes)")
        print(paste0("For each entry belonging to a population assessed via within population ", args$k_folds, "-fold cross-validation,"))
        print("the best model in terms of Person's correlation will be used to predict their unknown phenotypes.")
        vec_popns = c()
        vec_entry = c()
        vec_y_pred = c()
        vec_model = c()
        vec_corr = c()
        vec_corr_sd = c()
        vec_corr_pop = c()
        vec_pop_for_prediction_per_se = vec_pop[vec_pop %in% unique(list_G_list_y_pop_COVAR$list_y_pop$pop[idx_entries_for_genomic_prediction])]
        for (pop in vec_pop_for_prediction_per_se) {
            # pop = vec_pop_for_prediction_per_se[1]
            idx = which(list_G_list_y_pop_COVAR$list_y_pop$pop==pop)
            G = list_G_list_y_pop_COVAR$G[idx, ]
            y = matrix(list_G_list_y_pop_COVAR$list_y_pop$y[idx], ncol=1)
            rownames(y) = names(list_G_list_y_pop_COVAR$list_y_pop$y[idx])
            if (is.null(list_G_list_y_pop_COVAR$COVAR) == TRUE) {
                COVAR = NULL
            } else {
                COVAR = list_G_list_y_pop_COVAR$COVAR[idx, ]
            }
            idx_training = which(!is.na(y[,1]))
            idx_validation = which(is.na(y[,1]))
            idx_top_model = which(SUMMARY$pop == pop)[1]
            if (length(idx_top_model)==0) {
                ### Use the overall top model if the entries do not belong to any other populations
                idx_top_model = which(SUMMARY$pop == "overall")[1]
            }
            model = SUMMARY$model[idx_top_model]
            corr = SUMMARY$corr[idx_top_model]
            corr_sd = SUMMARY$corr_sd[idx_top_model]
            corr_pop = SUMMARY$pop[idx_top_model]

            if (grepl("Bayes", model)==TRUE) {
                other_params = list(nIter=12e3, burnIn=2e3, h2=0.5, out_prefix=paste0("bglr_", model, "-"), covariate=COVAR)
            } else {
                other_params = list(covariate=COVAR)
            }
            gp = eval(parse(text=paste0("fn_", model, "(G=G, y=y, idx_training=idx_training, idx_validation=idx_validation, other_params=other_params)")))
            print(gp)
            vec_popns = c(vec_popns, rep(pop, length(gp$y_pred)))
            if (is.null(names(gp$y_pred))) {
                vec_entry = c(vec_entry, rownames(gp$y_pred))
            } else {
                vec_entry = c(vec_entry, names(gp$y_pred))
            }
            vec_y_pred = c(vec_y_pred, gp$y_pred)
            vec_model = c(vec_model, rep(model, length(gp$y_pred)))
            vec_corr = c(vec_corr, rep(corr, length(gp$y_pred)))
            vec_corr_sd = c(vec_corr_sd, rep(corr_sd, length(gp$y_pred)))
            vec_corr_pop = c(vec_corr_pop, rep(corr_pop, length(gp$y_pred)))
        }
        GENOMIC_PREDICTIONS = data.frame(pop=vec_popns, entry=vec_entry, y_pred=vec_y_pred, top_model=vec_model, corr_from_kfoldcv=paste0(round(vec_corr, 2), "(±", round(vec_corr_sd, 2), "|", vec_corr_pop, ")"))
    }
    ### Output    
    out_per_phenotype = list(
        TRAIT_NAME                  = list_y_pop$trait_name,
        SUMMARY                     = SUMMARY,
        METRICS_WITHIN_POP          = METRICS_WITHIN_POP,
        YPRED_WITHIN_POP            = YPRED_WITHIN_POP,
        METRICS_ACROSS_POP_LOPO     = METRICS_ACROSS_POP_LOPO,
        YPRED_ACROSS_POP_LOPO       = YPRED_ACROSS_POP_LOPO,
        METRICS_ACROSS_POP_PAIRWISE = METRICS_ACROSS_POP_PAIRWISE,
        YPRED_ACROSS_POP_PAIRWISE   = YPRED_ACROSS_POP_PAIRWISE,
        GENOMIC_PREDICTIONS         = GENOMIC_PREDICTIONS
    )
    return(out_per_phenotype)
}

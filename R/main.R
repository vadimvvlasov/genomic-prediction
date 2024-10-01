#' Main genomic prediction cross-validation within and across populations
#'
#' @param args list of arguments or inputs
#'  - $fname_geno: filename of the input genotype file (see ?fn_load_genotype for details)
#'  - $fname_pheno: filename of the input phenotype file (see ?fn_load_phenotype for details)
#'  - $population: names of the population to used in within population k-fold cross-validation
#'  - $fname_covar: filename of the covariate data (TODO: implement fn_load_covariate)
#'  - $dir_output: output directory into which all temporary and final output files will be written
#'  - $geno_fname_snp_list: filename of the SNP list for filtering the genotype data
#'      (see ?fn_filter_genotype for details)
#'  - $geno_ploidy: expected ploidy level of the genotype file (see ?fn_load_genotype for details)
#'  - $geno_bool_force_biallelic: force all loci including multi-allelic ones to be biallelic?
#'      (see ?fn_load_genotype for details)
#'  - $geno_bool_retain_minus_one_alleles_per_locus: remove the trailing allele per locus?
#'      (see ?fn_load_genotype for details)
#'  - $geno_min_depth: minimum depth per allele (see ?fn_load_genotype for details)
#'  - $geno_max_depth: maximum depth per allele (see ?fn_load_genotype for details)
#'  - $geno_maf: minimum allele frequency (see ?fn_filter_genotype for details)
#'  - $geno_sdev_min: minimum allele frequency standard deviation (see ?fn_filter_genotype for details)
#'  - $geno_max_n_alleles: maximum number of alleles per locus (see ?fn_filter_genotype for details)
#'  - $geno_max_sparsity_per_locus: maximum sparsity or fraction of missing data per locus
#'      (see ?fn_filter_genotype for details)
#'  - $geno_frac_topmost_sparse_loci_to_remove: fraction of the total number of loci from which the 
#'      top most sparse loci will be removed (see ?fn_filter_genotype for details)
#'  - $geno_n_topmost_sparse_loci_to_remove: number of top most sparse loci to be removed
#'      (see ?fn_filter_genotype for details)
#'  - $geno_max_sparsity_per_sample: maximum sparsity or fraction of missing data per sample
#'      (see ?fn_filter_genotype for details)
#'  - $geno_frac_topmost_sparse_samples_to_remove: fraction of the total number of samples from which 
#'      the top most sparse samples will be removed (see ?fn_filter_genotype for details)
#'  - $geno_n_topmost_sparse_samples_to_remove: number of top most sparse samples to be removed
#'      (see ?fn_filter_genotype for details)
#'  - $pheno_sep: delimiter used in the phenotype file (see ?fn_load_phenotype for details)
#'  - $pheno_header: does the phenotype file have a header line? (see ?fn_load_phenotype for details)
#'  - $pheno_idx_col_id: column number in the phenotype file corresponding to the sample names
#'      (see ?fn_load_phenotype for details)
#'  - $pheno_idx_col_pop: column number in the phenotype file corresponding to the population/grouping names
#'      (see ?fn_load_phenotype for details)
#'  - $pheno_idx_col_y: column number in the phenotype file corresponding to the numeric phenotype data
#'      (see ?fn_load_phenotype for details)
#'  - $pheno_vec_na_strings: strings of characters corresponding to missing data in the phenotype file
#'      (see ?fn_load_phenotype for details)
#'  - $pheno_bool_remove_outliers: remove outliers from the phenotype file?
#'  - $pheno_bool_remove_NA: remove samples missing phenotype data in the phenotype file?
#'      (see ?fn_load_phenotype for details)
#'  - $bool_within: perform within population k-fold cross-validation?
#'      (see ?fn_cross_validation_within_population for details)
#'  - $bool_across: perform across populations cross-validations?
#'      (see ?fn_cross_validation_across_populations_bulk,
#'      ?fn_cross_validation_across_populations_bulk,
#'      ?fn_cross_validation_across_populations_pairwise, and 
#'      ?fn_cross_validation_across_populations_lopo for details)
#'  - $n_folds: number of folds or training and validation sets for within population k-fold cross validation
#'      (see ?fn_cross_validation_within_population for details)
#'  - $n_reps: number of replications of random shuffling samples for each k-fold cross-validation 
#'      in within population cross-validation (see ?fn_cross_validation_within_population for details)
#'  - $vec_models_to_test: genomic prediction models to use 
#'      (see ?fn_cross_validation_preparation for details, as well as the individual models:
#'      ?fn_ridge, ?fn_lasso, ?fn_elastic_net, ?fn_Bayes_A, ?fn_Bayes_B, ?fn_Bayes_C, and ?fn_gBLUP)
#'  - $bool_parallel: run the cross-validations in parallel (see the following functions for details:
#'      ?fn_cross_validation_within_population,
#'      ?fn_cross_validation_across_populations_bulk,
#'      ?fn_cross_validation_across_populations_pairwise, and
#'      ?fn_cross_validation_across_populations_lopo)
#'  - $max_mem_Gb: maximum amount of memory available in gigabytes
#'      (see ?fn_estimate_memory_footprint for details)
#'  - $n_threads: maximum number of computing threads available
#'      (see ?fn_estimate_memory_footprint for details)
#'  - $verbose: show messages?
#' @returns
#'  - Ok:
#'      + $TRAIT_NAME: name of the phenotype or trait used to fit single-trait genomic prediction models 
#'          (multi-trait models will be implemented in genomic_breeding)
#'      + $POPULATION: name of the population used in within population genomic prediction k-fold cross-validation
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
#'          - $var_additive: variance of predicted phenotype values (estimator of additive genetic variance)
#'          - $var_residual: variance of differnce between observed and predicted phenotype values (estimator of residual variance)
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
#'          - $var_additive: variance of predicted phenotype values (estimator of additive genetic variance)
#'          - $var_residual: variance of differnce between observed and predicted phenotype values (estimator of residual variance)
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
#'          - $var_additive: variance of predicted phenotype values (estimator of additive genetic variance)
#'          - $var_residual: variance of differnce between observed and predicted phenotype values (estimator of residual variance)
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
#'          - $var_additive: variance of predicted phenotype values (estimator of additive genetic variance)
#'          - $var_residual: variance of differnce between observed and predicted phenotype values (estimator of residual variance)
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
#'      + $GENOMIC_PREDICTIONS
#'          - $id: names of the samples/entries/pools which has genotype data but missing or outlying phenotype data
#'          - $pop: population or grouping the samples/entries/pools belong to
#'          - $y_true: missing observed phenotype data
#'          - $y_pred: predicted phenotype data using the best performing genomic prediction model
#'          - $model: best performing model based on within population cross-validation used to predict the missing phenotype data
#'      + $ADDITIVE_GENETIC_EFFECTS
#'          - $b: named numeric vector of SNP effects (or genotype effects for gBLUP), 
#'              where the names consist of tab-delimited chromosome/scaffold, position and if present, alleles.
#'          - $model: best performing genomic prediction model based on within population cross-validation 
#'              which yielded these effects.
#'  _ Err: gpError
#' @details
#' Functions list per module:
#'  1. io.R
#'      + 1.1. fn_G_extract_names
#'      + 1.2. fn_G_split_off_alternative_allele
#'      + 1.3. fn_G_numeric_to_non_numeric
#'      + 1.4. fn_G_non_numeric_to_numeric
#'      + 1.5. fn_G_to_vcf
#'      + 1.6. fn_vcf_to_G
#'      + 1.7. fn_classify_allele_frequencies
#'      + 1.8. fn_simulate_data
#'      + 1.9. fn_load_genotype
#'      + 1.10. fn_filter_genotype
#'      + 1.11. fn_save_genotype
#'      + 1.12. fn_load_phenotype
#'      + 1.13. fn_filter_phenotype
#'      + 1.14. fn_save_phenotype
#'      + 1.15. fn_merge_genotype_and_phenotype
#'      + 1.16. fn_subset_merged_genotype_and_phenotype
#'      + 1.17. fn_estimate_memory_footprint
#'  2. distances.R
#'      + 2.1. fn_grm_and_dist
#'  3. cross_validation.R
#'      + 3.1. fn_cv_1
#'      + 3.2. fn_cross_validation_preparation
#'      + 3.3. fn_cross_validation_within_population
#'      + 3.4. fn_cross_validation_across_populations_bulk
#'      + 3.5. fn_cross_validation_across_populations_pairwise
#'      + 3.6. fn_cross_validation_across_populations_lopo
#'  4. models.R
#'      + 4.1. fn_ols
#'      + 4.2. fn_ridge
#'      + 4.3. fn_lasso
#'      + 4.4. fn_elastic_net
#'      + 4.5. fn_Bayes_A
#'      + 4.6. fn_Bayes_B
#'      + 4.7. fn_Bayes_C
#'      + 4.8. fn_gBLUP
#'  5. metrics.R
#'      + 5.1. fn_prediction_performance_metrics
#' @examples
#' list_sim = fn_simulate_data(n=300, n_pop=3, verbose=TRUE)
#' df_pheno = read.delim(list_sim$fname_pheno_tsv, header=TRUE)
#' df_pheno$trait[which(df_pheno$pop=="pop_1")[1:3]] = NA
#' write.table(df_pheno, file=list_sim$fname_pheno_tsv,
#'     row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
#' args = list(
#'     fname_geno=list_sim$fname_geno_vcf,
#'     fname_pheno=list_sim$fname_pheno_tsv,
#'     population="pop_1",
#'     fname_covar=NULL,
#'     dir_output=NULL,
#'     geno_fname_snp_list=NULL,
#'     geno_ploidy=NULL,
#'     geno_bool_force_biallelic=TRUE,
#'     geno_bool_retain_minus_one_alleles_per_locus=TRUE,
#'     geno_min_depth=0,
#'     geno_max_depth=.Machine$integer.max,
#'     geno_maf=0.01,
#'     geno_sdev_min=0.0001,
#'     geno_max_n_alleles=NULL,
#'     geno_max_sparsity_per_locus=NULL,
#'     geno_frac_topmost_sparse_loci_to_remove=NULL,
#'     geno_n_topmost_sparse_loci_to_remove=NULL,
#'     geno_max_sparsity_per_sample=NULL,
#'     geno_frac_topmost_sparse_samples_to_remove=NULL,
#'     geno_n_topmost_sparse_samples_to_remove=NULL,
#'     pheno_sep="\t",
#'     pheno_header=TRUE,
#'     pheno_idx_col_id=1,
#'     pheno_idx_col_pop=2,
#'     pheno_idx_col_y=3,
#'     pheno_vec_na_strings=c("", "-", "NA", "na", "NaN", "missing", "MISSING"),
#'     pheno_bool_remove_outliers=TRUE,
#'     pheno_bool_remove_NA=FALSE,
#'     bool_within=TRUE,
#'     bool_across=TRUE,
#'     n_folds=2,
#'     n_reps=2,
#'     vec_models_to_test=c("ridge","lasso"),
#'     bool_parallel=TRUE,
#'     max_mem_Gb=15,
#'     n_threads=2,
#'     verbose=TRUE
#' )
#' fname_out_Rds = gp(args=args)
#' @export
gp = function(args) {
    ###################################################
    ### TEST
    # source("R/io.R")
    # source("R/metrics.R")
    # source("R/models.R")
    # source("R/cross_validation.R")
    # list_sim = fn_simulate_data(n=300, n_pop=3, verbose=TRUE)
    # df_pheno = read.delim(list_sim$fname_pheno_tsv, header=TRUE); df_pheno$trait[which(df_pheno$pop=="pop_1")[1:3]] = NA; write.table(df_pheno, file=list_sim$fname_pheno_tsv, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
    # args = list(
    #     fname_geno=list_sim$fname_geno_vcf,
    #     fname_pheno=list_sim$fname_pheno_tsv,
    #     population="pop_1",
    #     fname_covar=NULL,
    #     dir_output=NULL,
    #     geno_fname_snp_list=NULL,
    #     geno_ploidy=NULL,
    #     geno_bool_force_biallelic=TRUE,
    #     geno_bool_retain_minus_one_alleles_per_locus=TRUE,
    #     geno_min_depth=0,
    #     geno_max_depth=.Machine$integer.max,
    #     geno_maf=0.01,
    #     geno_sdev_min=0.0001,
    #     geno_max_n_alleles=NULL,
    #     geno_max_sparsity_per_locus=NULL,
    #     geno_frac_topmost_sparse_loci_to_remove=NULL,
    #     geno_n_topmost_sparse_loci_to_remove=NULL,
    #     geno_max_sparsity_per_sample=NULL,
    #     geno_frac_topmost_sparse_samples_to_remove=NULL,
    #     geno_n_topmost_sparse_samples_to_remove=NULL,
    #     pheno_sep="\t",
    #     pheno_header=TRUE,
    #     pheno_idx_col_id=1,
    #     pheno_idx_col_pop=2,
    #     pheno_idx_col_y=3,
    #     pheno_vec_na_strings=c("", "-", "NA", "na", "NaN", "missing", "MISSING"),
    #     pheno_bool_remove_outliers=FALSE,
    #     pheno_bool_remove_NA=FALSE,
    #     bool_within=TRUE,
    #     bool_across=TRUE,
    #     n_folds=2,
    #     n_reps=2,
    #     vec_models_to_test=c("ridge","lasso"),
    #     bool_parallel=TRUE,
    #     max_mem_Gb=15,
    #     n_threads=2,
    #     verbose=TRUE
    # )
    ###################################################
    ### Load genotype and phenotype data (TODO: covariate generation)
    G = fn_load_genotype(
        fname_geno=args$fname_geno,
        ploidy=args$geno_ploidy,
        force_biallelic=args$geno_bool_force_biallelic,
        retain_minus_one_alleles_per_locus=args$geno_bool_retain_minus_one_alleles_per_locus,
        min_depth=args$geno_min_depth,
        max_depth=args$geno_max_depth,
        verbose=args$verbose
    )
    if (methods::is(G, "gpError")) {return(G)}
    list_pheno = fn_load_phenotype(
        fname_pheno=args$fname_pheno,
        sep=args$pheno_sep,
        header=args$pheno_header,
        idx_col_id=args$pheno_idx_col_id,
        idx_col_pop=args$pheno_idx_col_pop,
        idx_col_y=args$pheno_idx_col_y,
        na_strings=args$pheno_vec_na_strings,
        verbose=args$verbose
    )
    if (methods::is(list_pheno, "gpError")) {return(list_pheno)}
    ##################################
    ### TODO: covariate generation ###
    ##################################
    ### Filter genotype and phenotype data
    G = fn_filter_genotype(
        G=G,
        maf=args$geno_maf,
        sdev_min=args$geno_sdev_min,
        fname_snp_list=args$geno_fname_snp_list,
        max_n_alleles=args$geno_max_n_alleles,
        max_sparsity_per_locus=args$geno_max_sparsity_per_locus,
        frac_topmost_sparse_loci_to_remove=args$geno_frac_topmost_sparse_loci_to_remove,
        n_topmost_sparse_loci_to_remove=args$geno_n_topmost_sparse_loci_to_remove,
        max_sparsity_per_sample=args$geno_max_sparsity_per_sample,
        frac_topmost_sparse_samples_to_remove=args$geno_frac_topmost_sparse_samples_to_remove,
        n_topmost_sparse_samples_to_remove=args$geno_n_topmost_sparse_samples_to_remove,
        verbose=args$verbose
    )
    if (methods::is(G, "gpError")) {
        error = chain(G, methods::new("gpError",
            code=100,
            message=paste0(
                "All loci were filtered out. ",
                "Please consider reducing --geno-min-depth (", args$geno_min_depth, ") and/or ",
                "increasing --geno-max-depth (", args$geno_max_depth, ") and/or ",
                "impute your input genotype data (", args$fname_geno, ")."
            )))
        return(error)
    }
    gc()
    list_pheno = fn_filter_phenotype(
        list_pheno=list_pheno,
        remove_outliers=args$pheno_bool_remove_outliers,
        remove_NA=args$pheno_bool_remove_NA,
        verbose=args$verbose
    )
    if (methods::is(list_pheno, "gpError")) {return(list_pheno)}
    gc()
    ### Merge filtered genotype, phenotype and covariate data
    list_merged = fn_merge_genotype_and_phenotype(
        G=G,
        list_pheno=list_pheno,
        COVAR=NULL,
        verbose=args$verbose
    )
    if (methods::is(list_merged, "gpError")) {return(list_merged)}
    ### Missing values are not allowed in the genotype data
    if (sum(rowSums(is.na(list_merged$G))) > 0) {
        error = methods::new("gpError",
            code=101,
            message=paste0(
                "Genotype data has missing values. ",
                "Please impute the missing data. ",
                "You may also consider reducing ", 
                "--geno-min-depth (", args$geno_min_depth, ") and/or increasing ",
                "--geno-max-depth (", args$geno_max_depth, "), if possible. ",
                "You may also set the following to zero: --geno-max-sparsity-per-locus (", args$geno_max_sparsity_per_locus, ") and ",
                "--geno-max-sparsity-per-sample (", args$geno_max_sparsity_per_sample, ")."
            ))
        return(error)
    }
    ### Clean-up
    rm("G")
    rm("list_pheno")
    gc(); gc(); gc(); gc(); gc(); gc(); gc(); gc()
    ### Extract the trait name
    trait_name = list_merged$list_pheno$trait_name
    ### Start with across population cross-validation using the full dataset
    if (args$bool_across) {
        ################################
        ### ACROSS POPULATIONS: BULK ###
        ################################
        fname_across_bulk_Rds = fn_cross_validation_across_populations_bulk(
            list_merged=list_merged,
            n_folds=args$n_folds,
            n_reps=args$n_reps,
            vec_models_to_test=args$vec_models_to_test,
            bool_parallel=args$bool_parallel,
            max_mem_Gb=args$max_mem_Gb,
            n_threads=args$n_threads,
            dir_output=args$dir_output,
            verbose=args$verbose
        )
        if (methods::is(fname_across_bulk_Rds, "gpError")) {
            return(fname_across_bulk_Rds)
        } else {
            list_across_bulk = readRDS(fname_across_bulk_Rds)
            METRICS_ACROSS_POP_BULK = list_across_bulk$METRICS_ACROSS_POP_BULK
            YPRED_ACROSS_POP_BULK = list_across_bulk$YPRED_ACROSS_POP_BULK
        }
        ########################################
        ### ACROSS POPULATIONS: PAIRWISE-POP ###
        ########################################
        fname_across_pairwise_Rds = fn_cross_validation_across_populations_pairwise(
            list_merged=list_merged,
            vec_models_to_test=args$vec_models_to_test,
            bool_parallel=args$bool_parallel,
            max_mem_Gb=args$max_mem_Gb,
            n_threads=args$n_threads,
            dir_output=args$dir_output,
            verbose=args$verbose
        )
        if (methods::is(fname_across_pairwise_Rds, "gpError")) {
            return(fname_across_pairwise_Rds)
        } else {
            list_across_pairwise = readRDS(fname_across_pairwise_Rds)
            METRICS_ACROSS_POP_PAIRWISE = list_across_pairwise$METRICS_ACROSS_POP_PAIRWISE
            YPRED_ACROSS_POP_PAIRWISE = list_across_pairwise$YPRED_ACROSS_POP_PAIRWISE
        }
        #############################################
        ### ACROSS POPULATIONS: LEAVE-ONE-POP-OUT ###
        #############################################
        fname_across_lopo_Rds = fn_cross_validation_across_populations_lopo(
            list_merged=list_merged,
            vec_models_to_test=args$vec_models_to_test,
            bool_parallel=args$bool_parallel,
            max_mem_Gb=args$max_mem_Gb,
            n_threads=args$n_threads,
            dir_output=args$dir_output,
            verbose=args$verbose
        )
        if (methods::is(fname_across_lopo_Rds, "gpError")) {
            return(fname_across_lopo_Rds)
        } else {
            list_across_lopo = readRDS(fname_across_lopo_Rds)
            METRICS_ACROSS_POP_LOPO = list_across_lopo$METRICS_ACROSS_POP_LOPO
            YPRED_ACROSS_POP_LOPO = list_across_lopo$YPRED_ACROSS_POP_LOPO
        }
    } else {
        METRICS_ACROSS_POP_BULK = NA
        YPRED_ACROSS_POP_BULK = NA
        METRICS_ACROSS_POP_LOPO = NA
        YPRED_ACROSS_POP_LOPO = NA
        METRICS_ACROSS_POP_PAIRWISE = NA
        YPRED_ACROSS_POP_PAIRWISE = NA
    }
    ### Finish-off with within population cross-validation
    ### Subset by a single population and clean-up
    if (args$bool_within) {
        list_merged = fn_subset_merged_genotype_and_phenotype(
            list_merged=list_merged,
            vec_idx=which(list_merged$list_pheno$pop==args$population),
            verbose=args$verbose
        )
        if (methods::is(list_merged, "gpError")) {return(list_merged)}
        gc(); gc(); gc(); gc(); gc(); gc(); gc(); gc()
        ### Filter subsetted genotype and phenotype data
        G = fn_filter_genotype(
            G=list_merged$G,
            maf=args$geno_maf,
            sdev_min=args$geno_sdev_min,
            fname_snp_list=args$geno_fname_snp_list,
            max_n_alleles=args$geno_max_n_alleles,
            max_sparsity_per_locus=args$geno_max_sparsity_per_locus,
            frac_topmost_sparse_loci_to_remove=args$geno_frac_topmost_sparse_loci_to_remove,
            n_topmost_sparse_loci_to_remove=args$geno_n_topmost_sparse_loci_to_remove,
            max_sparsity_per_sample=args$geno_max_sparsity_per_sample,
            frac_topmost_sparse_samples_to_remove=args$geno_frac_topmost_sparse_samples_to_remove,
            n_topmost_sparse_samples_to_remove=args$geno_n_topmost_sparse_samples_to_remove,
            verbose=args$verbose
        )
        if (methods::is(G, "gpError")) {return(G)}
        gc()
        list_pheno = fn_filter_phenotype(
            list_pheno=list_merged$list_pheno,
            remove_outliers=args$pheno_bool_remove_outliers,
            remove_NA=args$pheno_bool_remove_NA,
            verbose=args$verbose
        )
        if (methods::is(list_pheno, "gpError")) {return(list_pheno)}
        gc()
        ### Merge filtered genotype, phenotype and covariate data
        list_merged = fn_merge_genotype_and_phenotype(
            G=G,
            list_pheno=list_pheno,
            COVAR=NULL,
            verbose=args$verbose
        )
        if (methods::is(list_merged, "gpError")) {return(list_merged)}
        ### Clean-up
        rm("G")
        rm("list_pheno")
        gc(); gc(); gc(); gc(); gc(); gc(); gc(); gc()
    }
    ### Genomic prediction cross-validation
    if (args$bool_within) {
        #########################
        ### WITHIN POPULATION ###
        #########################
        fname_within_Rds = fn_cross_validation_within_population(
            list_merged=list_merged,
            n_folds=args$n_folds,
            n_reps=args$n_reps,
            vec_models_to_test=args$vec_models_to_test,
            bool_parallel=args$bool_parallel,
            max_mem_Gb=args$max_mem_Gb,
            n_threads=args$n_threads,
            dir_output=args$dir_output,
            verbose=args$verbose
        )
        if (methods::is(fname_within_Rds, "gpError")) {
            return(fname_within_Rds)
        } else {
            list_within = readRDS(fname_within_Rds)
            METRICS_WITHIN_POP = list_within$METRICS_WITHIN_POP
            YPRED_WITHIN_POP = list_within$YPRED_WITHIN_POP
        }
    } else {
        METRICS_WITHIN_POP = NA
        YPRED_WITHIN_POP = NA
    }
    ##################################
    ### GENOMIC PREDICTIONS PER SE ###
    ##################################
    if (args$bool_within) {
        ### Are there any genotypes with missing phenotype data?
        vec_idx_validation = which(is.na(list_merged$list_pheno$y))
        vec_idx_training = which(!is.na(list_merged$list_pheno$y))
        if (length(vec_idx_validation)==0) {
            GENOMIC_PREDICTIONS = NA
        } else {
            ### Find the best model in the args$population
            df_agg = stats::aggregate(corr ~ model, data=METRICS_WITHIN_POP, FUN=mean, na.rm=TRUE)
            idx = which(df_agg$corr == max(df_agg$corr, na.rm=TRUE))[1]
            model = df_agg$model[idx]
            ### Define additional model input/s
            if (grepl("Bayes", model)==TRUE) {
                ### Append the input prefix into the temporary file prefix generated by Bayesian models so we don't overwite these when performing parallel computations
                if (is.null(args$dir_output)) {
                    args$dir_output = dirname(tempfile())
                }
                other_params = list(nIter=12e3, burnIn=2e3, out_prefix=file.path(args$dir_output, paste0("bglr_", model)))
            } else {
                other_params = list(n_folds=10)
            }
            perf = eval(parse(text=paste0("fn_", model, "(list_merged=list_merged, vec_idx_training=vec_idx_training, vec_idx_validation=vec_idx_validation, other_params=other_params, verbose=args$verbose)")))
            GENOMIC_PREDICTIONS = data.frame(perf$df_y_validation, model=model)
        }
    } else {
        GENOMIC_PREDICTIONS = NA
    }
    ###############################################################################
    ### EXTRACT WITHIN POPULATION ADDITIVE GENETIC EFFECTS OF EACH MODEL TESTED ###
    ###############################################################################
    ### Extract the additive genetic effects from penalised and Bayesian regression models, and sample BLUPs from gBLUP
    if (args$bool_within) {
        vec_idx_training = which(!is.na(list_merged$list_pheno$y))
        ADDITIVE_GENETIC_EFFECTS = NULL
        for (model in args$vec_models_to_test) {
            # model = args$vec_models_to_test[1]
            if ((grepl("Bayes", model)==TRUE) | (grepl("gBLUP", model)==TRUE)) {
                if (is.null(args$dir_output)) {
                    args$dir_output = dirname(tempfile())
                }
                vec_idx_validation = c()
                other_params = list(nIter=12e3, burnIn=2e3, out_prefix=file.path(args$dir_output, paste0("bglr_", model)))
            } else {
                vec_idx_validation = vec_idx_training
                other_params = list(n_folds=10)
            }
            perf = eval(parse(text=paste0("fn_", model, "(list_merged=list_merged, vec_idx_training=vec_idx_training, vec_idx_validation=vec_idx_validation, other_params=other_params, verbose=FALSE)")))
            if (is.null(ADDITIVE_GENETIC_EFFECTS)) {
                ADDITIVE_GENETIC_EFFECTS = eval(parse(text=paste0("list(`", model, "`=list(b=perf$vec_effects))")))
            } else {
                eval(parse(text=paste0("ADDITIVE_GENETIC_EFFECTS$`", model, "`=list(b=perf$vec_effects)")))
            }
        }
    } else {
        ADDITIVE_GENETIC_EFFECTS = NA
    }
    ### Clean-up
    rm("list_merged")
    gc()
    ### Output
    if (is.null(args$dir_output)) {
        args$dir_output = dirname(tempfile())
    }
    time_rand_id = paste0(round(as.numeric(Sys.time())), sample.int(1e6, size=1))
    fname_out_Rds = file.path(args$dir_output, paste0("GENOMIC_PREDICTIONS_OUTPUT-", trait_name, "-", args$population, "-", time_rand_id, ".Rds"))
    saveRDS(list(
        TRAIT_NAME=trait_name,
        POPULATION=args$population,
        METRICS_WITHIN_POP=METRICS_WITHIN_POP,
        YPRED_WITHIN_POP=YPRED_WITHIN_POP,
        METRICS_ACROSS_POP_BULK=METRICS_ACROSS_POP_BULK,
        YPRED_ACROSS_POP_BULK=YPRED_ACROSS_POP_BULK,
        METRICS_ACROSS_POP_PAIRWISE=METRICS_ACROSS_POP_PAIRWISE,
        YPRED_ACROSS_POP_PAIRWISE=YPRED_ACROSS_POP_PAIRWISE,
        METRICS_ACROSS_POP_LOPO=METRICS_ACROSS_POP_LOPO,
        YPRED_ACROSS_POP_LOPO=YPRED_ACROSS_POP_LOPO,
        GENOMIC_PREDICTIONS=GENOMIC_PREDICTIONS,
        ADDITIVE_GENETIC_EFFECTS=ADDITIVE_GENETIC_EFFECTS
    ), file=fname_out_Rds)
    ### Clean-up
    if (exists("fname_within_Rds")) {unlink(fname_within_Rds)}
    if (exists("fname_across_bulk_Rds")) {unlink(fname_across_bulk_Rds)}
    if (exists("fname_across_pairwise_Rds")) {unlink(fname_across_pairwise_Rds)}
    if (exists("fname_across_lopo_Rds")) {unlink(fname_across_lopo_Rds)}
    ### Return the filename of the output list saved as an Rds file    
    return(fname_out_Rds)
}

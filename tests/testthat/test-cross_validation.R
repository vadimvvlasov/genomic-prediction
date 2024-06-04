# library(testthat)
# source("R/io.R")
# source("R/metrics.R")
# source("R/models.R")
# source("R/cross_validation.R")

test_that("fn_kfold_cross_validation", {
    set.seed(123)
    list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    COVAR = G %*% t(G)
    n = nrow(G)
    vec_models_to_test = c("ridge", "lasso", "elastic_net", "Bayes_A", "Bayes_B", "Bayes_C", "gBLUP")
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    mat_idx_shuffle = cbind(c(1:n), sample(1:n, size=n, replace=FALSE))
    vec_set_partition_groupings = kmeans(G, centers=5)$cluster
    n_reps = ncol(mat_idx_shuffle)
    k_folds = max(vec_set_partition_groupings)
    df_params = expand.grid(rep=c(1:n_reps), fold=c(1:k_folds), model=vec_models_to_test)
    i_clusteredFolds = which((df_params$rep == 1) & (df_params$model == "ridge"))[1]
    i_randomFolds = which((df_params$rep == 2) & (df_params$model == "ridge"))[1]
    prefix_tmp = tempfile()
    list_cv_clusteredFolds = fn_cv_1(
        list_merged=list_merged, 
        i=i_clusteredFolds, 
        df_params=df_params, 
        mat_idx_shuffle=mat_idx_shuffle, 
        vec_set_partition_groupings=vec_set_partition_groupings,
        prefix_tmp=prefix_tmp,
        verbose=TRUE)
    list_cv_randomFolds = fn_cv_1(
        list_merged=list_merged, 
        i=i_randomFolds, 
        df_params=df_params, 
        mat_idx_shuffle=mat_idx_shuffle, 
        vec_set_partition_groupings=vec_set_partition_groupings,
        prefix_tmp=prefix_tmp,
        verbose=TRUE)
    expect_equal(list_cv_randomFolds$df_metrics$corr > list_cv_clusteredFolds$df_metrics$corr, TRUE)
})
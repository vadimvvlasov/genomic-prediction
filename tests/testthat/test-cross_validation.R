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

test_that("fn_cross_validation_preparation", {
    set.seed(123)
    list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    COVAR = G %*% t(G)
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    list_cv_params_1 = fn_cross_validation_preparation(list_merged, cv_type=1)
    list_cv_params_2_err = fn_cross_validation_preparation(list_merged, cv_type=2)
    list_cv_params_3 = fn_cross_validation_preparation(list_merged, cv_type=3)
    expect_equal(nrow(list_cv_params_1$df_params), 7*10*10)
    expect_equal(sum(dim(list_cv_params_1$mat_idx_shuffle) == c(100, 10)), 2)
    expect_equal(sum(list_cv_params_1$vec_set_partition_groupings == rep(c(1:10), each=10)), 100)
    expect_equal(methods::is(list_cv_params_2_err, "gpError"), TRUE)
    expect_equal(nrow(list_cv_params_3$df_params), 7*3*1)
    expect_equal(sum(dim(list_cv_params_3$mat_idx_shuffle) == c(100, 1)), 2)
    expect_equal(sum(list_cv_params_3$vec_set_partition_groupings == as.numeric(as.factor(list_merged$list_pheno$pop))), 100)
    ### Force 2 populations only
    list_merged$list_pheno$pop = rep(c("popA", "popB"), times=length(list_merged$list_pheno$pop)/2)
    list_cv_params_2 = fn_cross_validation_preparation(list_merged, cv_type=2)
    list_cv_params_3_paired = fn_cross_validation_preparation(list_merged, cv_type=2)
    expect_equal(nrow(list_cv_params_2$df_params), 7*2*1)
    expect_equal(sum(dim(list_cv_params_2$mat_idx_shuffle) == c(100, 1)), 2)
    expect_equal(sum(list_cv_params_2$vec_set_partition_groupings == as.numeric(as.factor(list_merged$list_pheno$pop))), 100)
    expect_equal(sum(list_cv_params_2$df_params == list_cv_params_3_paired$df_params), prod(dim(list_cv_params_2$df_params)))
    expect_equal(sum(list_cv_params_2$mat_idx_shuffle == list_cv_params_3_paired$mat_idx_shuffle), prod(dim(list_cv_params_2$mat_idx_shuffle)))
    expect_equal(sum(list_cv_params_2$vec_set_partition_groupings == list_cv_params_3_paired$vec_set_partition_groupings), length(list_cv_params_2$vec_set_partition_groupings))
})

test_that("fn_cross_validation_within_population", {
    set.seed(123)
    list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    COVAR = G %*% t(G)
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    fname_within_Rds = fn_cross_validation_within_population(list_merged, n_folds=2, n_reps=1, vec_models_to_test=c("ridge","lasso"), verbose=TRUE)
    list_within = readRDS(fname_within_Rds)
    expect_equal(sum(dim(list_within$METRICS_WITHIN_POP) == c(3*2*1*2, 21)), 2)
    expect_equal(sum(dim(list_within$YPRED_WITHIN_POP) == c(100*2, 8)), 2)
    expect_equal(mean(list_within$METRICS_WITHIN_POP$corr) < 0.5, TRUE)
    expect_equal(cor(list_within$YPRED_WITHIN_POP$y_true, list_within$YPRED_WITHIN_POP$y_pred) < 0.5, TRUE)
})

test_that("fn_cross_validation_across_populations_bulk", {
    set.seed(123)
    list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    COVAR = G %*% t(G)
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    fname_across_bulk_Rds = fn_cross_validation_across_populations_bulk(list_merged, n_folds=2, n_reps=1, vec_models_to_test=c("ridge","lasso"), verbose=TRUE)
    list_across_bulk = readRDS(fname_across_bulk_Rds)
    expect_equal(sum(dim(list_across_bulk$METRICS_ACROSS_POP_BULK) == c(1*2*1*2, 21)), 2)
    expect_equal(sum(dim(list_across_bulk$YPRED_ACROSS_POP_BULK) == c(100*2, 8)), 2)
    expect_equal(mean(list_across_bulk$METRICS_ACROSS_POP_BULK$corr) < 0.5, TRUE)
    expect_equal(cor(list_across_bulk$YPRED_ACROSS_POP_BULK$y_true, list_across_bulk$YPRED_ACROSS_POP_BULK$y_pred) < 0.5, TRUE)
})

test_that("fn_cross_validation_across_populations_pairwise", {
    set.seed(123)
    list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    COVAR = G %*% t(G)
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    fname_across_pairwise_Rds = fn_cross_validation_across_populations_pairwise(list_merged, vec_models_to_test=c("ridge","lasso"), verbose=TRUE)
    list_across_pairwise = readRDS(fname_across_pairwise_Rds)
    expect_equal(sum(dim(list_across_pairwise$METRICS_ACROSS_POP_PAIRWISE) == c((3*(3-1))*1*2, 21)), 2)
    expect_equal(sum(dim(list_across_pairwise$YPRED_ACROSS_POP_PAIRWISE) == c((3*(3-1))*(100/3)*2, 8)), 2)
    expect_equal(mean(list_across_pairwise$METRICS_ACROSS_POP_PAIRWISE$corr) < 0.5, TRUE)
    expect_equal(cor(list_across_pairwise$YPRED_ACROSS_POP_PAIRWISE$y_true, list_across_pairwise$YPRED_ACROSS_POP_PAIRWISE$y_pred) < 0.5, TRUE)
})

test_that("fn_cross_validation_across_populations_lopo", {
    set.seed(123)
    list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    COVAR = G %*% t(G)
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    fname_across_lopo_Rds = fn_cross_validation_across_populations_lopo(list_merged, vec_models_to_test=c("ridge","lasso"), verbose=TRUE)
    list_across_lopo = readRDS(fname_across_lopo_Rds)
    expect_equal(sum(dim(list_across_lopo$METRICS_ACROSS_POP_LOPO) == c(3*1*2, 21)), 2)
    expect_equal(sum(dim(list_across_lopo$YPRED_ACROSS_POP_LOPO) == c(3*(100/3)*2, 8)), 2)
    expect_equal(mean(list_across_lopo$METRICS_ACROSS_POP_LOPO$corr) < 0.9, TRUE)
    expect_equal(cor(list_across_lopo$YPRED_ACROSS_POP_LOPO$y_true, list_across_lopo$YPRED_ACROSS_POP_LOPO$y_pred) < 0.9, TRUE)
})

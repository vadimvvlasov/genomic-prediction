# library(testthat)
# source("/group/pasture/Jeff/gp/R/load.R")
# source("/group/pasture/Jeff/gp/R/metrics.R")
# source("/group/pasture/Jeff/gp/R/models.R")

test_that("fn_ols", {
    set.seed(123)
    list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
    n = nrow(list_merged$G)
    vec_idx_training = sample(c(1:n), floor(n/2))
    vec_idx_validation = c(1:n)[!(c(1:n) %in% vec_idx_training)]
    list_ols = fn_ols(list_merged, vec_idx_training, vec_idx_validation, verbose=TRUE)
    expect_equal(list_ols$list_perf$corr < 0.5, TRUE)
    expect_equal(list_ols$list_perf$corr, cor(list_ols$df_y_validation$y_true, list_ols$df_y_validation$y_pred))
    expect_equal(nrow(list_ols$df_y_validation), length(vec_idx_validation))
    expect_equal(list_ols$n_non_zero, ncol(G))
})


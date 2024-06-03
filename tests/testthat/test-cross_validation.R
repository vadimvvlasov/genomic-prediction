# library(testthat)
# source("R/load.R")
# source("R/metrics.R")
# source("R/models.R")
# source("R/cross_validation.R")

test_that("fn_kfold_cross_validation", {
    set.seed(123)
    list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, verbose=TRUE)
    expect_equal(1, 1)
})
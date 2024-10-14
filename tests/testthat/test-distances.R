# library(testthat)
# source("R/io.R")
# source("R/distances.R")

test_that("fn_G_extract_names", {
    set.seed(123)
    list_sim = fn_simulate_data(verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_grm_and_dist = fn_grm_and_dist(G=G, verbose=TRUE)
    expect_equal(length(list_grm_and_dist), 10)
    expect_equal(unique(unlist(lapply(list_grm_and_dist, FUN=function(x){dim(x)}))), 100)
})

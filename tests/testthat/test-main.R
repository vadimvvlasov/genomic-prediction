# library(testthat)
# devtools::load_all()

test_that("gp", {
    set.seed(123)
    list_sim = fn_simulate_data(n=300, n_pop=3, verbose=TRUE)
    df_pheno = read.delim(list_sim$fname_pheno_tsv, header=TRUE); df_pheno$trait[which(df_pheno$pop=="pop_1")[1:3]] = NA; write.table(df_pheno, file=list_sim$fname_pheno_tsv, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
    args = list(
        fname_geno=list_sim$fname_geno_vcf,
        fname_pheno=list_sim$fname_pheno_tsv,
        population="pop_1",
        fname_covar=NULL,
        dir_output=NULL,
        geno_fname_snp_list=NULL,
        geno_ploidy=NULL,
        geno_bool_force_biallelic=TRUE,
        geno_bool_retain_minus_one_alleles_per_locus=TRUE,
        geno_min_depth=0,
        geno_max_depth=.Machine$integer.max,
        geno_maf=0.01,
        geno_sdev_min=0.0001,
        geno_max_n_alleles=NULL,
        geno_max_sparsity_per_locus=NULL,
        geno_frac_topmost_sparse_loci_to_remove=NULL,
        geno_n_topmost_sparse_loci_to_remove=NULL,
        geno_max_sparsity_per_sample=NULL,
        geno_frac_topmost_sparse_samples_to_remove=NULL,
        geno_n_topmost_sparse_samples_to_remove=NULL,
        pheno_sep="\t",
        pheno_header=TRUE,
        pheno_idx_col_id=1,
        pheno_idx_col_pop=2,
        pheno_idx_col_y=3,
        pheno_na_strings=c("", "-", "NA", "na", "NaN", "missing", "MISSING"),
        pheno_bool_remove_outliers=TRUE,
        pheno_bool_remove_NA=FALSE,
        bool_within=TRUE,
        bool_across=TRUE,
        n_folds=2,
        n_reps=2,
        vec_models_to_test=c("ridge","lasso"),
        bool_parallel=TRUE,
        max_mem_Gb=15,
        n_threads=2,
        verbose=TRUE
    )
    fname_out_Rds = gp(args=args)
    expect_equal(methods::is(fname_out_Rds, "gpError"), FALSE)
})

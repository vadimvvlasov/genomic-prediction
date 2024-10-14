# library(testthat)
# source("R/io.R")

test_that("fn_G_extract_names", {
    set.seed(123)
    n = 100
    n_alleles = 3
    mat_genotypes = simquantgen::fn_simulate_genotypes(n=n, n_alleles=n_alleles, verbose=TRUE)
    list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=mat_genotypes, verbose=TRUE)
    ### Assertions
    n_digits = length(unlist(strsplit(as.character(n), "")))
    bool_test_entry_names = list_ids_chr_pos_all$vec_ids == paste0("entry_", sprintf(paste0("%0", n_digits, "d"), 1:n))
    expect_equal(sum(bool_test_entry_names), n)
    expect_equal(length(list_ids_chr_pos_all$vec_chr), ncol(mat_genotypes))
    expect_equal(length(list_ids_chr_pos_all$vec_pos), ncol(mat_genotypes))
    expect_equal(length(list_ids_chr_pos_all$vec_all), ncol(mat_genotypes))
    expect_equal(length(unique(list_ids_chr_pos_all$vec_all)), n_alleles-1)
})

test_that("fn_G_split_off_alternative_allele", {
    set.seed(123)
    G_ref = simquantgen::fn_simulate_genotypes(verbose=TRUE)
    G_alt = 1 - G_ref; colnames(G_alt) = gsub("allele_1$", "allele_2", colnames(G_alt))
    G_refalt = cbind(G_ref, G_alt)
    list_G_G_alt = fn_G_split_off_alternative_allele(G=G_refalt, verbose=TRUE)
    expect_equal(sum(G_ref == list_G_G_alt$G), prod(dim(G_ref)))
    expect_equal(sum(G_alt == list_G_G_alt$G_alt), prod(dim(G_alt)))
})

test_that("fn_G_numeric_to_non_numeric", {
    set.seed(123)
    ploidy = 42
    G_numeric = simquantgen::fn_simulate_genotypes(ploidy=ploidy, n_alleles=52, verbose=TRUE)
    G_non_numeric = fn_G_numeric_to_non_numeric(G=G_numeric, ploidy=ploidy, verbose=TRUE)
    vec_counts = table(unlist(strsplit(G_non_numeric[1,1], ""))); names(vec_counts) = NULL
    expect_equal(G_numeric[1,1], vec_counts[1] / sum(vec_counts))
    expect_equal(sum(vec_counts), ploidy)
})

test_that("fn_G_non_numeric_to_numeric", {
    set.seed(123)
    ploidy = 42
    n_alleles = 2
    G_numeric = simquantgen::fn_simulate_genotypes(ploidy=ploidy, n_alleles=n_alleles, verbose=TRUE)
    G_non_numeric = fn_G_numeric_to_non_numeric(G=G_numeric, ploidy=ploidy, verbose=TRUE)
    G_numeric_back = fn_G_non_numeric_to_numeric(G_non_numeric=G_non_numeric, verbose=TRUE)
    expect_equal(sum(abs(G_numeric - G_numeric_back) < 1e-4), prod(dim(G_numeric)))
    ### The converted non-numeric to numeric matrix can have less loci-alleles than the original numeric matrix as fixed loci will be omitted
    expect_equal(ncol(G_numeric_back) <= ncol(G_numeric), TRUE)
})

test_that("fn_G_to_vcf", {
    set.seed(123)
    n = 123
    l = 456
    n_alleles = 2
    G = simquantgen::fn_simulate_genotypes(n=n, l=l, n_alleles=n_alleles, verbose=TRUE)
    vcf = fn_G_to_vcf(G, verbose=TRUE)
    expect_equal(dim(vcf@gt), c(l, n+1))
    ### Error catching
    G_triallelic = simquantgen::fn_simulate_genotypes(n=n, l=l, n_alleles=3, verbose=TRUE)
    vcf_error = fn_G_to_vcf(G_triallelic, verbose=TRUE)
    expect_equal(class(vcf_error)[1], "gpError")
})

test_that("fn_vcf_to_G", {
    set.seed(123)
    G = simquantgen::fn_simulate_genotypes(verbose=TRUE)
    vcf = fn_G_to_vcf(G=G, min_depth=1000, max_depth=1000, verbose=TRUE)
    list_G_D = fn_vcf_to_G(vcf=vcf, verbose=TRUE)
    G_back = list_G_D$G
    expect_equal(sum(colnames(G) == colnames(G_back)), ncol(G))
    expect_equal(sum(rownames(G) == rownames(G_back)), nrow(G))
    expect_equal(sum(abs(G_back-G) < 1e-7), prod(dim(G)))
})

test_that("fn_classify_allele_frequencies", {
    set.seed(123)
    ploidy = 4
    G = simquantgen::fn_simulate_genotypes(ploidy=ploidy, verbose=TRUE)
    G_classes = fn_classify_allele_frequencies(G=G, ploidy=ploidy, verbose=TRUE)
    G_classes_diploid = fn_classify_allele_frequencies(G=G, ploidy=2, verbose=TRUE)
    expect_equal(sum(G == G_classes), prod(dim(G)))
    expect_equal(length(unique(as.vector(G_classes))), ploidy+1)
    expect_equal(length(unique(as.vector(G_classes_diploid))), 2+1)
})

test_that("fn_simulate_data", {
    set.seed(123)
    list_sim = fn_simulate_data(verbose=TRUE)
    expect_equal(is.null(list_sim$fname_geno_vcf), FALSE)
    expect_equal(is.null(list_sim$fname_geno_tsv), TRUE)
    expect_equal(is.null(list_sim$fname_geno_rds), TRUE)
    expect_equal(is.null(list_sim$fname_pheno_tsv), FALSE)
    unlink(list_sim$fname_geno_vcf)
    unlink(list_sim$fname_pheno_tsv)
    list_sim = fn_simulate_data(min_depth=1000, max_depth=1000, save_geno_vcf=TRUE, save_geno_rds=TRUE, save_geno_tsv=TRUE, save_pheno_tsv=TRUE, verbose=TRUE)
    list_G_D = fn_vcf_to_G(vcf=vcfR::read.vcfR(list_sim$fname_geno_vcf))
    G_vcf = list_G_D$G
    df_tsv = utils::read.delim(list_sim$fname_geno_tsv, sep="\t", header=TRUE)
    G_tsv = as.matrix(t(df_tsv[, c(-1,-2,-3)])); rownames(G_tsv) = colnames(df_tsv)[c(-1,-2,-3)]; colnames(G_tsv) = paste(df_tsv$chr, df_tsv$pos, df_tsv$allele, sep="\t")
    G_rds = readRDS(list_sim$fname_geno_rds)
    df_pheno = utils::read.delim(list_sim$fname_pheno_tsv, sep="\t", header=TRUE)
    expect_equal(sum(abs(G_vcf - G_tsv) < 1e-7), prod(dim(G_vcf)))
    expect_equal(sum(abs(G_vcf - G_rds) < 1e-7), prod(dim(G_vcf)))
    expect_equal(sum(abs(G_rds - G_tsv) < 1e-7), prod(dim(G_vcf)))
    expect_equal(sum(df_pheno$id == rownames(G_vcf)), nrow(df_pheno))
    expect_equal(sum(df_pheno$id == rownames(G_tsv)), nrow(df_pheno))
    expect_equal(sum(df_pheno$id == rownames(G_rds)), nrow(df_pheno))
    unlink(list_sim$fname_geno_vcf)
    unlink(list_sim$fname_geno_tsv)
    unlink(list_sim$fname_geno_rds)
    unlink(list_sim$fname_pheno_tsv)
    list_sim = fn_simulate_data(save_geno_vcf=FALSE, save_geno_rds=TRUE, non_numeric_Rds=TRUE, verbose=TRUE)
    G_non_numeric = readRDS(list_sim$fname_geno_rds)
    expect_equal(is.numeric(G_non_numeric), FALSE)
    unlink(list_sim$fname_geno_rds)
    unlink(list_sim$fname_pheno_tsv)
})

test_that("fn_load_genotype", {
    set.seed(123)
    list_sim = fn_simulate_data(min_depth=1000, max_depth=1000, save_geno_vcf=TRUE, save_geno_tsv=TRUE, save_geno_rds=TRUE, save_pheno_tsv=TRUE, verbose=TRUE)
    G_vcf = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf, verbose=TRUE)
    G_tsv = fn_load_genotype(fname_geno=list_sim$fname_geno_tsv, verbose=TRUE)
    G_rds = fn_load_genotype(fname_geno=list_sim$fname_geno_rds, verbose=TRUE)
    expect_equal(sum(abs(G_vcf - G_tsv) < 1e-7), prod(dim(G_vcf)))
    expect_equal(sum(abs(G_vcf - G_rds) < 1e-7), prod(dim(G_vcf)))
    expect_equal(sum(abs(G_rds - G_tsv) < 1e-7), prod(dim(G_vcf)))
    unlink(list_sim$fname_geno_vcf)
    unlink(list_sim$fname_geno_tsv)
    unlink(list_sim$fname_geno_rds)
    unlink(list_sim$fname_pheno_tsv)
})

test_that("fn_filter_genotype", {
    set.seed(123)
    list_sim = fn_simulate_data(verbose=TRUE)
    maf = 0.05
    sdev_min = 0.0001
    max_sparsity_per_locus = 0.4
    frac_topmost_sparse_loci_to_remove = 0.01
    n_topmost_sparse_loci_to_remove = 100
    max_sparsity_per_sample = 0.3
    frac_topmost_sparse_samples_to_remove = 0.01
    n_topmost_sparse_samples_to_remove = 10
    verbose = TRUE
    ### Do not load the alternative alleles
    G = fn_load_genotype(list_sim$fname_geno_vcf)
    G_filtered_1 = fn_filter_genotype(G=G, verbose=TRUE)
    expect_equal(sum(dim(G) == dim(G_filtered_1)), 2)
    ### Simulate SNP list for filtering
    G = fn_load_genotype(list_sim$fname_geno_vcf, retain_minus_one_alleles_per_locus=FALSE)
    colnames(G) = gsub("allele_1", "A", colnames(G)) ### Rename allele_1 and allele_alt to A and T, respectively to allow filtering using a SNP list
    colnames(G) = gsub("allele_alt", "T", colnames(G)) ### Rename allele_1 and allele_alt to A and T, respectively to allow filtering using a SNP list
    n_sim_missing = 100
    mat_loci = matrix(unlist(strsplit(colnames(G), "\t")), byrow=TRUE, ncol=3)
    vec_loci = unique(paste0(mat_loci[,1], "\t", mat_loci[,2]))
    mat_loci = matrix(unlist(strsplit(vec_loci, "\t")), byrow=TRUE, ncol=2)
    df_snp_list = data.frame(CHROM=mat_loci[,1], POS=as.numeric(mat_loci[,2]), REF_ALT=paste0("A,T"))
    df_snp_list$REF_ALT[1:n_sim_missing] = "C,G"
    colnames(df_snp_list) = c("#CHROM", "POS", "REF,ALT")
    fname_snp_list = tempfile(fileext=".snplist")
    utils::write.table(df_snp_list, file=fname_snp_list, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    ### Filter with the alternative alleles
    G_filtered_2 = fn_filter_genotype(G=G, fname_snp_list=fname_snp_list, verbose=TRUE)
    expect_equal(ncol(G_filtered_2), ncol(G) - (2*n_sim_missing))
    ### Filter without the alternative alleles
    G = fn_load_genotype(list_sim$fname_geno_vcf)
    colnames(G) = gsub("allele_1", "A", colnames(G)) ### Rename allele_1 to A to allow filtering using a SNP list
    G_filtered_3 = fn_filter_genotype(G=G, maf=0.05, fname_snp_list=fname_snp_list, verbose=TRUE)
    expect_equal(ncol(G_filtered_3), ncol(G) - n_sim_missing)
    G_filtered_2_split_G = fn_G_split_off_alternative_allele(G=G_filtered_2, verbose=TRUE)$G
    expect_equal(sum(abs(G_filtered_3 - G_filtered_2_split_G) < 1e-12), prod(dim(G_filtered_3)))
    ### Using additional filtering by sparsity (only works if G has missing data)
    G = fn_load_genotype(list_sim$fname_geno_vcf, min_depth=10, max_depth=500)
    colnames(G) = gsub("allele_1", "A", colnames(G)) ### Rename allele_1 and allele_alt to A and T, respectively to allow filtering using a SNP list
    colnames(G) = gsub("allele_alt", "T", colnames(G)) ### Rename allele_1 and allele_alt to A and T, respectively to allow filtering using a SNP list
    G_filtered_4 = fn_filter_genotype(G=G, maf=0.05, max_n_alleles=1, fname_snp_list=fname_snp_list,
        max_sparsity_per_locus=0.5,
        frac_topmost_sparse_loci_to_remove=0.01,
        n_topmost_sparse_loci_to_remove=100,
        max_sparsity_per_sample=0.5,
        frac_topmost_sparse_samples_to_remove=0.3,
        n_topmost_sparse_samples_to_remove=10,
        verbose=TRUE)
    expect_equal(nrow(G_filtered_4) < nrow(G), TRUE)
    expect_equal(ncol(G_filtered_4) < ncol(G), TRUE)
    G_filtered_5 = fn_filter_genotype(G=G, maf=0.05,
        frac_topmost_sparse_loci_to_remove=0.1,
        verbose=TRUE)
    G_filtered_6 = fn_filter_genotype(G=G, maf=0.05,
        n_topmost_sparse_loci_to_remove=100,
        verbose=TRUE)
    expect_equal(sum((G_filtered_5 - G_filtered_6) < 1e-7, na.rm=TRUE), sum(!is.na(G_filtered_5)))
    G_filtered_7 = fn_filter_genotype(G=G, maf=0.05,
        frac_topmost_sparse_samples_to_remove=0.1,
        verbose=TRUE)
    G_filtered_8 = fn_filter_genotype(G=G, maf=0.05,
        n_topmost_sparse_samples_to_remove=10,
        verbose=TRUE)
    expect_equal(sum((G_filtered_7 - G_filtered_8) < 1e-7, na.rm=TRUE), sum(!is.na(G_filtered_7)))
    ### Clean-up
    unlink(list_sim$fname_geno_vcf)
    unlink(list_sim$fname_pheno_tsv)
    unlink(fname_snp_list)
})

test_that("fn_save_genotype", {
    set.seed(123)
    list_sim = fn_simulate_data(min_depth=1000, max_depth=1000, verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf, verbose=TRUE)
    fname_Rds = tempfile(fileext=".Rds")
    fname_tsv = tempfile(fileext=".tsv")
    expect_equal(fname_Rds, fn_save_genotype(G=G, fname=fname_Rds, file_type="RDS"))
    expect_equal(fname_tsv, fn_save_genotype(G=G, fname=fname_tsv, file_type="TSV"))
    G_reloaded_Rds = fn_load_genotype(fname_geno=fname_Rds, verbose=TRUE)
    G_reloaded_tsv = fn_load_genotype(fname_geno=fname_tsv, verbose=TRUE)
    expect_equal(G, G_reloaded_Rds)
    expect_equal(G, G_reloaded_tsv)
    unlink(list_sim$fname_geno_vcf)
    unlink(list_sim$fname_pheno_tsv)
    unlink(fname_Rds)
    unlink(fname_tsv)
})

test_that("fn_load_phenotype", {
    set.seed(123)
    list_sim = fn_simulate_data(verbose=TRUE)
    df_y = utils::read.table(list_sim$fname_pheno_tsv, header=TRUE)
    fname_csv = gsub(".tsv$", ".csv", list_sim$fname_pheno_tsv)
    fname_ssv_no_header = gsub(".tsv$", "-no_header.ssv", list_sim$fname_pheno_tsv)
    utils::write.table(df_y, file=fname_csv, sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)
    utils::write.table(df_y, file=fname_ssv_no_header, sep=";", row.names=FALSE, col.names=FALSE, quote=FALSE)
    list_pheno_tsv = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv, verbose=TRUE)
    list_pheno_csv = fn_load_phenotype(fname_pheno=fname_csv, sep=",", verbose=TRUE)
    list_pheno_ssv_no_header = fn_load_phenotype(fname_pheno=fname_ssv_no_header, sep=";", header=FALSE, verbose=TRUE)
    expect_equal(sum(list_pheno_tsv$y - list_pheno_csv$y), 0)
    expect_equal(sum(list_pheno_tsv$y - list_pheno_ssv_no_header$y), 0)
    unlink(list_sim$fname_geno_vcf)
    unlink(list_sim$fname_pheno_tsv)
    unlink(fname_csv)
    unlink(fname_ssv_no_header)
})

test_that("fn_filter_phenotype", {
    set.seed(123)
    list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    list_pheno$y[1] = Inf
    list_pheno$y[2] = NA
    list_pheno_filtered = fn_filter_phenotype(list_pheno, remove_NA=TRUE, verbose=TRUE)
    expect_equal(sum(!is.na(list_pheno_filtered$y)) <= 98, TRUE)
    unlink(list_sim$fname_geno_vcf)
    unlink(list_sim$fname_pheno_tsv)
})

test_that("fn_merge_genotype_genotype", {
    set.seed(123)
    list_sim = gp::fn_simulate_data(verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    ### Union row-wise
    G1 = G[1:ceiling(0.5*nrow(G)), ]
    G2 = G[(ceiling(0.5*nrow(G))+1):nrow(G), ]
    G_merged = fn_merge_genotype_genotype(G1=G1, G2=G2)
    expect_equal(dim(G_merged), dim(G))
    expect_equal(sum(is.na(G_merged)), 0)
    ### Union column-wise
    G1 = G[, 1:ceiling(0.5*ncol(G))]
    G2 = G[, (ceiling(0.5*ncol(G))+1):ncol(G)]
    G_merged = fn_merge_genotype_genotype(G1=G1, G2=G2)
    expect_equal(dim(G_merged), dim(G))
    expect_equal(sum(is.na(G_merged)), 0)
    ### Union row and column-wise with 50% skipped
    G1 = G[seq(from=1, to=nrow(G), by=2), seq(from=1, to=ncol(G), by=2)]
    G2 = G[seq(from=2, to=nrow(G), by=2), seq(from=2, to=ncol(G), by=2)]
    G_merged = fn_merge_genotype_genotype(G1=G1, G2=G2)
    expect_equal(dim(G_merged), dim(G))
    expect_equal(sum(is.na(G_merged)), prod(dim(G))/2)
    ### Intersection row-wise where G2 data is divided by 2
    idx_0 = 1
    idx_1 = ceiling(0.5*nrow(G))
    idx_2 = nrow(G)
    G2 = G[idx_0:idx_1, ] / 2
    G_merged_1 = fn_merge_genotype_genotype(G1=G, G2=G2, str_conflict_resolution="1")
    expect_equal(sum(G_merged_1 - G), 0)
    G_merged_2 = fn_merge_genotype_genotype(G1=G, G2=G2, str_conflict_resolution="2")
    expect_equal(sum(G_merged_2 - rbind(G2[idx_0:idx_1, ], G[(idx_1+1):idx_2, ])), 0)
    G_merged_3 = fn_merge_genotype_genotype(G1=G, G2=G2, str_conflict_resolution="3")
    expect_equal(sum(G_merged_3 - rbind((G[idx_0:idx_1, ] + G2[idx_0:idx_1, ])/2, G[(idx_1+1):idx_2, ])), 0)
    ### Intersection column-wise where G2 data is divided by 2
    idx_0 = 1
    idx_1 = ceiling(0.5*ncol(G))
    idx_2 = ncol(G)
    G2 = G[, idx_0:idx_1] / 2
    G_merged_1 = fn_merge_genotype_genotype(G1=G, G2=G2, str_conflict_resolution="1")
    expect_equal(sum(G_merged_1 - G), 0)
    G_merged_2 = fn_merge_genotype_genotype(G1=G, G2=G2, str_conflict_resolution="2")
    expect_equal(sum(G_merged_2 - cbind(G2[, idx_0:idx_1], G[, (idx_1+1):idx_2])), 0)
    G_merged_3 = fn_merge_genotype_genotype(G1=G, G2=G2, str_conflict_resolution="3")
    expect_equal(sum(G_merged_3 - cbind((G[, idx_0:idx_1] + G2[, idx_0:idx_1])/2, G[, (idx_1+1):idx_2])), 0)
    ### Clean-up
    unlink(list_sim$fname_geno_vcf)
    unlink(list_sim$fname_pheno_tsv)
})

test_that("fn_save_phenotype", {
    set.seed(123)
    list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    fname_tsv = tempfile(fileext=".tsv")
    fname_csv = tempfile(fileext=".csv")
    expect_equal(fname_tsv, fn_save_phenotype(list_pheno=list_pheno, fname=fname_tsv, sep="\t"))
    expect_equal(fname_csv, fn_save_phenotype(list_pheno=list_pheno, fname=fname_csv, sep=","))
    list_pheno_reloaded_tsv = fn_load_phenotype(fname_pheno=fname_tsv)
    list_pheno_reloaded_csv = fn_load_phenotype(fname_pheno=fname_csv, sep=",")
    expect_equal(list_pheno, list_pheno_reloaded_tsv)
    expect_equal(list_pheno, list_pheno_reloaded_csv)

    unlink(list_sim$fname_geno_vcf)
    unlink(list_sim$fname_pheno_tsv)
    unlink(fname_tsv)
    unlink(fname_csv)
})

test_that("fn_merge_genotype_and_phenotype", {
    set.seed(123)
    list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    rownames(G)[1] = "entry_exclude_me"
    rownames(G)[2] = "entry_exclude_me_too"
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    COVAR = matrix(stats::rnorm(n=(10*nrow(G))), nrow=nrow(G))
    rownames(COVAR) = rownames(G); colnames(COVAR) = paste0("covariate_", 1:ncol(COVAR))
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    expect_equal(nrow(list_merged$G), length(list_merged$list_pheno$y))
    expect_equal(nrow(list_merged$G), nrow(list_merged$COVAR))
    expect_equal(sum(is.na(list_merged$list_pheno$y)), 2)
    unlink(list_sim$fname_geno_vcf)
    unlink(list_sim$fname_pheno_tsv)
})

test_that("fn_subset_merged_genotype_and_phenotype", {
    set.seed(123)
    list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    COVAR = matrix(stats::rnorm(n=(10*nrow(G))), nrow=nrow(G))
    rownames(COVAR) = rownames(G); colnames(COVAR) = paste0("covariate_", 1:ncol(COVAR))
    list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    vec_idx = which(list_merged$list_pheno$pop == list_merged$list_pheno$pop[1])
    list_merged_subset = fn_subset_merged_genotype_and_phenotype(list_merged=list_merged, vec_idx=vec_idx, verbose=TRUE)
    expect_equal(nrow(list_merged_subset$G), length(vec_idx))
    expect_equal(length(list_merged_subset$list_pheno$y), length(vec_idx))
    expect_equal(nrow(list_merged_subset$COVAR), length(vec_idx))
    expect_equal(unique(list_merged_subset$list_pheno$pop), list_merged$list_pheno$pop[1])
    unlink(list_sim$fname_geno_vcf)
    unlink(list_sim$fname_pheno_tsv)
})

test_that("fn_estimate_memory_footprint", {
    set.seed(123)
    X = matrix(0.0, nrow=500, ncol=500e3)
    list_mem = fn_estimate_memory_footprint(X=X, verbose=TRUE)
    expect_equal(list_mem$size_X, object.size(X))
    expect_equal(list_mem$size_total > list_mem$size_X, TRUE)
    expect_equal(list_mem$n_threads, 5)
})

test_that(
    "fn_extract_allele_frequencies", {
        print("fn_extract_allele_frequencies:")
        vcf = vcfR::read.vcfR(paste0(dirname_functions, "/../tests/test.vcf"), verbose=FALSE)
        F = fn_extract_allele_frequencies(vcf)
        expect_equal(sum(F - matrix(c(0.0, 0.25, 0.7,
                                        0.1, 0.80, 0.1,
                                        0.2, 0.50, 0.9,
                                        0.0, 1.00, 0.5,
                                        0.3, 0.30, 0.9), byrow=TRUE, nrow=5, ncol=3)), 0)
    }
)

test_that(
    "fn_classify_allele_frequencies", {
        print("fn_classify_allele_frequencies:")
        mat_genotypes = matrix(c(0.0, 0.3, 0.7,
                                    0.1, 0.8, 0.1,
                                    0.2, 0.5, 0.9,
                                    0.0, 1.0, 0.5,
                                    0.3, 0.3, 0.9), byrow=TRUE, nrow=5, ncol=3)
        mx1 = fn_classify_allele_frequencies(mat_genotypes, ploidy=2, strict_boundaries=TRUE)
        mx2 = fn_classify_allele_frequencies(mat_genotypes, ploidy=2, strict_boundaries=FALSE)
        mx3 = fn_classify_allele_frequencies(mat_genotypes, ploidy=3, strict_boundaries=TRUE)
        mx4 = fn_classify_allele_frequencies(mat_genotypes, ploidy=3, strict_boundaries=FALSE)
        expect_equal(mx1, matrix(c(0.0, 0.5, 0.5,
                                    0.5, 0.5, 0.5,
                                    0.5, 0.5, 0.5,
                                    0.0, 1.0, 0.5,
                                    0.5, 0.5, 0.5), byrow=TRUE, nrow=5, ncol=3))
        expect_equal(mx2, matrix(c(0.0, 0.5, 0.5,
                                    0.0, 1.0, 0.0,
                                    0.0, 0.5, 1.0,
                                    0.0, 1.0, 0.5,
                                    0.5, 0.5, 1.0), byrow=TRUE, nrow=5, ncol=3))
        expect_equal(mx3, matrix(c(0.0, 1/3, 2/3,
                                    1/3, 2/3, 1/3,
                                    1/3, 1/3, 2/3,
                                    0.0, 1.0, 1/3,
                                    1/3, 1/3, 2/3), byrow=TRUE, nrow=5, ncol=3))
        expect_equal(mx4, matrix(c(0.0, 1/3, 2/3,
                                    0.0, 2/3, 0.0,
                                    1/3, 2/3, 1.0,
                                    0.0, 1.0, 2/3,
                                    1/3, 1/3, 1.0), byrow=TRUE, nrow=5, ncol=3))
    }
)

test_that(
    "fn_load_genotype", {
        print("fn_load_genotype:")
        G_vcf = fn_load_genotype(fname_rds_or_vcf=paste0(dirname_functions, "/../tests/test.vcf"))
        G_rds = fn_load_genotype(fname_rds_or_vcf=paste0(dirname_functions, "/../tests/test.rds"))
        expect_equal(sum(G_vcf - t(matrix(c(0.0, 0.25, 0.7,
                                            0.1, 0.80, 0.1,
                                            0.2, 0.50, 0.9,
                                            0.0, 1.00, 0.5,
                                            0.3, 0.30, 0.9), byrow=TRUE, nrow=5, ncol=3))), 0)
        expect_equal(G_vcf, G_rds)
    }
)

test_that(
    "fn_filter_loci", {
        print("fn_filter_loci:")
        G = fn_load_genotype(fname_rds_or_vcf=paste0(dirname_functions, "/../tests/test.vcf"))
        G_filtered = fn_filter_loci(G, maf=0.001, sdev_min=0.001)
        expect_equal(G, G_filtered)
    }
)

test_that(
    "fn_load_phenotype", {
        print("fn_load_phenotype:")
        list_y_pop = fn_load_phenotype(fname_csv_txt=paste0(dirname_functions, "/../tests/test_pheno.csv"), header=TRUE, idx_col_id=1, idx_col_pop=2, idx_col_y=3)
        y = list_y_pop$y
        yex = c(0.32, 0.67, 0.93); names(yex) = paste0("Entry-", 1:3)
        expect_equal(y, yex)
        expect_equal(list_y_pop$pop, c("pop-A", "pop-A", "pop-A"))
        expect_equal(list_y_pop$trait_name, "yield")
    }
)

test_that(
    "fn_filter_outlying_phenotypes", {
        print("fn_filter_outlying_phenotypes:")
        list_y_pop = fn_load_phenotype(fname_csv_txt=paste0(dirname_functions, "/../tests/test_pheno.csv"), header=TRUE, idx_col_id=1, idx_col_pop=2, idx_col_y=3)
        list_y_pop_filtered = fn_filter_outlying_phenotypes(list_y_pop)
        expect_equal(sum(list_y_pop$y - c(0.32, 0.67, 0.93)), 0)
        expect_equal(list_y_pop, list_y_pop_filtered)
    }
)

test_that(
    "fn_merge_genotype_and_phenotype", {
        print("fn_merge_genotype_and_phenotype:")
        G = fn_load_genotype(fname_rds_or_vcf=paste0(dirname_functions, "/../tests/test.vcf"))
        list_y_pop = fn_load_phenotype(fname_csv_txt=paste0(dirname_functions, "/../tests/test_pheno.csv"), header=TRUE, idx_col_id=1, idx_col_pop=2, idx_col_y=3)
        COVAR = readRDS(paste0(dirname_functions, "/../tests/test_covar.rds"))
        merged1 = fn_merge_genotype_and_phenotype(G, list_y_pop, COVAR=NULL)
        merged2 = fn_merge_genotype_and_phenotype(G, list_y_pop, COVAR=COVAR)
        out1 = list(G=G, list_y_pop=list_y_pop, COVAR=NULL)
        out2 = list(G=G, list_y_pop=list_y_pop, COVAR=COVAR)
        expect_equal(merged1, out1)
        expect_equal(merged2, out2)
    }
)

suppressWarnings(suppressPackageStartupMessages(library(vcfR)))
suppressWarnings(suppressPackageStartupMessages(library(txtplot)))
suppressWarnings(suppressPackageStartupMessages(library(testthat)))

dirname_functions = dirname(sys.frame(1)$ofile)
# dirname_functions = "/group/pasture/Jeff/genomic_selection/src"; source(file.path(dirname_functions, "load.R"))

### Convert the vcf data into allele frequencies (Output: p x n matrix)
### where we have p biallelic loci and n entries, and allele frequencies refer to the frequency of the reference allele
fn_extract_allele_frequencies = function(vcf) {
    vec_loci_names = paste(vcfR::getCHROM(vcf), vcfR::getPOS(vcf), vcfR::getREF(vcf), sep="_")
    vec_pool_names = colnames(vcf@gt)[-1]
    vec_elements = unique(vcf@gt[, 1])
    if (length(vec_elements) > 1) {
        print("Please check the format of your input vcf file.")
        print(paste0("     - The same fields across loci is required."))
        print(paste0("     - Please pick one field architecture from the following: ", paste(paste0("'", vec_elements, "'"), collapse=","), "."))
        print(paste0("     - Reformat your vcf file to have the same fields across loci."))
        quit()
    }
    if (!grepl("AD", vec_elements) & !grepl("GT", vec_elements)) {
        print("Please check the format of your input vcf file.")
        print(paste0("     - Make sure the 'AD' and/or 'GT' fields are present."))
        print(paste0("     - These are the fields present in your vcf file: ", gsub(":", ", ", vec_elements), "."))
        print(paste0("     - Regenerate your vcf file to include the 'AD' field and/or 'GT' field."))
        quit()
    }
    if (grepl("AD", vec_elements)) {
        mat_allele_counts = vcfR::extract.gt(vcf, element="AD")
        mat_ref_counts = vcfR::masplit(mat_allele_counts, delim=',', record=1, sort=0)
        mat_alt_counts = vcfR::masplit(mat_allele_counts, delim=',', record=2, sort=0)
        ### Set missing allele counts to 0, if the other allele is non-missing and non-zero
        idx_for_ref = which(!is.na(mat_alt_counts) & (mat_alt_counts != 0) & is.na(mat_ref_counts))
        idx_for_alt = which(!is.na(mat_ref_counts) & (mat_ref_counts != 0) & is.na(mat_alt_counts))
        mat_ref_counts[idx_for_ref] = 0
        mat_alt_counts[idx_for_alt] = 0
        ### Calculate reference allele frequencies
        mat_genotypes = mat_ref_counts / (mat_ref_counts + mat_alt_counts)
    } else if (grepl("GT", vec_elements)) {
        GT = vcfR::extract.gt(vcf, element="GT")
        mat_genotypes = matrix(NA, nrow=length(vec_loci_names), ncol=length(vec_pool_names))
        mat_genotypes[(GT == "0/0") | (GT == "0|0")] = 1.0
        mat_genotypes[(GT == "1/1") | (GT == "1|1")] = 0.0
        mat_genotypes[(GT == "0/1") | (GT == "0|1") | (GT == "1|0")] = 0.5        
    } else {
        print("You need to have the 'AD' and/or 'GT' fields present in your vcf file!")
        quit()
    }
    ### Label the loci and pools
    rownames(mat_genotypes) = vec_loci_names
    colnames(mat_genotypes) = vec_pool_names
    ### Output
    return(mat_genotypes)
}

### Genotype classification function assuming biallelic loci (Output: p x n matrix)
### This uses the ploidy level of the species to define these genotype classes, 
### e.g. for diploids we expect 3 genotype classes - AA, AB/BA, and BB, while for tetraploids we expect 5 genotype classes - AAAA, AAAB, AABB, ABBB, and BBBB.
### The default behaviour is to define strict boundaries for the extreme genotypes, i.e. we only consider genotypes to be homozygotes if the allele depth is fixed for one allele.
### Non-strict boundaries is reserved for imputation where we use weighted means and hences the probability of the imputed genotype to belong to one class or another is not strictly bounded at the extremes.
fn_classify_allele_frequencies = function(mat_genotypes, ploidy, strict_boundaries=TRUE) {
    if (ploidy < 1) {
        print(paste0("Are you sure the ploidy is ", ploidy, "X?"))
        print(paste0("  - How on this beautiful universe does that work?"))
        print(paste0("  - Please pick a positive integer!"))
        quit()
    } else if (ploidy != round(ploidy)) {
        print(paste0("Are you sure the ploidy is ", ploidy, "X?"))
        print(paste0("  - How on this beautiful universe does that work?"))
        print(paste0("  - Please pick a positive integer!"))
        quit()
    } else if (ploidy > 1440) {
        print(paste0("Are you sure the ploidy is ", ploidy, "X?"))
        print(paste0("  - What on this beautiful universe are you working on?"))
        print(paste0("  - What species has that many chromosomes?"))
        print(paste0("  - It has more chromosomes that Adder's tongue fern (Ophioglossum reticulatum), last I checked."))
        print(paste0("  - Revise this upper limit if biology begs to differ."))
        quit()
    }
    if (strict_boundaries) {
        ### We are only setting the unfixed frequencies, i.e., (q!=0.0) & (q!=1.0)
        vec_expected_frequencies = c(0, c(0:ploidy))/ploidy
        for (i in 3:(ploidy+1)) {
            # i = 3
            dfreq0 = (vec_expected_frequencies[i-1] - vec_expected_frequencies[i-2])/2
            dfreq1 = (vec_expected_frequencies[i-0] - vec_expected_frequencies[i-1])/2
            idx = (!is.na(mat_genotypes) & 
                (mat_genotypes >  (vec_expected_frequencies[i-1] + dfreq0)) & 
                (mat_genotypes <= (vec_expected_frequencies[i-0] + dfreq1)))
            if (i == ploidy+1) {
                idx = (!is.na(mat_genotypes) & 
                    (mat_genotypes > (vec_expected_frequencies[i-1] + dfreq0)) & 
                    (mat_genotypes < 1.00))
            }
            mat_genotypes[idx] = vec_expected_frequencies[i]
        }
    } else {
        mat_genotypes = round(mat_genotypes * ploidy) / ploidy
    }
    return(mat_genotypes)
}

### Outputs a $n$ entries x $l$ loci x (n_alleles-1) matrix
fn_load_genotype = function(fname_rds_or_vcf, retain_minus_one_alleles_per_locus=TRUE) {
    # fname_rds_or_vcf = "tests/test.rds"
    # fname_rds_or_vcf = "tests/test.vcf"
    # retain_minus_one_alleles_per_locus = TRUE
    ### Load the genotype matrix (n x p)
    G = tryCatch(readRDS(fname_rds_or_vcf), error=function(e){
        vcf = vcfR::read.vcfR(fname_rds_or_vcf, verbose=TRUE)
        return(t(fn_extract_allele_frequencies(vcf)))
    })
    ### If the input genotype matrix is non-numeric, then we assume biallelic loci, e.g. "AA", "AB", and "BB".
    ### Convert them into numerics setting the first allele as 0 and the second as 1, such that "AA" = 0.0, "AB" = 1.0, and "BB" = 2.0
    ### Assuming "A" is the reference allele, then what we end up with is the alternative allele counts
    if (is.numeric(G)==FALSE) {
        print("The input genotype data is non-numeric. We are converting them into a numeric matrix assuming biallelic loci.")
        G_numeric = matrix(0, nrow=nrow(G), ncol=ncol(G))
        rownames(G_numeric) = rownames(G)
        colnames(G_numeric) = colnames(G)
        pb = txtProgressBar(min=0, max=ncol(G), style=3)
        for (j in 1:ncol(G)) {
            # j = 1
            alleles = sort(unique(unlist(strsplit(unique(G[, j]), ""))))
            if (length(alleles) == 2) {
                G_numeric[, j] = unlist(lapply(strsplit(gsub(alleles[2], 1, gsub(alleles[1], 0, G[, j])), ""), FUN=function(x){as.numeric(x[1]) + as.numeric(x[2])}))
            }
            if (length(alleles) > 2) {
                print(paste0("Error reading genotype file: ", fname_rds_or_vcf, "."))
                print(paste0("      - Are you certain that your file is biallelic diploid, i.e. it has a maximum of 2 alleles per locus?"))
                print(paste0("      - Then explain these alleles: ", paste(paste0("'", alleles, "'"), collapse=","), "."))
                print(paste0("      - Do these look like your genotype data: ", paste(head(G_numeric[, j]), collapse=","), "?"))
                print(paste0("      - These too: ", paste(tail(G_numeric[, j]), collapse=","), "?"))
                quit()
            }
            setTxtProgressBar(pb, j)
        }
        close(pb)
        G = G_numeric
    }
    ### Retain a-1 allele/s per locus, i.e. remove duplicates assuming all loci are biallelic
    if (retain_minus_one_alleles_per_locus==TRUE) {
        n = nrow(G)
        p = ncol(G)
        pool_names = rownames(G)
        alleles_across_loci_names = colnames(G)
        first_locus_for_assertions = unlist(strsplit(alleles_across_loci_names[1], "_"))
        if (length(first_locus_for_assertions) < 3) {
            print("Assuming biallelic loci without allele names and one of the alleles are already omitted across all loci. We are now omitting loci or columns with duplicate names.")
            vec_idx_retained = which(!duplicated(colnames(G)))
        } else {
            if (is.na(suppressWarnings(as.numeric(first_locus_for_assertions[length(first_locus_for_assertions)-1])))) {
                print("Assuming biallelic loci without allele names and one of the alleles are already omitted across all loci. We are now omitting loci or columns with duplicate names.")
                vec_idx_retained = which(!duplicated(colnames(G)))
            } else {
                loci_names_tmp = unlist(lapply(alleles_across_loci_names, FUN=function(x){y=unlist(strsplit(x, "_")); l=length(y); chr=paste(y[1:(l-2)], collapse="_"); pos=y[(l-1)]; return(paste0(chr, "_", pos))}))
                loci_names = sort(unique(loci_names_tmp))
                vec_idx_retained = which(!duplicated(loci_names_tmp))
            }
        }
        G = G[, vec_idx_retained]
    }
    return(G)
}

### Filter genotypes by minimum allele frequency and minimum variation within locus (in an attempt to prevent duplicating the intercept)
fn_filter_loci = function(G, maf=0.001, sdev_min=0.001) {
    # fname_rds_or_vcf = "/group/pasture/Jeff/genomic_selection/tests/test.vcf"
    # G = fn_load_genotype(fname_rds_or_vcf)
    # maf = 0.001
    # sdev_min = 0.001
    freqs = colMeans(G)
    sdev = apply(G, MARGIN=2, FUN=sd)
    idx = which((freqs >= maf) & (freqs <= (1-maf)) & (sdev >= sdev_min))
    G = G[, idx]
    gc()
    return(G)
}

### Load phenotype data from a comma-delimited file
### Outputs a named vector with $n$ entries elements
fn_load_phenotype = function(fname_csv_txt, sep=",", header=TRUE, idx_col_id=1, idx_col_pop=2, idx_col_y=3, na.strings=c("", "-", "NA", "na", "NaN", "missing", "MISSING")) {
    # fname_csv_txt = "/group/pasture/Jeff/genomic_selection/tests/test_pheno.csv"; sep=","; header=TRUE; idx_col_id=1; idx_col_pop=2; idx_col_y=3; na.strings=c("", "-", "NA", "na", "NaN", "missing", "MISSING")
    # fname_csv_txt = "/group/pasture/Jeff/genomic_selection/tests/grape_pheno.txt"; sep="\t"; header=TRUE; idx_col_id=1; idx_col_pop=2; idx_col_y=4; na.strings=c("", "-", "NA", "na", "NaN", "missing", "MISSING")
    df = read.table(fname_csv_txt, sep=sep, header=header, na.strings=na.strings)
    entry = as.character(df[, idx_col_id])
    pop = as.character(df[, idx_col_pop])
    y = df[, idx_col_y]
    if (is.numeric(y)==FALSE) {
        print(paste0("Phenotype file: ", fname_csv_txt, ", contains non-numeric data at column ", idx_col_y, "."))
        print(paste0("      - Are you certain that your file is separated by: '", sep, "'?"))
        print(paste0("      - Is the phenotype data column really at column '", idx_col_y, "'?"))
        print(paste0("      - Are your missing data encoded as any of these: ", paste(paste0("'", na.strings, "'"), collapse=", "), "?"))
        print(paste0("      - Do these look like numbers to you: ", paste(head(y), collapse=", "), "?"))
        print(paste0("      - These too: ", paste(tail(y), collapse=", "), "?"))
        quit()
    }
    names(y) = entry
    if (header==TRUE) {
        trait_name = colnames(df)[idx_col_y]
    } else {
        trait_name = paste0("trait_", idx_col_y)
    }
    return(list(y=y, pop=pop, trait_name=trait_name))
}

fn_filter_outlying_phenotypes = function(list_y_pop, verbose=FALSE) {
    # fname_csv_txt = "/group/pasture/Jeff/genomic_selection/tests/test_pheno.csv"; header=TRUE; idx_col_id=1; idx_col_pop=2; idx_col_y=3
    # y = fn_load_phenotype(fname_csv_txt)$y
    # sigma = 1
    y = list_y_pop$y
    pop = list_y_pop$pop
    trait_name = list_y_pop$trait_name
    ### Identify outliers with boxplot, i.e. values beyond -2.698 standard deviations (definition of R::boxplot whiskers)
    b = boxplot(y, plot=FALSE)
    idx = which(!(y %in% b$out))
    if (verbose) {
        print("Before removing outlier/s:")
        print(paste0("n=", length(y)))
        txtplot::txtdensity(y[!is.na(y)])
        print("After removing outlier/s:")
        print(paste0("n=", length(y[idx])))
        txtplot::txtdensity(y[idx][!is.na(y[idx])])
    }
    return(list(y=y[idx], pop=pop[idx], trait_name=trait_name))
}

### Merge genotype and phenotype data by their entry names, i.e. rownames for G and names for y
fn_merge_genotype_and_phenotype = function(G, list_y_pop, COVAR=NULL, verbose=FALSE) {
    # G = fn_load_genotype(fname_rds_or_vcf=paste0("/group/pasture/Jeff/genomic_selection/tests/test.vcf"))
    # list_y_pop = fn_load_phenotype(fname_csv_txt=paste0("/group/pasture/Jeff/genomic_selection/tests/test_pheno.csv"), header=TRUE, idx_col_id=1, idx_col_pop=2, idx_col_y=3)
    # COVAR=NULL; verbose=FALSE
    ### All samples with genotype data will be included and samples without phenotype data will be set to NA (all.x=TRUE)
    M = merge(data.frame(id=rownames(G), G), data.frame(id=names(list_y_pop$y), pop=list_y_pop$pop, y=list_y_pop$y), by="id", all.x=TRUE)
    idx_G = which(!(colnames(M) %in% c("id", "pop", "y")))
    if (!is.null(COVAR)) {
        vec_colnames = colnames(COVAR)
        M = merge(M, data.frame(id=rownames(COVAR), COVAR), by="id", all.x=TRUE)
        idx_G = which(!(colnames(M) %in% c("id", "pop", "y", gsub("-", ".", vec_colnames)))) ### merging converts all the dashes into dots in the column names
        idx_C = which(colnames(M) %in% gsub("-", ".", vec_colnames)) ### merging converts all the dashes into dots in the column names
        COVAR = as.matrix(M[, idx_C])
        rownames(COVAR) = M$id
        colnames(COVAR) = vec_colnames ### Revert to original column names
    }
    if (sum(!is.na(M$y)) == 0) {
        print("Genotype and phenotype data do not match, i.e. none of the entries are common across these datasets.")
        print("IDs in the phenotype data:")
        print(paste(c(head(names(list_y_pop$y)), "..."), collapse=", "))
        print("IDs in the genotype data:")
        print(paste(c(head(rownames(G)), "..."), collapse=", "))
        return(NULL)
    }
    G = as.matrix(M[, idx_G])
    rownames(G) = M$id
    y = M$y
    names(y) = M$id
    pop = M$pop
    trait_name = list_y_pop$trait_name
    if (verbose) {
        print("Genotype distribution:")
        print(paste0("p=", ncol(G)))
        txtplot::txtdensity(G)
        print("Phenotype distribution:")
        print(paste0("n=", length(y)))
        txtplot::txtdensity(y[!is.na(y)])
        if (is.null(COVAR)==FALSE) {
            print("Covariate distribution:")
            print(paste0("m=", ncol(COVAR)))
            txtplot::txtdensity(COVAR)
        } else {
            print("Covariate is null.")
        }
    }
    return(list(G=G, list_y_pop=list(y=y, pop=pop, trait_name=trait_name), COVAR=COVAR))
}

#####################################################################
############################### Tests ###############################
#####################################################################
tests_load = function() {
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
}

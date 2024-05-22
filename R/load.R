suppressWarnings(suppressPackageStartupMessages(library(simquantgen)))
suppressWarnings(suppressPackageStartupMessages(library(vcfR)))
suppressWarnings(suppressPackageStartupMessages(library(txtplot)))
suppressWarnings(suppressPackageStartupMessages(library(testthat)))

# dirname_functions = dirname(sys.frame(1)$ofile)
# dirname_functions = "/group/pasture/Jeff/genomic_selection/src"; source(file.path(dirname_functions, "load.R"))

### Create an error class
setClass("gpError", representation(code="numeric", message="character"), prototype(code=0, message="Empty error message. Please define."))
### Create an error chaining method
setGeneric("chain", function(x, y){
    standardGeneric("chain")
})
setMethod(f="chain",
  signature=c(x="gpError", y="gpError"),
  function(x, y) {
    codes = c(x@code, y@code)
    messages = c(x@message, y@message)
    return(new("gpError", code=codes, message=messages))
  }
)
# err1 = new("gpError", code=1, message="message 1")
# err2 = new("gpError", code=2, message="message 2")
# err3 = chain(err1, err2)
# err4 = chain(err1, err3)



#' Simple wrapper of simquantgen simulation of 10 standard normally distributed QTL additive effects on 5-chromosome genome and a single trait at 50% heritability
#' 
#' @param n number of samples [Default=100]
#' @param l number of loci [Default=1000]
#' @param ploidy ploidy level which can refer to the number of haploid genomes to simulate pools [Default=42]
#' @param n_alleles macimum number of alleles per locus [Default=2]
#' @param depth maximum depth per locus [Default=100]
#' @param seed randomisation seed for replicability [Default=12345]
#' @param save_geno_vcf save the genotype data as a vcf file [Default=FALSE]
#' @param save_geno_tsv save the genotype data as a tab-delimited allele frequency table file [Default=FALSE]
#' @param save_geno_rds save the named genotype matrix as an Rds file [Default=FALSE]
#' @param save_pheno_tsv save the phenotype data as a tab-delimited file [Default=FALSE]
#' @param verbose show simulation messages? [Default=FALSE]
#' @returns
#' vcf: simulated genotype data as a vcfR object
#' df: simulated phenotype data as a data frame
#' fname_geno_vcf: filename of the simulated genotype data as a vcf file
#' fname_geno_tsv: filename of the simulated genotype data as a tab-delimited allele frequency table file
#' fname_geno_rds: filename of the simulated named genotype matrix as an Rds file
#' fname_pheno_tsv: filename of the simulated phenotype data as a tab-delimited file
#' @examples
#' list_sim = fn_simulate_data(verbose=TRUE, save_geno_vcf=TRUE, save_geno_tsv=TRUE, save_geno_rds=TRUE, save_pheno_tsv=TRUE)
#' @export
fn_simulate_data = function(n=100, l=1000, ploidy=42, n_alleles=2, depth=100, seed=12345, save_geno_vcf=FALSE, save_geno_tsv=FALSE, save_geno_rds=FALSE, save_pheno_tsv=FALSE, verbose=FALSE) {
    ###################################################
    ### TEST
    # n = 100
    # l = 1000
    # ploidy = 42
    # n_alleles = 2
    # depth = 100
    # pheno_reps = 1
    # seed = 12345
    # save_geno_vcf = TRUE
    # save_geno_tsv = TRUE
    # save_geno_rds = TRUE
    # save_pheno_tsv = TRUE
    # verbose = TRUE
    ###################################################
    set.seed(seed)
    G = simquantgen::fn_simulate_genotypes(n=n, l=l, ploidy=ploidy, n_alleles=n_alleles, verbose=verbose)
    list_Y_b_E_b_epi = simquantgen::fn_simulate_phenotypes(G=G, n_alleles=n_alleles, dist_effects="norm", n_effects=10, h2=0.5, pheno_reps=1, verbose=FALSE)
    y = list_Y_b_E_b_epi$Y[,1]
    n = nrow(G)
    p = ncol(G)
    vec_ids = rownames(G)
    vec_loci = colnames(G)
    vec_chr = unlist(lapply(strsplit(vec_loci, "-"), FUN=function(x){x[1]}))
    vec_pos = unlist(lapply(strsplit(vec_loci, "-"), FUN=function(x){x[2]}))
    vec_ref = unlist(lapply(strsplit(vec_loci, "-"), FUN=function(x){x[3]}))
    vec_alt = rep(paste0("allele_", n_alleles), times=p)
    META = c("##fileformat=VCFv", paste0("##", vec_chr), "##Extracted from text file.")
    FIX = cbind(vec_chr, vec_pos, vec_loci, vec_ref, vec_alt, rep(NA, each=p), rep("PASS", each=p), rep(NA, each=p))
    colnames(FIX) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    GT_AD = matrix("", nrow=p, ncol=(n+1))
    GT_AD[,1] = "GT:AD"
    colnames(GT_AD) = c("FORMAT", vec_ids)
    if (verbose) {pb = txtProgressBar(min=0, max=n, style=3)}
    for (i in 1:n) {
        for (j in 1:p) {
            ### We're assuming that the counts in G are for the reference alleles
            g = round(2 * G[i, j])
            if (is.na(g)) {
                GT_AD[j, (i+1)] = "./."
            } else {
                if (g == 0) {
                    GT_AD[j, (i+1)] = "1/1"
                } else if (g == 1) {
                    GT_AD[j, (i+1)] = "0/1"
                } else if (g == 2) {
                    GT_AD[j, (i+1)] = "0/0"
                } else {
                    GT_AD[j, (i+1)] = "./."
                }
            }
            tot_depth = sample(c(ceiling(depth/2):depth), size=1)
            ref_depth = round(tot_depth*G[i, j])
            alt_depth = tot_depth - ref_depth
            GT_AD[j, (i+1)] = paste0(GT_AD[j, (i+1)], ":", ref_depth, ",", alt_depth)
        }
        if (verbose) {setTxtProgressBar(pb, i)}
    }
    if (verbose) {close(pb)}
    ### Load dummy vcf file from pinfsc50 package which comes with vcfR as one of its dependencies
    vcf_dummy = vcfR::read.vcfR(system.file("extdata", "pinf_sc50.vcf.gz", package = "pinfsc50"), verbose=FALSE)
    ### Create the new vcfR object
    vcf = vcf_dummy
    vcf@meta = META
    vcf@fix = FIX
    vcf@gt = GT_AD
    if (verbose) {print(vcf)}
    ### Define the phenotype data frame
    df = data.frame(id=names(y), trait=y); rownames(df) = NULL
    if (verbose) {
        print("Simulated allele frequency distribution:")
        print(paste0(c("q_min=", "q_max="), round(range(G), 4)))
        txtplot::txtdensity(G)
        print("Simulated phenotype distribution:")
        print(paste0(c("y_min=", "y_max="), round(range(df$trait), 4)))
        txtplot::txtdensity(df$trait)
    }
    ### Save into files
    date = gsub("-", "", gsub("[.]", "", gsub(":", "", gsub(" ", "", as.character(Sys.time())))))
    fname_geno_vcf = NULL
    fname_geno_tsv = NULL
    fname_geno_rds = NULL
    fname_pheno_tsv = NULL
    if (save_geno_vcf) {
        fname_geno_vcf = file.path(getwd(), paste0("simulated_genotype-", date, ".vcf.gz"))
        vcfR::write.vcf(vcf, file=fname_geno_vcf)
        print("Output genotype file (sync.gz):")
        print(fname_geno_vcf)
    }
    if (save_geno_tsv | save_geno_rds) {
        df_geno = data.frame(chr=vec_chr, pos=vec_pos, allele=vec_ref, t(G))
        fname_geno_tsv = file.path(getwd(), paste0("simulated_genotype-", date, ".tsv"))
        write.table(df_geno, file=fname_geno_tsv, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        print("Output genotype file (tsv):")
        print(fname_geno_tsv)
        if (save_geno_rds) {
            fname_geno_rds = file.path(getwd(), paste0("simulated_genotype-", date, ".Rds"))
            saveRDS(G, file=fname_geno_rds)
            print("Output genotype file (Rds):")
            print(fname_geno_rds)
        }
    }
    if (save_pheno_tsv) {
        fname_pheno_tsv = file.path(getwd(), paste0("simulated_phenotype-", date, ".tsv"))
        write.table(df, file=fname_pheno_tsv, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        print("Output phenotype file (tsv):")
        print(fname_pheno_tsv)
    }
    return(list(
        vcf=vcf,
        df=df,
        fname_geno_vcf=fname_geno_vcf,
        fname_geno_tsv=fname_geno_tsv,
        fname_geno_rds=fname_geno_rds,
        fname_pheno_tsv=fname_pheno_tsv
    ))
}

#' Convert the vcf data into allele frequencies
#' where we have p biallelic loci and n entries, and allele frequencies refer to the frequency of the reference allele
#' 
#' @param vcf simulated genotype data as a vcfR object
#' @param verbose show format conversion messages? [Default=FALSE]
#' @returns
#' named n samples x p loci-alleles matrix
#' @examples
#' list_sim = fn_simulate_data(verbose=TRUE)
#' mat_genotypes = fn_extract_allele_frequencies(vcf=list_sim$vcf, verbose=TRUE)
#' @export
fn_extract_allele_frequencies = function(vcf, verbose=FALSE) {
    ###################################################
    ### TEST
    # vcf = fn_simulate_data(verbose=TRUE, save_geno_vcf=FALSE)$vcf
    # verbose = TRUE
    ###################################################
    ### Check input type
    if (class(vcf) != "vcfR") {
        error = list(code=101, message="Error in load::fn_extract_allele_frequencies(vcf): vcf is not a vcfR object.")
        return(error)
    }
    ### Extract loci and pool/sample names
    vec_loci_names = paste(vcfR::getCHROM(vcf), vcfR::getPOS(vcf), vcfR::getREF(vcf), sep="_")
    vec_pool_names = colnames(vcf@gt)[-1]
    vec_elements = unique(vcf@gt[, 1])
    if (length(vec_elements) > 1) {
        error = list(code=102, 
            message=paste0(
                "Error in load::fn_extract_allele_frequencies(vcf): please check the format of your input vcf file. ",
                "The same fields across loci is required. ",
                "Please pick one field architecture from the following: ", paste(paste0("'", vec_elements, "'"), collapse=","), ". ",
                "Reformat your vcf file to have the same fields across loci."))
        return(error)
    }
    if (!grepl("AD", vec_elements) & !grepl("GT", vec_elements)) {
        error = list(code=103, 
            message=paste0(
                "Error in load::fn_extract_allele_frequencies(vcf): please check the format of your input vcf file. ",
                "Make sure the 'AD' and/or 'GT' fields are present. ",
                "These are the fields present in your vcf file: ", gsub(":", ", ", vec_elements), ". ",
                "Regenerate your vcf file to include the 'AD' field and/or 'GT' field."))
        return(error)
    }
    ### Extract genotype data where the AD field takes precedence over the GT field
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
        error = list(code=104, message="Error in load::fn_extract_allele_frequencies(vcf): vcf needs to have the 'AD' and/or 'GT' fields present.")
        return(error)
    }
    ### Label the loci and pools
    rownames(mat_genotypes) = vec_loci_names
    colnames(mat_genotypes) = vec_pool_names
    ### Allele frequency distribution
    if (verbose) {
        print("Allele frequency distribution:")
        print(paste0(c("q_min=", "q_max="), round(range(mat_genotypes), 4)))
        txtplot::txtdensity(mat_genotypes)
    }
    ### Output n x p matrix of allele frequencies
    return(t(mat_genotypes))
}

#' Genotype classification function assuming biallelic loci (Output: p x n matrix)
#' This uses the ploidy level of the species to define these genotype classes, 
#' e.g. for diploids we expect 3 genotype classes - AA, AB/BA, and BB, while for tetraploids we expect 5 genotype classes - AAAA, AAAB, AABB, ABBB, and BBBB.
#' The default behaviour is to define strict boundaries for the extreme genotypes, i.e. we only consider genotypes to be homozygotes if the allele depth is fixed for one allele.
#' Non-strict boundaries is reserved for imputation where we use weighted means and hences the probability of the imputed genotype to belong to one class or another is not strictly bounded at the extremes.
#' 
#' @param vcf simulated genotype data as a vcfR object
#' @param verbose show format conversion messages? [Default=FALSE]
#' @returns
#' named n samples x p loci-alleles matrix
#' @examples
#' list_sim = fn_simulate_data(ploidy=10, verbose=TRUE)
#' mat_genotypes = fn_extract_allele_frequencies(vcf=list_sim$vcf, verbose=TRUE)
#' mat_genotypes = fn_classify_allele_frequencies(mat_genotypes=mat_genotypes, ploidy=4, strict_boundaries=FALSE, verbose=TRUE)
#' @export
fn_classify_allele_frequencies = function(mat_genotypes, ploidy, strict_boundaries=FALSE, verbose=FALSE) {
    ###################################################
    ### TEST
    # vcf = fn_simulate_data(verbose=TRUE, save_geno_vcf=FALSE)$vcf
    # mat_genotypes = fn_extract_allele_frequencies(vcf=vcf, verbose=TRUE)
    # ploidy = 2
    # strict_boundaries = FALSE
    ###################################################
    if (ploidy < 1) {
        error = list(code=105, message=paste0(
            "Error in load::fn_classify_allele_frequencies(...): Are you sure the ploidy is ", ploidy, "X?",
            "How on this beautiful universe does that work?",
            "Please pick a positive integer!"))
        return(error)
    } else if (ploidy != round(ploidy)) {
        error = list(code=106, message=paste0(
            "Error in load::fn_classify_allele_frequencies(...): Are you sure the ploidy is ", ploidy, "X?",
            "How on this beautiful universe does that work?",
            "Please pick a positive integer!"))
        return(error)
    } else if (ploidy > 1440) {
        error = list(code=107, message=paste0(
            "Error in load::fn_classify_allele_frequencies(...): Are you sure the ploidy is ", ploidy, "X?",
            "What on this beautiful universe are you working on?",
            "What species has that many chromosomes?",
            "It has more chromosomes that Adder's tongue fern (Ophioglossum reticulatum; 2n=1,440), last I checked.",
            "Revise this upper limit if biology begs to differ."))
        return(error)
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
    ### Genotype classses distribution
    if (verbose) {
        print("Genotype classes distribution:")
        print(paste0(c("q_min=", "q_max="), round(range(mat_genotypes), 4)))
        txtplot::txtdensity(mat_genotypes)
        print("Genotype classes:")
        print(table(mat_genotypes))
    }
    ### Output n x p matrix of genotype classes
    return(mat_genotypes)
}

### Outputs a $n$ entries x $l$ loci x (n_alleles-1) matrix
fn_load_genotype = function(fname_genotype, retain_minus_one_alleles_per_locus=TRUE) {
    ###################################################
    ### TEST
    fname_genotype = fn_simulate_data(verbose=TRUE, save_geno_vcf=FALSE)$fname_genotype_vcf
    fname_genotype = fn_simulate_data(verbose=TRUE, save_geno_vcf=FALSE)$fname_genotype_tsv
    fname_genotype = fn_simulate_data(verbose=TRUE, save_geno_vcf=FALSE)$fname_genotype_rds
    retain_minus_one_alleles_per_locus = TRUE
    ###################################################
    ### Load the genotype matrix (n x p)
    G = tryCatch(readRDS(fname_genotype), error=function(e) {
        tryCatch(readRDS(fname_genotype), error=function(e) {
            vcf = vcfR::read.vcfR(fname_genotype, verbose=TRUE)
            mat_genotypes = fn_extract_allele_frequencies(vcf)
            if (class(mat_genotypes) == "gpError") {
                error = chain(mat_genotypes, new("gpError", code=108, message="Error in load::fn_load_genotype: error loading the vcf file."))
                return(error)
            } else {
                return(mat_genotypes)
            }
        })
    })
    ### If the input genotype matrix is non-numeric, then we assume biallelic loci, e.g. "AA", "AB", and "BB".
    ### Convert them into numerics setting the first allele as 0 and the second as 1, such that "AA" = 0.0, "AB" = 1.0, and "BB" = 2.0
    ### Assuming "A" is the reference allele, then what we end up with is the alternative allele counts
    if (is.numeric(G)==FALSE) {
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
                print(paste0("Error reading genotype file: ", fname_genotype, "."))
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

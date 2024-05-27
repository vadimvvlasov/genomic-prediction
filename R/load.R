suppressWarnings(suppressPackageStartupMessages(library(simquantgen)))
suppressWarnings(suppressPackageStartupMessages(library(vcfR)))
suppressWarnings(suppressPackageStartupMessages(library(txtplot)))
suppressWarnings(suppressPackageStartupMessages(library(testthat)))

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

#' Extract names of samples/entries/pools and loci from a numeric (allele frequency) or non-numeric (genotype classes) matrix
#' 
#' @param mat_genotypes n samples x p loci matrix of:
#'  allele frequencies (numeric ranging from 0 to 1) or 
#'  genotype classes (e.g. "AA", "AB", "AAB", "AAAB") 
#'  with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#' @param verbose show names extraction messages? [Default=FALSE]
#' @returns
#' Ok:
#'  $vec_ids: vector of sample/entry/pool names
#'  $vec_chr: vector of chromosome or scaffold names
#'  $vec_pos: vector of numeric positions per chromosome
#'  $vec_all: vector of allele names
#' Err: grError
#' @examples
#' mat_genotypes = simquantgen::fn_simulate_genotypes(verbose=TRUE)
#' list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=mat_genotypes, verbose=TRUE)
#' @export
fn_G_extract_names = function(mat_genotypes, verbose=FALSE) {
    ###################################################
    ### TEST
    # mat_genotypes = simquantgen::fn_simulate_genotypes()
    # verbose = TRUE
    ###################################################
    ### Extract the names of the entries/samples/pools and loci
    n = nrow(mat_genotypes)
    p = ncol(mat_genotypes)
    vec_ids = rownames(mat_genotypes)
    vec_loci = colnames(mat_genotypes)
    n_identifiers = length(unlist(strsplit(vec_loci[1], "\t"))) ### Number of loci identifiers where we expect at least 2
    if (n_identifiers < 2) {
        error = new("gpError",
            code=000,
            message=paste0(
                "Error in load::fn_G_extract_names(...). ",
                "The loci names (column names) are not correctly formatted. ",
                "We expect column names to be tab-delimited ('\t'), where ",
                "the first element refers to the chromosome or scaffold name, ",
                "the second should be numeric which refers to the position in the chromosome/scaffold, and ",
                "subsequent elements are optional which may refer to the allele identifier and other identifiers."))
        return(error)
    }
    mat_loci_ids = matrix(unlist(strsplit(vec_loci, "\t")), byrow=TRUE, ncol=n_identifiers)
    vec_chr = mat_loci_ids[,1]
    vec_pos = as.numeric(mat_loci_ids[,2])
    if (sum(is.na(vec_pos)) > 0) {
        error = new("gpError",
            code=000,
            message=paste0(
                "Error in load::fn_G_extract_names(...). ",
                "The second element of the tab-delimited loci names should be numeric position."))
        return(error)
    }
    if (n_identifiers == 2) {
        if (verbose) {print("The loci are identified by the chromosome and position, hence we are assuming a diploid dataset.")}
        vec_all = rep("allele1", times=p)
    } else {
        if (verbose) {print("The loci are identified by the chromosome, position, allele, and other identifiers which will not be used.")}
        vec_all = mat_loci_ids[,3]
    }
    ### Return list of vectors of chromosome names, positions, and alleles
    return(list(
        vec_ids=vec_ids,
        vec_chr=vec_chr,
        vec_pos=vec_pos,
        vec_all=vec_all
    ))
}

#' Convert a numeric (allele frequency) genotype matrix into a non-numeric (genotype classes) genotype matrix
#' 
#' @param G n samples x p loci matrix of allele frequencies (numeric ranging from 0 to 1) with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#' @param ploidy ploidy level which can refer to the number of haploid genomes to simulate pools [Default=2]
#' @param verbose show numeric to non-numeric genotype data conversion messages? [Default=FALSE]
#' @returns
#' Ok: n samples x p loci matrix of genotype classes
#' Err: grError
#' @examples
#' ploidy = 42
#' G_numeric = simquantgen::fn_simulate_genotypes(ploidy=ploidy, n_alleles=52, verbose=TRUE)
#' G_non_numeric = fn_G_numeric_to_non_numeric(G=G_numeric, ploidy=ploidy, verbose=TRUE)
#' @export
fn_G_numeric_to_non_numeric = function(G, ploidy=2, verbose=FALSE) {
    ###################################################
    ### TEST
    # G = simquantgen::fn_simulate_genotypes(ploidy=42, n_alleles=52, verbose=TRUE)
    # ploidy = 42
    # verbose = TRUE
    ###################################################
    ### Make sure G contains allele frequencies
    if (sum(((G < 0) | (G > 1) | is.infinite(G)), na.rm=TRUE) > 0) {
        error = new("gpError",
            code=000,
            message=paste0(
                "Error in load::fn_G_numeric_to_non_numeric(...). ",
                "We are expecting a matrix allele frequencies but ",
                "we are getting negative values and/or values greater than 1 and/or infinite values."))
        return(error)
    }
    ### Extract sample/entry/pool, and loci names
    list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=G, verbose=verbose)
    ### Extract the unique loci and sort per chromosome
    df_loci = as.data.frame(matrix(unlist(strsplit(unique(paste0(list_ids_chr_pos_all$vec_chr, "\t", list_ids_chr_pos_all$vec_pos)), "\t")), byrow=TRUE, ncol=2))
    colnames(df_loci) = c("chr", "pos")
    df_loci$pos = as.numeric(df_loci$pos)
    idx = order(df_loci$chr, df_loci$pos, decreasing=FALSE)
    df_loci = df_loci[idx, ]
    vec_loci = paste0(df_loci$chr, "\t", df_loci$pos)
    ### Populate the non-numeric genotype matrix
    p = length(vec_loci)
    n = nrow(G)
    G_non_numeric = matrix("", nrow=n, ncol=p)
    rownames(G_non_numeric) = rownames(G)
    colnames(G_non_numeric) = vec_loci
    ### Iterate across sorted loci
    if (verbose) {pb = txtProgressBar(min=0, max=p, style=3)}
    for (j in 1:p) {
        # j = 1
        vec_chr_pos = unlist(strsplit(vec_loci[j], "\t"))
        vec_idx = which(
            (list_ids_chr_pos_all$vec_chr == vec_chr_pos[1]) &
            (list_ids_chr_pos_all$vec_pos == as.numeric(vec_chr_pos[2])))
        if (sum(rowSums(G[, vec_idx, drop=FALSE], na.rm=TRUE) == 1) == n) {
            ### All alleles are represented
            n_alleles = length(vec_idx)
        } else {
            ### The trailing alleles was omitted
            n_alleles = length(vec_idx) + 1
        }
        vec_allele_choices = c(LETTERS, letters)[1:n_alleles]
        for (i in 1:nrow(G)) {
            # i = 1
            str_code = ""
            for (k in 1:n_alleles) {
                # k = 1
                if ((k == n_alleles) & (length(vec_idx) < n_alleles)) {
                    ### Add the omitted trailing alternative allele
                    count = round((1-sum(G[i, vec_idx]))*ploidy)
                    str_code = paste0(str_code, paste(rep(vec_allele_choices[k], times=count), collapse=""))
                } else {
                    count = round(G[i, vec_idx[k]]*ploidy)
                    str_code = paste0(str_code, paste(rep(vec_allele_choices[k], times=count), collapse=""))   
                }
            }
            G_non_numeric[i, j] = str_code
        }
        if (verbose) {setTxtProgressBar(pb, j)}
    }
    if (verbose) {close(pb)}
    return(G_non_numeric)
}

#' Convert a non-numeric genotype classes matrix into a numeric (allele frequency) genotype matrix
#' 
#' @param G n samples x p loci matrix of genotype classes with non-null row and column names.
#'  The genotype classes are represented by uppercase and/or lowercase letters,
#'  where the total number of letters correspond to the ploidy level,
#'  and the number of unique letters refer to the number of alleles per locus.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  and the second should be numeric which refers to the position in the chromosome/scaffold.
#' @param verbose show non-numeric to numeric genotype data conversion messages? [Default=FALSE]
#' @returns
#' Ok: n samples x p loci matrix of genotype classes (numeric ranging from 0 to 1)
#' Err: grError
#' @examples
#' ploidy = 42
#' G_numeric = simquantgen::fn_simulate_genotypes(ploidy=ploidy, n_alleles=52, verbose=TRUE)
#' G_non_numeric = fn_G_numeric_to_non_numeric(G=G_numeric, ploidy=ploidy, verbose=TRUE)
#' G_numeric_back = fn_G_non_numeric_to_numeric(G=G_non_numeric, verbose=TRUE)
#' @export
fn_G_non_numeric_to_numeric = function(G_non_numeric, verbose=FALSE) {
    ###################################################
    ### TEST
    # G = simquantgen::fn_simulate_genotypes(ploidy=42, n_alleles=52, verbose=TRUE)
    # G_non_numeric = fn_G_numeric_to_non_numeric(G=G, ploidy=42, verbose=TRUE)
    # verbose = TRUE
    ###################################################
    ### Extract ploidy where we assume the same ploidy across the entire data set
    ploidy = length(unlist(strsplit(G_non_numeric[1,1], "")))
    for (i in sample(1:nrow(G_non_numeric), size=100, replace=FALSE)) {
        for (j in sample(1:ncol(G_non_numeric), size=100, replace=FALSE)) {
            if (length(unlist(strsplit(G_non_numeric[i, j], ""))) != ploidy) {
                error = new("gpError",
                    code=000,
                    message=paste0(
                        "Error in load::fn_G_non_numeric_to_numeric(...). ",
                        "The ploidy level is not consistent across the data set. ",
                        "From the first locus of the first sample we expected a ploidy level of ", ploidy, "X."))
                return(error)
            }
        }
    }
    ### Extract the loci and alleles per locus
    vec_loci_names = colnames(G_non_numeric)
    list_allele_names = list()
    p = 0
    if (verbose) {pb = txtProgressBar(min=0, max=ncol(G_non_numeric), style=3)}
    for (j in 1:ncol(G_non_numeric)) {
        # j = 1
        vec_alleles = sort(unique(unlist(strsplit(G_non_numeric[, j], ""))))
        eval(parse(text=paste0("list_allele_names$`", vec_loci_names[j], "` = vec_alleles")))
        p = p + length(vec_alleles)
        if (verbose) {setTxtProgressBar(pb, j)}
    }
    if (verbose) {close(pb)}
    ### Instantiate the numeric allele frequency genotype matrix
    n = nrow(G_non_numeric)
    G = matrix(0, nrow=n, ncol=p)
    rownames(G) = rownames(G_non_numeric)
    ### Populate
    vec_colnames = c()
    idx_locus_allle = 0
    if (verbose) {pb = txtProgressBar(min=0, max=ncol(G_non_numeric), style=3)}
    for (j in 1:ncol(G_non_numeric)) {
        # j = 1
        vec_alleles = list_allele_names[[j]]
        list_geno_classes = strsplit(G_non_numeric[, j], "")
        for (k in 1:length(vec_alleles)) {
            # k = 1
            idx_locus_allle = idx_locus_allle + 1
            vec_colnames = c(vec_colnames, paste0(vec_loci_names[j], "\t", vec_alleles[k]))
            for (i in 1:n) {
                # i = 1
                G[i, idx_locus_allle] = sum(list_geno_classes[[i]] == vec_alleles[k]) / ploidy
            }
        }
        if (verbose) {setTxtProgressBar(pb, j)}
    }
    if (verbose) {close(pb)}
    colnames(G) = vec_colnames
    ### Return numeric allele frequency genotype matrix
    return(G)
}
#' Remove trailing alternative allele from a numeric (allele frequency) matrix if it exists
#' 
#' @param G n samples x p loci matrix of allele frequencies (numeric ranging from 0 to 1) with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#' @param verbose show genotype data alternative allele splitting messages? [Default=FALSE]
#' @returns
#' Ok: n samples x p loci matrix of allele frequencies (numeric ranging from 0 to 1)
#' Err: grError
#' @examples
#' G_ref = simquantgen::fn_simulate_genotypes(verbose=TRUE)
#' G_alt = 1 - G_ref; colnames(G_alt) = gsub("allele_1$", "allele_2", colnames(G_alt))
#' G_refalt = cbind(G_ref, G_alt)
#' list_G_G_alt = fn_G_split_off_alternative_allele(G=G_refalt, verbose=TRUE)
#' @export
fn_G_split_off_alternative_allele = function(G, verbose=FALSE) {
    ###################################################
    ### TEST
    # G_ref = simquantgen::fn_simulate_genotypes(verbose=TRUE)
    # G_alt = 1 - G_ref; colnames(G_alt) = gsub("allele_1$", "allele_2", colnames(G_alt))
    # G = cbind(G_ref, G_alt)
    # verbose = TRUE
    ###################################################
    ### Make sure G contains allele frequencies
    if (sum(((G < 0) | (G > 1) | is.infinite(G)), na.rm=TRUE) > 0) {
        error = new("gpError",
            code=000,
            message=paste0(
                "Error in load::fn_G_split_off_alternative_allele(...). ",
                "We are expecting a matrix allele frequencies but ",
                "we are getting negative values and/or values greater than 1 and/or infinite values."))
        return(error)
    }
    ### Extract sample/entry/pool, and loci names
    list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=G, verbose=verbose)
    ### Iterate across loci removing the trailing alternative allele
    vec_loci = sort(unique(paste0(list_ids_chr_pos_all$vec_chr, "\t", list_ids_chr_pos_all$vec_pos)))
    p = length(vec_loci)
    n = nrow(G)
    vec_idx_loci = c()
    vec_idx_alt = c()
    if (ncol(G) > p) {
        if (verbose) {pb = txtProgressBar(min=0, max=p, style=3)}
        for (j in 1:p) {
            # j = 1
            vec_chr_pos = unlist(strsplit(vec_loci[j], "\t"))
            idx = which(
                (list_ids_chr_pos_all$vec_chr == vec_chr_pos[1]) &
                (list_ids_chr_pos_all$vec_pos == as.numeric(vec_chr_pos[2]))
            )
            if (sum(rowSums(G[, idx, drop=FALSE], na.rm=TRUE) == 1) == n) {
                ### Omit the trailing (alternative) allele
                vec_idx_loci = c(vec_idx_loci, idx[1:(length(idx)-1)])
                vec_idx_alt = c(vec_idx_alt, tail(idx, n=1))
            } else {
                ### Otherwise keep all the loci if all frequencies do not sum up to 1
                vec_idx_loci = c(vec_idx_loci, idx)
            }
            if (verbose) {setTxtProgressBar(pb, j)}
        }
        if (verbose) {close(pb)}
        ### Return filtered matrix of allele frequencies
        return(list(
            G=G[, sort(vec_idx_loci), drop=FALSE],
            G_alt=G[, sort(vec_idx_alt), drop=FALSE]))
    } else {
        if (verbose) {"Removal of trailing alternative allele is unnecessary as only the one allele represents each locus."}
        return(list(
            G=G,
            G_alt=NULL))
    }
}

#' Convert numeric allele frequency matrix into a vcfR::vcf object
#' 
#' @param G numeric n samples x p biallelic loci matrix of allele frequencies with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'  Note that this function only accepts biallelic loci.
#' @param min_depth minimum depth per locus [Default=100]
#' @param max_depth maximum depth per locus [Default=1000]
#' @param verbose show allele frequency genotype matrix to vcf conversion messages? [Default=FALSE]
#' @returns
#' Ok: simulated genotype data as a vcfR object
#' Err: grError
#' @examples
#' G = simquantgen::fn_simulate_genotypes(verbose=TRUE)
#' vcf = fn_G_to_vcf(G=G, verbose=TRUE)
#' @export
fn_G_to_vcf = function(G, min_depth=100, max_depth=1000, verbose=FALSE) {
    ###################################################
    ### TEST
    # G = simquantgen::fn_simulate_genotypes()
    # min_depth = 100
    # max_depth = 1000
    # verbose = TRUE
    ###################################################
    ### Make sure G contains allele frequencies
    if (sum(((G < 0) | (G > 1) | is.infinite(G)), na.rm=TRUE) > 0) {
        error = new("gpError",
            code=000,
            message=paste0(
                "Error in load::fn_G_to_vcf(...). ",
                "We are expecting a matrix allele frequencies but ",
                "we are getting negative values and/or values greater than 1 and/or infinite values."))
        return(error)
    }
    ### Split-off the alternative allele from the reference allele
    list_G_G_alt = fn_G_split_off_alternative_allele(G=G, verbose=verbose)
    G_alt = list_G_G_alt$G_alt
    G = list_G_G_alt$G
    n = nrow(G)
    p = ncol(G)
    ### Extract the names of the entries/samples/pools and loci
    list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=G, verbose=verbose)
    ### Make sure the input matrix is biallelic
    ### Also, convert the tabs in the loci-alleles names into dashes so as not to interfere with the VCF format
    vec_loci_names = paste0(list_ids_chr_pos_all$vec_chr, "_", list_ids_chr_pos_all$vec_pos)
    vec_allele_counts_per_locus = table(vec_loci_names)
    if (((sum(vec_allele_counts_per_locus > 2) > 0) & is.null(G_alt)) | ((sum(vec_allele_counts_per_locus >= 2) > 0) & !is.null(G_alt))) {
        error = new("gpError",
            code=000,
            message=paste0(
                "Error in load::fn_G_to_vcf(...). ",
                "Apologies, as this function at the moment can only convert biallelic allele frequency matrices into VCF format."))
        return(error)
    }
    ### Extract the names of the alternative alleles
    if (is.null(G_alt)) {
        vec_alt = rep("allele_alt", times=ncol(G))
    } else {
        ### Convert the tabs in the loci-alleles names into dashes so as not to interfere with the VCF format
        vec_alt = colnames(gsub("\t", "_", list_G_G_alt$G_alt))
    }
    ### Populate required vcf fields with loci identities
    META = c("##fileformat=VCFv", paste0("##", list_ids_chr_pos_all$vec_chr), "##Extracted from text file.")
    FIX = cbind(
        list_ids_chr_pos_all$vec_chr, 
        list_ids_chr_pos_all$vec_pos, 
        vec_loci_names, 
        list_ids_chr_pos_all$vec_all, 
        vec_alt, 
        rep(NA, each=p),
        rep("PASS", each=p),
        rep(NA, each=p))
    colnames(FIX) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    ### Population the GT:AD matrix with strings of genotype class and allele depths
    GT_AD = matrix("", nrow=p, ncol=(n+1))
    GT_AD[,1] = "GT:AD"
    colnames(GT_AD) = c("FORMAT", list_ids_chr_pos_all$vec_ids)
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
            if (min_depth == max_depth) {
                vec_range = c(min_depth, max_depth)
            } else {
                vec_range = c(min_depth:max_depth)
            }
            tot_depth = sample(x=vec_range, size=1) ### Sample the depth from a range from half the input max_depth to the input max_depth
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
    return(vcf)
}

#' Convert biallelic vcf data into allele frequencies
#' 
#' @param vcf biallelic genotype data as a vcfR object
#' @param minimum_depth minimum sequencing depth [Default=5]
#' @param retain_minus_one_alleles_per_locus omit the alternative or trailing allele per locus? [Default=TRUE]
#' @param verbose show vcf to allele frequency genotype matrix conversion messages? [Default=FALSE]
#' @returns
#' Ok: named n samples x p biallelic loci matrix.
#'  Row names can be any string of characters.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#' Err: gpError
#' @examples
#' G = simquantgen::fn_simulate_genotypes(verbose=TRUE)
#' vcf = fn_G_to_vcf(G=G, verbose=TRUE)
#' G_back = fn_vcf_to_G(vcf=vcf, verbose=TRUE)
#' @export
fn_vcf_to_G = function(vcf, minimum_depth=5, retain_minus_one_alleles_per_locus=TRUE, verbose=FALSE) {
    ###################################################
    ### TEST
    # G = simquantgen::fn_simulate_genotypes(verbose=TRUE)
    # vcf = fn_G_to_vcf(G=G, verbose=TRUE)
    # verbose = TRUE
    ###################################################
    ### Check input type
    if (class(vcf) != "vcfR") {
        error = new("gpError", 
            code=000,
            message="Error in load::fn_vcf_to_G(vcf): vcf is not a vcfR object.")
        return(error)
    }
    ### Extract loci and pool/sample names
    vec_loci_ref_names = paste(vcfR::getCHROM(vcf), vcfR::getPOS(vcf), vcfR::getREF(vcf), sep="\t")
    if (!retain_minus_one_alleles_per_locus) {
        vec_loci_alt_names = paste(vcfR::getCHROM(vcf), vcfR::getPOS(vcf), vcfR::getALT(vcf), sep="\t")
    }
    vec_pool_names = colnames(vcf@gt)[-1]
    vec_elements = unique(vcf@gt[, 1])
    if (length(vec_elements) > 1) {
        error = new("gpError", 
            code=000,
            message=paste0(
                "Error in load::fn_vcf_to_G(vcf): please check the format of your input vcf file. ",
                "The same fields across loci is required. ",
                "Please pick one field architecture from the following: ", paste(paste0("'", vec_elements, "'"), collapse=","), ". ",
                "Reformat your vcf file to have the same fields across loci."))
        return(error)
    }
    if (!grepl("AD", vec_elements) & !grepl("GT", vec_elements)) {
        error = new("gpError", 
            code=000,
            message=paste0(
                "Error in load::fn_vcf_to_G(vcf): please check the format of your input vcf file. ",
                "Make sure the 'AD' and/or 'GT' fields are present. ",
                "These are the fields present in your vcf file: ", gsub(":", ", ", vec_elements), ". ",
                "Regenerate your vcf file to include the 'AD' field and/or 'GT' field."))
        return(error)
    }

    ### TODO: handle multi-row multi-allelic loci

    ### TODO: Convert loci into missing if depth is below minimum_depth
    mat_depth = vcfR::extract.gt(vcf, element="DP", as.numeric=TRUE) ### WILL NEED TO ADD THIS DP FIELD ABOVE in fn_G_to_vcf(...)


    ### Extract genotype data where the AD field takes precedence over the GT field
    if (grepl("AD", vec_elements)) {
        mat_allele_counts = vcfR::extract.gt(vcf, element="AD")
        if (length(unlist(strsplit(mat_allele_counts[1,1], ","))) > 2) {
            error = new("gpError",
                code=000,
                message=paste0(
                    "Error in load::fn_vcf_to_G(...). ",
                    "Apologies because at the moment we can only convert biallelic VCF files into an allele frequency matrix."))
            return(error)
        }
        mat_ref_counts = vcfR::masplit(mat_allele_counts, delim=',', record=1, sort=0)
        mat_alt_counts = vcfR::masplit(mat_allele_counts, delim=',', record=2, sort=0)
        ### Set missing allele counts to 0, if the other allele is non-missing and non-zero
        idx_for_ref = which(!is.na(mat_alt_counts) & (mat_alt_counts != 0) & is.na(mat_ref_counts))
        idx_for_alt = which(!is.na(mat_ref_counts) & (mat_ref_counts != 0) & is.na(mat_alt_counts))
        mat_ref_counts[idx_for_ref] = 0
        mat_alt_counts[idx_for_alt] = 0
        ### Calculate reference allele frequencies
        G = mat_ref_counts / (mat_ref_counts + mat_alt_counts)
    } else if (grepl("GT", vec_elements)) {
        GT = vcfR::extract.gt(vcf, element="GT")
        G = matrix(NA, nrow=length(vec_loci_ref_names), ncol=length(vec_pool_names))
        G[(GT == "0/0") | (GT == "0|0")] = 1.0
        G[(GT == "1/1") | (GT == "1|1")] = 0.0
        G[(GT == "0/1") | (GT == "0|1") | (GT == "1|0")] = 0.5
    } else {
        error = new("gpError", 
            code=000,
            message="Error in load::fn_vcf_to_G(vcf): vcf needs to have the 'AD' and/or 'GT' fields present.")
        return(error)
    }
    ### Column-bind reference and alternative allele frequencies
    if (!retain_minus_one_alleles_per_locus) {
        G = rbind(G, 1-G)
    }
    ### Label the loci and pools
    if (!retain_minus_one_alleles_per_locus) {
        rownames(G) = c(vec_loci_ref_names, vec_loci_alt_names)
    } else {
        rownames(G) = vec_loci_ref_names
    }
    colnames(G) = vec_pool_names
    ### Sort by locus, i.e. put the reference and alternative alleles next to each other
    if (!retain_minus_one_alleles_per_locus) {
        vec_idx = order(c(seq(from=1, to=nrow(G), by=2), seq(from=2, to=nrow(G), by=2)), decreasing=FALSE)
        G = G[vec_idx, , drop=FALSE]
    }
    ### Allele frequency distribution
    if (verbose) {
        print("Allele frequency distribution:")
        print(paste0(c("q_min=", "q_max="), round(range(G), 4)))
        txtplot::txtdensity(G)
    }
    ### Output n x p matrix of allele frequencies
    return(t(G))
}

#' Classify or bin allele frequencies into genotype classes.
#' This uses the ploidy level of the species to define these genotype classes, 
#' e.g. for diploids we expect 3 genotype classes - AA, AB/BA, and BB, 
#'while for tetraploids we expect 5 genotype classes - AAAA, AAAB, AABB, ABBB, and BBBB.
#' 
#' @param G numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#' @param ploidy ploidy level which can refer to the number of haploid genomes to simulate pools [Default=2]
#' @param verbose show genotype binning or classification messages? [Default=FALSE]
#' @returns
#' Ok: named n samples x p loci-alleles matrix of numeric genotype classes
#' Err: gpError
#' @examples
#' ploidy = 4
#' G = simquantgen::fn_simulate_genotypes(ploidy=ploidy, verbose=TRUE)
#' G_classes = fn_classify_allele_frequencies(G=G, ploidy=ploidy, verbose=TRUE)
#' @export
fn_classify_allele_frequencies = function(G, ploidy=2, verbose=FALSE) {
    ###################################################
    ### TEST
    # ploidy = 4
    # G = simquantgen::fn_simulate_genotypes(ploidy=ploidy, verbose=TRUE)
    # verbose = TRUE
    ###################################################
    if (ploidy < 1) {
        error = new("gpError", 
            code=000,
            message=paste0(
                "Error in load::fn_classify_allele_frequencies(...): Are you sure the ploidy is ", ploidy, "X?",
                "How on this beautiful universe does that work?",
                "Please pick a positive integer!"))
        return(error)
    } else if (ploidy != round(ploidy)) {
        error = new("gpError", 
            code=000,
            message=paste0(
                "Error in load::fn_classify_allele_frequencies(...): Are you sure the ploidy is ", ploidy, "X?",
                "How on this beautiful universe does that work?",
                "Please pick a positive integer!"))
        return(error)
    } else if (ploidy > 1440) {
        error = new("gpError", 
            code=000,
            message=paste0(
                "Error in load::fn_classify_allele_frequencies(...): Are you sure the ploidy is ", ploidy, "X?",
                "What on this beautiful universe are you working on?",
                "What species has that many chromosomes?",
                "It has more chromosomes that Adder's tongue fern (Ophioglossum reticulatum; 2n=1,440), last I checked.",
                "Revise this upper limit if biology begs to differ."))
        return(error)
    }
    G_classes = round(G * ploidy) / ploidy
    ### Genotype classses distribution
    if (verbose) {
        print("Genotype classes distribution:")
        print(paste0(c("q_min=", "q_max="), round(range(G_classes), 4)))
        txtplot::txtdensity(G_classes)
        print("Genotype classes:")
        print(table(G_classes))
    }
    ### Output n x p matrix of genotype classes
    return(G_classes)
}

#' Simple wrapper of simquantgen simulation of 10 standard normally distributed QTL additive effects 
#' on 5-chromosome genome and a single trait at 50% heritability
#' 
#' @param n number of samples [Default=100]
#' @param l number of loci [Default=1000]
#' @param ploidy ploidy level which can refer to the number of haploid genomes to simulate pools [Default=42]
#' @param n_alleles macimum number of alleles per locus [Default=2]
#' @param seed randomisation seed for replicability [Default=12345]
#' @param save_pheno_tsv save the phenotype data as a tab-delimited file? [Default=TRUE]
#' @param save_geno_vcf save the genotype data as a vcf file? [Default=TRUE]
#' @param save_geno_tsv save the genotype data as a tab-delimited allele frequency table file? [Default=FALSE]
#' @param save_geno_rds save the named genotype matrix as an Rds file? [Default=FALSE]
#' @param non_numeric_Rds save non-numeric Rds genotype file? [Default=FALSE]
#' @param verbose show genotype and phenotype data simulation messages? [Default=FALSE]
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
fn_simulate_data = function(n=100, l=1000, ploidy=42, n_alleles=2, depth=100, seed=12345, save_pheno_tsv=TRUE, save_geno_vcf=TRUE, save_geno_tsv=FALSE, save_geno_rds=FALSE, non_numeric_Rds=FALSE, verbose=FALSE) {
    ###################################################
    ### TEST
    # n = 100
    # l = 1000
    # ploidy = 42
    # n_alleles = 2
    # depth = 100
    # pheno_reps = 1
    # seed = 12345
    # save_pheno_tsv = TRUE
    # save_geno_vcf = TRUE
    # save_geno_tsv = TRUE
    # save_geno_rds = TRUE
    # non_numeric_Rds = TRUE
    # verbose = TRUE
    ###################################################
    ### Set a randomisation seed for repeatability
    set.seed(seed)
    ### Simulate genotype matrix
    G = simquantgen::fn_simulate_genotypes(n=n, l=l, ploidy=ploidy, n_alleles=n_alleles, verbose=verbose)
    ### Simulate phenotype data frame
    list_Y_b_E_b_epi = simquantgen::fn_simulate_phenotypes(G=G, n_alleles=n_alleles, dist_effects="norm", n_effects=10, h2=0.5, pheno_reps=1, verbose=FALSE)
    y = list_Y_b_E_b_epi$Y[,1]
    df = data.frame(id=names(y), trait=y); rownames(df) = NULL
    ### Report distribution of simulation data
    if (verbose) {
        print("Simulated allele frequency distribution:")
        print(paste0(c("q_min=", "q_max="), round(range(G), 4)))
        txtplot::txtdensity(G)
        print("Simulated phenotype distribution:")
        print(paste0(c("y_min=", "y_max="), round(range(df$trait), 4)))
        txtplot::txtdensity(df$trait)
    }
    ### Instantiate output file names
    fname_geno_vcf = NULL
    fname_geno_tsv = NULL
    fname_geno_rds = NULL
    fname_pheno_tsv = NULL
    ### Define date-time-random-number identifier so that we minimise the possibility of unintentional over-writing of the output file/s
    date = gsub("-", "", gsub("[.]", "", gsub(":", "", gsub(" ", "", as.character(Sys.time())))))
    ### Save phenotype file
    if (save_pheno_tsv) {
        fname_pheno_tsv = file.path(getwd(), paste0("simulated_phenotype-", date, ".tsv"))
        write.table(df, file=fname_pheno_tsv, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        if (verbose) {
            print("Output phenotype file (tsv):")
            print(fname_pheno_tsv)
        }
    }
    if (save_geno_vcf) {
        vcf = fn_G_to_vcf(G=G, min_depth=1e4, max_depth=1e4, verbose=verbose)
        if (class(vcf) == "gpError") {
            error = new("gpError",
                code=000,
                message=paste0(
                    "Error in load::fn_simulate_data(...). ",
                    "Please set n_alleles=2 as we can only convert biallelic loci into VCF format at the moment."
                ))
            error = chain(vcf, error)
            return(error)
        }
        fname_geno_vcf = file.path(getwd(), paste0("simulated_genotype-", date, ".vcf.gz"))
        vcfR::write.vcf(vcf, file=fname_geno_vcf)
        if (verbose) {
            print("Output genotype file (sync.gz):")
            print(fname_geno_vcf)
        }
    }
    if (save_geno_tsv)  {
        list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=G, verbose=verbose)
        df_geno = data.frame(chr=list_ids_chr_pos_all$vec_chr, pos=list_ids_chr_pos_all$vec_pos, allele=list_ids_chr_pos_all$vec_all, t(G))
        fname_geno_tsv = file.path(getwd(), paste0("simulated_genotype-", date, ".tsv"))
        write.table(df_geno, file=fname_geno_tsv, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        if (verbose) {
            print("Output genotype file (tsv):")
            print(fname_geno_tsv)
        }
    }
    if (save_geno_rds) {
        ### Convert loci names so that the chromosome, position and allele are tab-delimited
        fname_geno_rds = file.path(getwd(), paste0("simulated_genotype-", date, ".Rds"))
        ### Revert to numeric if we have more than 52 alleles (limit of the Latin alphabet)
        if (!non_numeric_Rds) {
            saveRDS(G, file=fname_geno_rds)
        } else {
            G_non_numeric = fn_G_numeric_to_non_numeric(G=G, ploidy=ploidy, verbose=verbose)
            saveRDS(G_non_numeric, file=fname_geno_rds)
        }
        if (verbose) {
            print("Output genotype file (Rds):")
            print(fname_geno_rds)
        }
    }
    return(list(
        fname_geno_vcf=fname_geno_vcf,
        fname_geno_tsv=fname_geno_tsv,
        fname_geno_rds=fname_geno_rds,
        fname_pheno_tsv=fname_pheno_tsv
    ))
}

#' Load genotype data as an allele frequency matrix
#' Load genotype data as an allele frequency matrix
#' 
#' @param fname_geno filename of the genotype file
#'  This may be in:
#'      - VCF format, 
#'      - allele frequency table saved as a tab-delimited file with a header line and the first 3 columns refer to the
#'          chromosome (chr), position (pos), and allele (allele),
#'          with subsequent columns referring to the allele frequencies of a sample, entry or pool.
#'          Names of the samples, entries, or pools in the header line can be any unique string of characters.
#' @param ploidy ploidy level which will generate genotype classes instead of continuous allele frequencies.
#'  If NULL, then continuous allele frequencies and no binning or classification into genotype classes
#'  will be performed [Default=NULL].
#' @param retain_minus_one_alleles_per_locus omit the alternative or trailing allele per locus? [Default=TRUE]
#' @param verbose show genotype loading messages? [Default=FALSE]
#' @returns
#' Ok: numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#' Err: gpError
#' @examples
#' list_sim = fn_simulate_data(verbose=TRUE, save_geno_vcf=TRUE, save_geno_tsv=TRUE, save_geno_rds=TRUE, save_pheno_tsv=TRUE)
#' G_vcf = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf, verbose=TRUE)
#' G_tsv = fn_load_genotype(fname_geno=list_sim$fname_geno_tsv, verbose=TRUE)
#' G_rds = fn_load_genotype(fname_geno=list_sim$fname_geno_rds, verbose=TRUE)
#' @export
fn_load_genotype = function(fname_geno, ploidy=NULL, retain_minus_one_alleles_per_locus=TRUE, verbose=FALSE) {
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(ploidy=8, verbose=TRUE, save_geno_vcf=TRUE, save_geno_tsv=TRUE, save_geno_rds=TRUE, save_pheno_tsv=TRUE)
    # # list_sim = fn_simulate_data(ploidy=8, verbose=TRUE, save_geno_vcf=TRUE, save_geno_tsv=TRUE, save_geno_rds=TRUE, save_pheno_tsv=TRUE, non_numeric_Rds=TRUE)
    # fname_geno = list_sim$fname_geno_vcf
    # fname_geno = list_sim$fname_geno_tsv
    # fname_geno = list_sim$fname_geno_rds
    # ploidy = 4
    # retain_minus_one_alleles_per_locus = TRUE
    # verbose = TRUE
    ###################################################
    ### Load the genotype matrix (n x p)
    ### TryCatch Rds, vcf, then tsv
    G = tryCatch({
        ################
        ### RDS file ###
        ################
        G = readRDS(fname_geno)
        if (is.numeric(G)==FALSE) {
            fn_G_non_numeric_to_numeric(G_non_numeric=G, verbose=verbose)
        }
        return(G)
    }, 
    error=function(e) {
        tryCatch({
            ################
            ### VCF file ###
            ################
            vcf = vcfR::read.vcfR(fname_geno, verbose=TRUE)
            G = fn_vcf_to_G(vcf, retain_minus_one_alleles_per_locus=retain_minus_one_alleles_per_locus, verbose=verbose)
            if (class(G)[1] == "gpError") {
                error = chain(G, new("gpError",
                    code=000,
                                        message=paste0(
                        "Error in load::fn_load_genotype(...).", 
                        "Error loading the vcf file: ", fname_geno, ".")))
                return(error)
            } else {
                return(G)
            }
        }, 
        error=function(e) {
            ########################################
            ### TSV: allele frequency table file ###
            ########################################
            df = read.delim(fname_geno, sep="\t", header=TRUE)
            if (sum(colnames(df)[1:3] == c("chr", "pos", "allele"), na.rm=TRUE) < 3) {
                error = new("gpError", 
                    code=000,
                    message=paste0(
                        "Error in load::fn_load_genotype(...). ",
                        "The file: ", fname_geno, " is not in allele frequency table format as described in the README.md. ",
                        "The first 3 columns do not correspond to 'chr', 'pos', and 'allele'."))
                return(error)
            }
            vec_loci_names = paste(df$chr, df$pos, df$allele, sep="\t")
            vec_entries = colnames(df)[c(-1:-3)]
            G = as.matrix(t(df[, c(-1:-3)]))
            rownames(G) = vec_entries
            colnames(G) = vec_loci_names
            return(G)
        })
    })
    ### Retain a-1 allele/s per locus, i.e. remove duplicates assuming all loci are biallelic
    if (retain_minus_one_alleles_per_locus==TRUE) {
        G = fn_G_split_off_alternative_allele(G=G, verbose=verbose)$G
    }
    return(G)
}

#' Filter genotypes by minimum allele frequency and minimum variation within locus (in an attempt to prevent duplicating the intercept)
#' and optionally by loci identities using a SNP which specifies the SNP coordinates as well as the reference and alternative alleles
#'
#' @param G numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#' @param maf minimum allele frequency [Default=0.01]
#' @param sdev_min minimum allele frequency standard deviation (constant allele frequency across samples is just another intercept) [Default=0.0001]
#' @param fname_snp_list name of the file containing the list of expected SNPs including their coordinates and alleles.
#'  This is a tab-delimited file with 3 columns: '#CHROM', 'POS', 'REF,ALT', corresponding to [Default=NULL]
#'  chromosome names (e.g. 'chr1' & 'chrom_1'), 
#'  numeric positions (e.g. 12345 & 100001), and 
#' reference-alternative allele strings separated by a comma (e.g. 'A,T' & 'allele_1,allele_alt') [Default=NULL]
#' @returns
#' Ok: numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#' Err: gpError
#' @examples
#' list_sim = fn_simulate_data(verbose=TRUE)
#' G = fn_load_genotype(list_sim$fname_geno_vcf, retain_minus_one_alleles_per_locus=FALSE)
#' maf = 0.05
#' sdev_min = 0.0001
#' verbose = TRUE
#' ### Simulate SNP list for filtering
#' mat_loci = matrix(unlist(strsplit(colnames(G), "\t")), byrow=TRUE, ncol=3)
#' vec_loci = unique(paste0(mat_loci[,1], "\t", mat_loci[,2]))
#' mat_loci = matrix(unlist(strsplit(vec_loci, "\t")), byrow=TRUE, ncol=2)
#' df_snp_list = data.frame(CHROM=mat_loci[,1], POS=as.numeric(mat_loci[,2]), REF_ALT=paste0("allele_1,allele_alt"))
#' df_snp_list$REF_ALT[1:100] = "allele_2,allele_4"
#' colnames(df_snp_list) = c("#CHROM", "POS", "REF,ALT")
#' fname_snp_list = "tmp_snp_list.txt"
#' write.table(df_snp_list, file=fname_snp_list, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
#' ### Filter
#' G_filtered = fn_filter_loci(G=G, maf=0.05, fname_snp_list=fname_snp_list, verbose=TRUE)
#' @export
fn_filter_loci = function(G, maf=0.01, sdev_min=0.0001, fname_snp_list=NULL, verbose=FALSE) {
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(verbose=TRUE)
    # G = fn_load_genotype(list_sim$fname_geno_vcf)
    # maf = 0.05
    # sdev_min = 0.0001
    # fname_snp_list = NULL
    # # ### Simulate SNP list for filtering
    # # G = fn_load_genotype(list_sim$fname_geno_vcf, retain_minus_one_alleles_per_locus=FALSE)
    # # mat_loci = matrix(unlist(strsplit(colnames(G), "\t")), byrow=TRUE, ncol=3)
    # # vec_loci = unique(paste0(mat_loci[,1], "\t", mat_loci[,2]))
    # # mat_loci = matrix(unlist(strsplit(vec_loci, "\t")), byrow=TRUE, ncol=2)
    # # df_snp_list = data.frame(CHROM=mat_loci[,1], POS=as.numeric(mat_loci[,2]), REF_ALT=paste0("allele_1,allele_alt"))
    # # df_snp_list$REF_ALT[1:100] = "allele_2,allele_4"
    # # colnames(df_snp_list) = c("#CHROM", "POS", "REF,ALT")
    # # fname_snp_list = "tmp_snp_list.txt"
    # # write.table(df_snp_list, file=fname_snp_list, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    # verbose = TRUE
    ###################################################
    ### Make sure the input thresholds are sensible
    if ((maf < 0.0) | (maf > 1.0)) {
        error = new("gpError",
            code=000,
            message=paste0(
                "Error in load::fn_filter_loci(...). ",
                "Please use a minimum allele frequency (maf) between 0 and 1."
            ))
        return(error)
    }
    if ((sdev_min < 0.0) | (sdev_min > 1.0)) {
        error = new("gpError",
            code=000,
            message=paste0(
                "Error in load::fn_filter_loci(...). ",
                "Please use a minimum standard deviation in allele frequency (sdev_min) between 0 and 1."
            ))
        return(error)
    }
    ### Extract the mean and standard deviation in allele frequencies per locus
    vec_freqs = colMeans(G, na.rm=TRUE)
    vec_sdevs = apply(G, MARGIN=2, FUN=sd)
    if (verbose) {
        print("Distribution of mean allele frequencies per locus:")
        txtplot::txtdensity(vec_freqs)
        print("Distribution of allele frequency standard deviation per locus:")
        txtplot::txtdensity(vec_sdevs)
    }
    ### Filter by minimum allele frequency and minimum allele frequency variance (constant allele frequency across samples is just another intercept)
    vec_idx = which(
        (vec_freqs >= maf) &
        (vec_freqs <= (1-maf)) &
        (vec_sdevs >= sdev_min))
    if (length(vec_idx) < ncol(G)) {
        G = G[, vec_idx, drop=FALSE]
    } else {
        if (verbose) {print("All loci passed the minimum allele frequency and standard deviation thresholds.")}
    }
    ### Filter using a SNP list
    if (!is.null(fname_snp_list)) {
        df = read.delim(fname_snp_list, sep="\t", header=TRUE, check.names=FALSE)
        if (sum(colnames(df) != c("#CHROM", "POS", "REF,ALT")) > 0) {
            error = new("gpError",
                code=000,
                message=paste0(
                    "Error in load::fn_filter_loci(...). ",
                    "The SNP list file: ", fname_snp_list, " does not have the expected field names: ",
                    "'#CHROM', 'POS', and 'REF,ALT'."
                ))
            return(error)
        }
        mat_ref_alt = matrix(unlist(strsplit(df$`REF,ALT`, ",")), byrow=TRUE, ncol=2)
        df$REF = mat_ref_alt[,1]
        df$ALT = mat_ref_alt[,2]
        vec_expected_snps = c(
            paste(df$`#CHROM`, df$POS, df$REF, sep="\t"), 
            paste(df$`#CHROM`, df$POS, df$ALT, sep="\t"))
        vec_colnames = colnames(G)
        p = length(vec_colnames)
        vec_idx = c()
        if (verbose) {pb = txtProgressBar(min=0, max=p, style=3)}
        for (j in 1:p) {
            # j = 1
            if (sum(vec_expected_snps %in% vec_colnames[j]) == 1) {
                vec_idx = c(vec_idx, j)
            }
            if (verbose) {setTxtProgressBar(pb, j)}
        }
        if (verbose) {close(pb)}
        if (length(vec_idx) < p) {
            G = G[, vec_idx, drop=FALSE]
        } else {
            if (verbose) {print("All loci passed the minimum allele frequency and standard deviation thresholds.")}
        }
    }
    ### Return filtered allele frequency matrix
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

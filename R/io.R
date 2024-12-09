suppressWarnings(suppressPackageStartupMessages(library(simquantgen)))
suppressWarnings(suppressPackageStartupMessages(library(vcfR)))
suppressWarnings(suppressPackageStartupMessages(library(txtplot)))

### Create an error class
methods::setClass("gpError", representation(code="numeric", message="character"), prototype(code=0, message="Empty error message. Please define."))
### Create an error chaining method
methods::setGeneric("chain", function(x, y){
    standardGeneric("chain")
})
methods::setMethod(f="chain",
  signature=c(x="gpError", y="gpError"),
  function(x, y) {
    codes = c(x@code, y@code)
    messages = c(x@message, y@message)
    return(methods::new("gpError", code=codes, message=messages))
  }
)
# err1 = methods::new("gpError", code=1, message="message 1")
# err2 = methods::new("gpError", code=2, message="message 2")
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
#' @param verbose show names extraction messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'      + $vec_ids: vector of sample/entry/pool names (same order as in G)
#'      + $vec_chr: vector of chromosome or scaffold names (same order as in G)
#'      + $vec_pos: vector of numeric positions per chromosome (same order as in G)
#'      + $vec_all: vector of allele names (same order as in G)
#'  - Err: grError
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
    if (sum(duplicated(vec_ids)) > 0) {
        error = methods::new("gpError",
            code=200,
            message=paste0(
                "Error in io::fn_G_extract_names(...). ",
                "The sample names (row names) have duplicates. ",
                "We expect unique samples in the genotype file. ",
                "Please remove or generate consensus among duplicated genotypes."))
        return(error)
    }
    if (sum(duplicated(vec_loci)) > 0) {
        error = methods::new("gpError",
            code=201,
            message=paste0(
                "Error in io::fn_G_extract_names(...). ",
                "The loci names (column names) have duplicates. ",
                "Please remove the duplicated loci."))
        return(error)
    }
    n_identifiers = length(unlist(strsplit(vec_loci[1], "\t"))) ### Number of loci identifiers where we expect at least 2
    if (n_identifiers < 2) {
        error = methods::new("gpError",
            code=202,
            message=paste0(
                "Error in io::fn_G_extract_names(...). ",
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
        error = methods::new("gpError",
            code=203,
            message=paste0(
                "Error in io::fn_G_extract_names(...). ",
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

#' Split trailing alternative allele from a numeric (allele frequency) matrix if it exists
#' 
#' @param G n samples x p loci matrix of allele frequencies (numeric ranging from 0 to 1) with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#' @param verbose show genotype data alternative allele splitting messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'      + $G: n samples x p loci matrix of allele frequencies excluding trailing allele (numeric ranging from 0 to 1)
#'      + $G_alt: n samples x p loci matrix of allele frequencies of the trailing allele (numeric ranging from 0 to 1)
#'  - Err: grError
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
    ### Input sanity check
    if (methods::is(G, "gpError")) {
        error = chain(G, methods::new("gpError",
            code=204,
            message=paste0(
                "Error in io::fn_G_split_off_alternative_allele(...). ",
                "Input G is an error type."
            )))
        return(error)
    }
    ### Make sure G contains allele frequencies
    if (sum(((G < 0) | (G > 1) | is.infinite(G)), na.rm=TRUE) > 0) {
        error = methods::new("gpError",
            code=205,
            message=paste0(
                "Error in io::fn_G_split_off_alternative_allele(...). ",
                "We are expecting a matrix allele frequencies but ",
                "we are getting negative values and/or values greater than 1 and/or infinite values."))
        return(error)
    }
    ### Extract sample/entry/pool, and loci names
    list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=G, verbose=verbose)
    if (methods::is(list_ids_chr_pos_all, "gpError")) {
        error = chain(list_ids_chr_pos_all, methods::new("gpError",
            code=206,
            message=paste0(
                "Error in io::fn_G_split_off_alternative_allele(...). ",
                "Error type returned by list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=G, verbose=verbose)."
            )))
        return(error)
    }
    ### Iterate across loci removing the trailing alternative allele
    vec_loci = sort(unique(paste0(list_ids_chr_pos_all$vec_chr, "\t", list_ids_chr_pos_all$vec_pos)))
    p = length(vec_loci)
    n = nrow(G)
    vec_idx_loci = c()
    vec_idx_alt = c()
    if (ncol(G) > p) {
        if (verbose) {print("Removing trailing alternative alleles per locus:")}
        if (verbose) {pb = utils::txtProgressBar(min=0, max=p, style=3)}
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
                vec_idx_alt = c(vec_idx_alt, utils::tail(idx, n=1))
            } else {
                ### Otherwise keep all the loci if all frequencies do not sum up to 1
                vec_idx_loci = c(vec_idx_loci, idx)
            }
            if (verbose) {utils::setTxtProgressBar(pb, j)}
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

#' Convert a numeric (allele frequency) genotype matrix into a non-numeric (genotype classes) genotype matrix
#' 
#' @param G n samples x p loci matrix of allele frequencies (numeric ranging from 0 to 1) with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#' @param ploidy ploidy level which can refer to the number of haploid genomes to simulate pools (Default=2)
#' @param verbose show numeric to non-numeric genotype data conversion messages? (Default=FALSE)
#' @returns
#'  - Ok: n samples x p loci matrix of genotype classes
#'  - Err: grError
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
    ### Input sanity check
    if (methods::is(G, "gpError")) {
        error = chain(G, methods::new("gpError",
            code=207,
            message=paste0(
                "Error in io::fn_G_numeric_to_non_numeric(...). ",
                "Input G is an error type."
            )))
        return(error)
    }
    ### Make sure G contains allele frequencies
    if (sum(((G < 0) | (G > 1) | is.infinite(G)), na.rm=TRUE) > 0) {
        error = methods::new("gpError",
            code=208,
            message=paste0(
                "Error in io::fn_G_numeric_to_non_numeric(...). ",
                "We are expecting a matrix allele frequencies but ",
                "we are getting negative values and/or values greater than 1 and/or infinite values."))
        return(error)
    }
    if (ploidy < 1) {
        error = methods::new("gpError",
            code=209,
            message=paste0(
                "Error in io::fn_G_numeric_to_non_numeric(...). ",
                "Ploidy cannot be less than 1."))
        return(error)
    }
    if (ploidy != round(ploidy)) {
        error = methods::new("gpError", 
            code=210,
            message=paste0(
                "Error in io::fn_G_numeric_to_non_numeric(...). ",
                "Please pick a positive integer as the ploidy instead of ", ploidy, "."))
        return(error)
    }
    ### Extract sample/entry/pool, and loci names
    list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=G, verbose=verbose)
    if (methods::is(list_ids_chr_pos_all, "gpError")) {
        error = chain(list_ids_chr_pos_all, methods::new("gpError",
            code=211,
            message=paste0(
                "Error in io::fn_G_numeric_to_non_numeric(...). ",
                "Error type returned by list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=G, verbose=verbose)."
            )))
        return(error)
    }
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
    if (verbose) {pb = utils::txtProgressBar(min=0, max=p, style=3)}
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
        if (verbose) {utils::setTxtProgressBar(pb, j)}
    }
    if (verbose) {close(pb)}
    return(G_non_numeric)
}

#' Convert a non-numeric genotype classes matrix into a numeric (allele frequency) genotype matrix which includes all the alleles
#' 
#' @param G_non_numeric n samples x p loci matrix of genotype classes with non-null row and column names.
#'  The genotype classes are represented by uppercase and/or lowercase letters,
#'  where the total number of letters correspond to the ploidy level,
#'  and the number of unique letters refer to the number of alleles per locus.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  and the second should be numeric which refers to the position in the chromosome/scaffold.
#' @param retain_minus_one_alleles_per_locus omit the alternative or trailing allele per locus? (Default=TRUE)
#' @param verbose show non-numeric to numeric genotype data conversion messages? (Default=FALSE)
#' @returns
#'  - Ok: n samples x p loci-alleles matrix of genotype classes (numeric ranging from 0 to 1). If the column names do not have
#'      allele IDs then the alleles will be extracted from the non-numeric genotype encodings.
#'  - Err: grError
#' @examples
#' ploidy = 42
#' G_numeric = simquantgen::fn_simulate_genotypes(ploidy=ploidy, n_alleles=52, verbose=TRUE)
#' G_non_numeric = fn_G_numeric_to_non_numeric(G=G_numeric, ploidy=ploidy, verbose=TRUE)
#' G_numeric_back = fn_G_non_numeric_to_numeric(G=G_non_numeric, verbose=TRUE)
#' @export
fn_G_non_numeric_to_numeric = function(G_non_numeric, retain_minus_one_alleles_per_locus=TRUE, verbose=FALSE) {
    ###################################################
    ### TEST
    # G = simquantgen::fn_simulate_genotypes(ploidy=42, n_alleles=52, verbose=TRUE)
    # G_non_numeric = gp::fn_G_numeric_to_non_numeric(G=G, ploidy=42, verbose=TRUE)
    # verbose = TRUE
    ###################################################
    ### Input sanity check
    if (methods::is(G_non_numeric, "gpError")) {
        error = chain(G_non_numeric, methods::new("gpError",
            code=212,
            message=paste0(
                "Error in io::fn_G_non_numeric_to_numeric(...). ",
                "Input G_non_numeric is an error type."
            )))
        return(error)
    }
    ### Extract ploidy where we assume the same ploidy across the entire data set
    ploidy = length(unlist(strsplit(G_non_numeric[1,1], "")))
    for (i in sample(1:nrow(G_non_numeric), size=min(c(100, nrow(G_non_numeric))), replace=FALSE)) {
        for (j in sample(1:ncol(G_non_numeric), size=min(c(100, ncol(G_non_numeric))), replace=FALSE)) {
            if (length(unlist(strsplit(G_non_numeric[i, j], ""))) != ploidy) {
                error = methods::new("gpError",
                    code=213,
                    message=paste0(
                        "Error in io::fn_G_non_numeric_to_numeric(...). ",
                        "The ploidy level is not consistent across the data set. ",
                        "From the first locus of the first sample we expected a ploidy level of ", ploidy, "X."))
                return(error)
            }
        }
    }
    ### Do not update the allele IDs in the column names if there are already allele names, i.e. there are 3 elements per column name after splitting by tabs
    bool_update_allele_ids = TRUE
    if (length(unlist(strsplit(colnames(G_non_numeric)[1], "\t"))) == 3) {
        bool_update_allele_ids = FALSE
    }
    ### Extract the loci and alleles per locus
    vec_loci_names = colnames(G_non_numeric)
    list_allele_names = list()
    p = 0
    if (verbose) {pb = utils::txtProgressBar(min=0, max=ncol(G_non_numeric), style=3)}
    for (j in 1:ncol(G_non_numeric)) {
        # j = 1
        vec_alleles = sort(unique(unlist(strsplit(G_non_numeric[, j], ""))))
        eval(parse(text=paste0("list_allele_names$`", vec_loci_names[j], "` = vec_alleles")))
        p = p + length(vec_alleles)
        if (verbose) {utils::setTxtProgressBar(pb, j)}
    }
    if (verbose) {close(pb)}
    ### Instantiate the numeric allele frequency genotype matrix
    n = nrow(G_non_numeric)
    G = matrix(0, nrow=n, ncol=p)
    rownames(G) = rownames(G_non_numeric)
    ### Populate
    vec_colnames = c()
    idx_locus_allele = 0
    if (verbose) {pb = utils::txtProgressBar(min=0, max=ncol(G_non_numeric), style=3)}
    for (j in 1:ncol(G_non_numeric)) {
        # j = 1
        vec_alleles = list_allele_names[[j]]
        list_geno_classes = strsplit(G_non_numeric[, j], "")
        for (k in 1:length(vec_alleles)) {
            # k = 1
            idx_locus_allele = idx_locus_allele + 1
            ### Add the allele ID if the column names lack them
            ### But note that we are dividing the loci into its constituent alleles hence we are only properly labelling the first allele.
            if (bool_update_allele_ids | (k > 1)) {
                if (bool_update_allele_ids) {
                    vec_colnames = c(vec_colnames, paste0(vec_loci_names[j], "\t", vec_alleles[k]))
                } else {
                    ### Make sure we are extracting the chromosom and position info only for the column names with allele IDs
                    chrom_pos = paste(unlist(strsplit(vec_loci_names[j], "\t"))[1:2], collapse="\t")
                    vec_colnames = c(vec_colnames, paste0(chrom_pos, "\t", vec_alleles[k]))
                }
            } else {
                vec_colnames = c(vec_colnames, vec_loci_names[j])
            }
            for (i in 1:n) {
                # i = 1
                G[i, idx_locus_allele] = sum(list_geno_classes[[i]] == vec_alleles[k]) / ploidy
            }
        }
        if (verbose) {utils::setTxtProgressBar(pb, j)}
    }
    if (verbose) {close(pb)}
    ### Update the loci names
    colnames(G) = vec_colnames
    ### Return numeric allele frequency genotype matrix with or without all the alleles
    if (retain_minus_one_alleles_per_locus) {
        list_G_G_alt = fn_G_split_off_alternative_allele(G=G, verbose=verbose)
        if (methods::is(list_G_G_alt, "gpError")) {
            error = chain(list_G_G_alt, methods::new("gpError",
                code=214,
                message=paste0(
                    "Error in io::fn_G_non_numeric_to_numeric(...). ",
                    "Error type returned by list_G_G_alt = fn_G_split_off_alternative_allele(G=G, verbose=verbose)."
                )))
            return(error)
        }
    }
    ### Output
    return(list_G_G_alt$G)
}

#' Convert numeric allele frequency matrix into a vcfR::vcf object with randomly sampled depths (for fn_simulate_data(...) below)
#' 
#' @param G numeric n samples x p biallelic loci matrix of allele frequencies with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'  Note that this function only accepts biallelic loci.
#' @param min_depth minimum depth per locus for sampling simulated depths (Default=100)
#' @param max_depth maximum depth per locus for sampling simulated depths (Default=1000)
#' @param verbose show allele frequency genotype matrix to vcf conversion messages? (Default=FALSE)
#' @returns
#'  - Ok: simulated genotype data as a vcfR object with GT, AD and DP fields
#'  - Err: grError
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
    ### Input sanity check
    if (methods::is(G, "gpError")) {
        error = chain(G, methods::new("gpError",
            code=215,
            message=paste0(
                "Error in io::fn_G_to_vcf(...). ",
                "Input G is an error type."
            )))
        return(error)
    }
    ### Make sure G contains allele frequencies
    if (sum(((G < 0) | (G > 1) | is.infinite(G)), na.rm=TRUE) > 0) {
        error = methods::new("gpError",
            code=216,
            message=paste0(
                "Error in io::fn_G_to_vcf(...). ",
                "We are expecting a matrix allele frequencies but ",
                "we are getting negative values and/or values greater than 1 and/or infinite values."))
        return(error)
    }
    if (min_depth > max_depth) {
        error = methods::new("gpError",
            code=217,
            message=paste0(
                "Error in io::fn_G_to_vcf(...). ",
                "Minimum depth (", min_depth, ") cannot be greater than the maximum depth (", max_depth, ")."))
        return(error)
    }
    ### Split-off the alternative allele from the reference allele
    list_G_G_alt = fn_G_split_off_alternative_allele(G=G, verbose=verbose)
    if (methods::is(list_G_G_alt, "gpError")) {
        error = chain(list_G_G_alt, methods::new("gpError",
            code=218,
            message=paste0(
                "Error in io::fn_G_to_vcf(...). ",
                "Error type returned by list_G_G_alt = fn_G_split_off_alternative_allele(G=G, verbose=verbose)."
            )))
        return(error)
    }
    G_alt = list_G_G_alt$G_alt
    G = list_G_G_alt$G
    n = nrow(G)
    p = ncol(G)
    ### Extract the names of the entries/samples/pools and loci
    list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=G, verbose=verbose)
    if (methods::is(list_ids_chr_pos_all, "gpError")) {
        error = chain(list_ids_chr_pos_all, methods::new("gpError",
            code=219,
            message=paste0(
                "Error in io::fn_G_to_vcf(...). ",
                "Error type returned by list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=G, verbose=verbose)."
            )))
        return(error)
    }
    ### Make sure the input matrix is biallelic
    ### Also, convert the tabs in the loci-alleles names into dashes so as not to interfere with the VCF format
    vec_loci_names = paste0(list_ids_chr_pos_all$vec_chr, "_", list_ids_chr_pos_all$vec_pos)
    vec_allele_counts_per_locus = table(vec_loci_names)
    if (((sum(vec_allele_counts_per_locus > 2) > 0) & is.null(G_alt)) | ((sum(vec_allele_counts_per_locus >= 2) > 0) & !is.null(G_alt))) {
        error = methods::new("gpError",
            code=220,
            message=paste0(
                "Error in io::fn_G_to_vcf(...). ",
                "Apologies, as this function at the moment can only convert biallelic allele frequency matrices into VCF format."))
        return(error)
    }
    ### Extract the names of the alternative alleles
    if (is.null(G_alt)) {
        # print("WARNING! Alternative alleles are unknown. Making psuedo-alleles.")
        # vec_alleles_list = c("A", "T", "C", "G")
        # vec_alt = unlist(lapply(list_ids_chr_pos_all$vec_all, FUN=function(x){ sample(size=1, x=vec_alleles_list[grep(x, vec_alleles_list, ignore.case=TRUE, invert=TRUE)]) }))
        vec_alt = rep("N", times=ncol(G))
        # vec_alt = rep("allele_alt", times=ncol(G))
    } else {
        ### Convert the tabs in the loci-alleles names into dashes so as not to interfere with the VCF format
        vec_alt = colnames(gsub("\t", "_", list_G_G_alt$G_alt))
    }
    ### Populate required vcf fields with loci identities
    META = c(
        "##fileformat=VCFv", 
        paste0("##", unique(list_ids_chr_pos_all$vec_chr)), 
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">',
        '##Extracted from text file.')
    FIX = cbind(
        list_ids_chr_pos_all$vec_chr, 
        list_ids_chr_pos_all$vec_pos, 
        vec_loci_names, 
        list_ids_chr_pos_all$vec_all, 
        vec_alt, 
        rep(".", each=p),
        rep("PASS", each=p),
        rep(".", each=p))
    colnames(FIX) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    ### Population the GT:AD:DP matrix with strings of genotype class and allele depths
    GT_AD_DP = matrix("", nrow=p, ncol=(n+1))
    GT_AD_DP[,1] = "GT:AD:DP"
    colnames(GT_AD_DP) = c("FORMAT", list_ids_chr_pos_all$vec_ids)
    if (verbose) {pb = utils::txtProgressBar(min=0, max=n, style=3)}
    for (i in 1:n) {
        for (j in 1:p) {
            ### We're assuming that the counts in G are for the reference alleles
            g = round(2 * G[i, j])
            if (is.na(g)) {
                GT_AD_DP[j, (i+1)] = "./."
            } else {
                if (g == 0) {
                    GT_AD_DP[j, (i+1)] = "1/1"
                } else if (g == 1) {
                    GT_AD_DP[j, (i+1)] = "0/1"
                } else if (g == 2) {
                    GT_AD_DP[j, (i+1)] = "0/0"
                } else {
                    GT_AD_DP[j, (i+1)] = "./."
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
            GT_AD_DP[j, (i+1)] = paste0(GT_AD_DP[j, (i+1)], ":", ref_depth, ",", alt_depth, ":", tot_depth)
        }
        if (verbose) {utils::setTxtProgressBar(pb, i)}
    }
    if (verbose) {close(pb)}
    ### Load dummy vcf file from pinfsc50 package which comes with vcfR as one of its dependencies
    vcf_dummy = vcfR::read.vcfR(system.file("extdata", "pinf_sc50.vcf.gz", package = "pinfsc50"), verbose=FALSE)
    ### Create the new vcfR object
    vcf = vcf_dummy
    vcf@meta = META
    vcf@fix = FIX
    vcf@gt = GT_AD_DP
    if (verbose) {print(vcf)}
    return(vcf)
}

#' Convert biallelic vcf data into allele frequencies where loci beyond minimum and maximum depths are set to missing
#' 
#' @param vcf biallelic genotype data as a vcfR object with GT and/or AD and DP fields (where AD takes precedence over GT)
#' @param min_depth minimum depth per locus beyond which will be set to missing data (Default=0)
#' @param max_depth maximum depth per locus beyond which will be set to missing data (Default=Inf)
#' @param force_biallelic assume biallelic loci thereby dropping extra rows in vcf files corresponding to multi-allelic loci.
#'  Generate the vcf via `bcftools call -m ...` for multi-allelic and rare-variant calling,
#'  followed by `bcftools -m - ...` to  split multi-allelic sites into biallelic records (-) (Default=TRUE)
#' @param retain_minus_one_alleles_per_locus omit the alternative or trailing allele per locus? (Default=TRUE)
#' @param verbose show vcf to allele frequency genotype matrix conversion messages? (Default=FALSE)
#' @returns
#'  - Ok: 
#'      + G: named n samples x p biallelic loci matrix of allele frequencies.
#'      Row names can be any string of characters.
#'      Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'      the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'      subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'      + D: named n samples x p biallelic loci matrix of sequencing depth per locus with the same
#'      row and column naming conventions as G but only includes the names of the reference alleles per locus.
#'  - Err: gpError
#' @examples
#' G = simquantgen::fn_simulate_genotypes(verbose=TRUE)
#' vcf = fn_G_to_vcf(G=G, verbose=TRUE)
#' list_G_D = fn_vcf_to_G(vcf=vcf, verbose=TRUE)
#' G_back = list_G_D$G
#' @export
fn_vcf_to_G = function(vcf, min_depth=0, max_depth=.Machine$integer.max, force_biallelic=TRUE, retain_minus_one_alleles_per_locus=TRUE, verbose=FALSE) {
    ###################################################
    ### TEST
    # G = simquantgen::fn_simulate_genotypes(verbose=TRUE)
    # vcf = fn_G_to_vcf(G=G, verbose=TRUE)
    # min_depth = 10
    # max_depth = 500
    # force_biallelic = TRUE
    # retain_minus_one_alleles_per_locus = TRUE
    # verbose = TRUE
    ###################################################
    ### Check input type
    if (!methods::is(vcf, "vcfR")) {
        error = methods::new("gpError", 
            code=221,
            message="Error in io::fn_vcf_to_G(vcf): vcf is not a vcfR object.")
        return(error)
    }
    ### Extract loci names
    vec_chr = vcfR::getCHROM(vcf)
    vec_pos = vcfR::getPOS(vcf)
    ### Remove additional alleles from multi-allelic sites if we're assuming biallelic loci
    if (force_biallelic) {
        ### Removes instances of the same locus position except for the first instance,
        ### i.e. does not remove the multi-allelic loci but rather forces biallelic-ness at all loci
        vec_idx = which(!duplicated(paste(vec_chr, vec_pos, sep="\t")))
        if (length(vec_idx) == 0) {
            error = methods::new("gpError", 
                code=222,
                message=paste0(
                    "Error in io::fn_vcf_to_G(vcf). ",
                    "This is impossible as the first instance of duplicated elements will be included."
                ))
            return(error)
        } else if (length(vec_idx) < length(vec_chr)) {
            vcf = vcf[vec_idx, , ]
            if (verbose) {
                print(paste0("Removing the third, fourth, fifth, etc alleles in multi-allelic loci."))
                print(paste0("Omitting a total of ", length(vec_chr)-length(vec_idx), " alleles across multi-allelic loci."))
                print(vcf)
            }
            ### Re-extract loci names after omitting extra alleles
            vec_chr = vcfR::getCHROM(vcf)
            vec_pos = vcfR::getPOS(vcf)
        } else {
            if (verbose) {print("All loci are biallelic.")}
        }
    }
    ### Extract loci and pool/sample names
    vec_loci_ref_names = paste(vec_chr, vec_pos, vcfR::getREF(vcf), sep="\t")
    if (!retain_minus_one_alleles_per_locus) {
        vec_loci_alt_names = paste(vec_chr, vec_pos, vcfR::getALT(vcf), sep="\t")
    }
    vec_pool_names = colnames(vcf@gt)[-1]
    vec_elements = unique(vcf@gt[, 1])
    ### Make sure that all the loci consistently have either the AD or GT field across all type of field combinations, e.g. c(GT:AD, AD, GT:PL:AD) where the AD field is consistenly present
    bool_consistent_AD = sum(grepl("AD", vec_elements)) == length(vec_elements)
    bool_consistent_GT = sum(grepl("GT", vec_elements)) == length(vec_elements)
    if (!(bool_consistent_AD | bool_consistent_GT)) {
        error = methods::new("gpError", 
            code=223,
            message=paste0(
                "Error in io::fn_vcf_to_G(vcf): please check the format of your input vcf file. ",
                "The same fields across loci is required. ",
                "Please pick one field architecture from the following =", paste(paste0("'", vec_elements, "'"), collapse=","), ". ",
                "Reformat your vcf file to have the same fields across loci."))
        return(error)
    }
    ### Also make sure that the GP field is present so that we can filter by depth
    bool_consistent_DP = sum(grepl("DP", vec_elements)) == length(vec_elements)
    if (!bool_consistent_DP) {
        error = methods::new("gpError", 
            code=224,
            message=paste0(
                "Error in io::fn_vcf_to_G(vcf): please check the format of your input vcf file. ",
                "Make sure the 'AD' and/or 'GT' and 'DP' fields are present. ",
                "These are the fields present in your vcf file: ", gsub(":", ", ", vec_elements), ". ",
                "Regenerate your vcf file to include the 'AD' field and/or 'GT' field."))
        return(error)
    }
    ### Set loci into missing if depth is below min_depth or above max_depth
    mat_depth = vcfR::extract.gt(vcf, element="DP", as.numeric=TRUE)
    if (verbose) {
        print("Distribution of allele depths:")
        vec_sample_depths = sample(mat_depth, size=min(1e4, c(prod(dim(mat_depth)))))
        txtplot::txtdensity(vec_sample_depths[!is.na(vec_sample_depths)])
    }
    mat_idx = (mat_depth < min_depth) | ((mat_depth > max_depth))
    if (sum(mat_idx, na.rm=TRUE) > 0) {
        if (verbose) {
            print(paste0("Before defining missing data by minimum and maximum depths: ", min_depth, " - ", max_depth, ":"))
            print(vcf)
        }
        vcf@gt[, -1][mat_idx] = NA
        mat_depth = vcfR::extract.gt(vcf, element="DP", as.numeric=TRUE)
        if (verbose) {
            print(paste0("After defining missing data by minimum and maximum depths:"))
            print(vcf)
            print("Distribution of allele depths after defining missing data by minimum and maximum depths:")
            vec_sample_depths = sample(mat_depth, size=min(1e4, c(prod(dim(mat_depth)))))
            txtplot::txtdensity(vec_sample_depths[!is.na(vec_sample_depths)])
        }
    }
    ### Extract genotype data where the AD field takes precedence over the GT field
    if (sum(grepl("AD", vec_elements)) == length(vec_elements)) {
        mat_allele_counts = vcfR::extract.gt(vcf, element="AD")
        if (length(unlist(strsplit(mat_allele_counts[1,1], ","))) > 2) {
            error = methods::new("gpError",
                code=225,
                message=paste0(
                    "Error in io::fn_vcf_to_G(...). ",
                    "Apologies because at the moment we can only convert biallelic VCF files into an allele frequency matrix."))
            return(error)
        }
        mat_ref_counts = vcfR::masplit(mat_allele_counts, delim=',', record=1, sort=0)
        mat_alt_counts = vcfR::masplit(mat_allele_counts, delim=',', record=2, sort=0)
        if (verbose) {
            print("Distribution of reference allele depths")
            txtplot::txtdensity(mat_ref_counts[!is.na(mat_ref_counts)])
            print("Distribution of alternative allele depths")
            txtplot::txtdensity(mat_alt_counts[!is.na(mat_alt_counts)])
        }
        ### Set missing allele counts to 0, if the other allele is non-missing and non-zero
        idx_for_ref = which(!is.na(mat_alt_counts) & (mat_alt_counts != 0) & is.na(mat_ref_counts))
        idx_for_alt = which(!is.na(mat_ref_counts) & (mat_ref_counts != 0) & is.na(mat_alt_counts))
        mat_ref_counts[idx_for_ref] = 0
        mat_alt_counts[idx_for_alt] = 0
        ### Calculate reference allele frequencies
        G = mat_ref_counts / (mat_ref_counts + mat_alt_counts)
    } else if (sum(grepl("GT", vec_elements)) == length(vec_elements)) {
        GT = vcfR::extract.gt(vcf, element="GT")
        G = matrix(NA, nrow=length(vec_loci_ref_names), ncol=length(vec_pool_names))
        G[(GT == "0/0") | (GT == "0|0")] = 1.0
        G[(GT == "1/1") | (GT == "1|1")] = 0.0
        G[(GT == "0/1") | (GT == "0|1") | (GT == "1|0")] = 0.5
    } else {
        error = methods::new("gpError", 
            code=226,
            message="Error in io::fn_vcf_to_G(vcf): vcf needs to have the 'AD' and/or 'GT' fields present.")
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
        print(paste0(c("q_min=", "q_max="), round(range(G, na.rm=TRUE), 4)))
        txtplot::txtdensity(G[!is.na(G)])
    }
    ### Output n x p matrix of allele frequencies and matrix of depths
    rownames(mat_depth) = vec_loci_ref_names
    colnames(mat_depth) = vec_pool_names
    return(list(
        G=t(G),
        D=t(mat_depth)))
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
#' @param ploidy ploidy level which can refer to the number of haploid genomes to simulate pools (Default=2)
#' @param verbose show genotype binning or classification messages? (Default=FALSE)
#' @returns
#'  - Ok: named n samples x p loci-alleles matrix of numeric genotype classes
#'  - Err: gpError
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
    ### Input sanity check
    if (methods::is(G, "gpError")) {
        error = chain(G, methods::new("gpError",
            code=227,
            message=paste0(
                "Error in io::fn_classify_allele_frequencies(...). ",
                "Input G is an error type."
            )))
        return(error)
    }
    ### Make sure G contains allele frequencies
    if (sum(((G < 0) | (G > 1) | is.infinite(G)), na.rm=TRUE) > 0) {
        error = methods::new("gpError",
            code=228,
            message=paste0(
                "Error in io::fn_classify_allele_frequencies(...). ",
                "We are expecting a matrix allele frequencies but ",
                "we are getting negative values and/or values greater than 1 and/or infinite values."))
        return(error)
    }
    if (ploidy < 1) {
        error = methods::new("gpError", 
            code=229,
            message=paste0(
                "Error in io::fn_classify_allele_frequencies(...): Are you sure the ploidy is ", ploidy, "X?",
                "How on this beautiful universe does that work?",
                "Please pick a positive integer!"))
        return(error)
    }
    if (ploidy != round(ploidy)) {
        error = methods::new("gpError", 
            code=230,
            message=paste0(
                "Error in io::fn_classify_allele_frequencies(...): Are you sure the ploidy is ", ploidy, "X?",
                "How on this beautiful universe does that work?",
                "Please pick a positive integer!"))
        return(error)
    }
    G_classes = round(G * ploidy) / ploidy
    ### Genotype classses distribution
    if (verbose) {
        print("Genotype classes distribution:")
        print(paste0(c("q_min=", "q_max="), round(range(G_classes), 4)))
        txtplot::txtdensity(G_classes[!is.na(G_classes)])
        print("Genotype classes:")
        print(table(G_classes))
    }
    ### Output n x p matrix of genotype classes
    return(G_classes)
}

#' Simple wrapper of simquantgen simulation of 10 standard normally distributed QTL additive effects 
#' on 5-chromosome genome and a single trait at 50% heritability
#' 
#' @param n number of samples (Default=100)
#' @param l number of loci (Default=1000)
#' @param ploidy ploidy level which can refer to the number of haploid genomes to simulate pools (Default=2)
#' @param n_alleles maximum number of alleles per locus (Default=2)
#' @param min_depth minimum depth per locus (Default=5)
#' @param max_depth maximum depth per locus (Default=5000)
#' @param n_pop number of randomly assigned population groupings (Default=1)
#' @param n_effects number of additive QTL effects to simulate (Default=10)
#' @param h2 narrow-sense heritability to simulate (Default=0.5)
#' @param seed randomisation seed for replicability (Default=12345)
#' @param save_pheno_tsv save the phenotype data as a tab-delimited file? (Default=TRUE)
#' @param save_geno_vcf save the genotype data as a vcf file? (Default=TRUE)
#' @param save_geno_tsv save the genotype data as a tab-delimited allele frequency table file? (Default=FALSE)
#' @param save_geno_rds save the named genotype matrix as an Rds file? (Default=FALSE)
#' @param non_numeric_Rds save non-numeric Rds genotype file? (Default=FALSE)
#' @param verbose show genotype and phenotype data simulation messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'      + $fname_geno_vcf: filename of the simulated genotype data as a vcf file
#'      + $fname_geno_tsv: filename of the simulated genotype data as a tab-delimited allele frequency table file
#'      + $fname_geno_rds: filename of the simulated named genotype matrix as an Rds file
#'      + $fname_pheno_tsv: filename of the simulated phenotype data as a tab-delimited file
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(verbose=TRUE, 
#'  save_geno_vcf=TRUE, save_geno_tsv=TRUE, save_geno_rds=TRUE, save_pheno_tsv=TRUE)
#' @export
fn_simulate_data = function(n=100, l=1000, ploidy=2, n_alleles=2, min_depth=5, max_depth=1000, n_pop=1, 
    n_effects=10, h2=0.5, seed=12345, 
    save_pheno_tsv=TRUE, save_geno_vcf=TRUE, save_geno_tsv=FALSE, save_geno_rds=FALSE, non_numeric_Rds=FALSE, verbose=FALSE) {
    ###################################################
    ### TEST
    # n = 100
    # l = 1000
    # ploidy = 42
    # n_alleles = 2
    # min_depth = 5
    # max_depth = 1000
    # pheno_reps = 1
    # n_effects = 10
    # h2 = 0.5
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
    list_Y_b_E_b_epi = simquantgen::fn_simulate_phenotypes(G=G, n_alleles=n_alleles, dist_effects="norm", n_effects=n_effects, h2=h2, pheno_reps=1, verbose=FALSE)
    y = list_Y_b_E_b_epi$Y[,1]
    df = data.frame(id=names(y), pop=sample(paste0("pop_", 1:n_pop), size=length(y), replace=TRUE), trait=y); rownames(df) = NULL
    ### Report distribution of simulation data
    if (verbose) {
        print("Simulated allele frequency distribution:")
        print(paste0(c("q_min=", "q_max="), round(range(G), 4)))
        txtplot::txtdensity(G[!is.na(G)])
        print("Simulated phenotype distribution:")
        print(paste0(c("y_min=", "y_max="), round(range(df$trait), 4)))
        txtplot::txtdensity(df$trait[!is.na(df$trait)])
    }
    ### Instantiate output file names
    fname_geno_vcf = NULL
    fname_geno_tsv = NULL
    fname_geno_rds = NULL
    fname_pheno_tsv = NULL
    ### Define date-time-random-number identifier so that we minimise the possibility of unintentional over-writing of the output file/s
    # date = gsub("-", "", gsub("[.]", "", gsub(":", "", gsub(" ", "", as.character(Sys.time())))))
    ### Save phenotype file
    if (save_pheno_tsv) {
        # fname_pheno_tsv = file.path(getwd(), paste0("simulated_phenotype-", date, ".tsv"))
        fname_pheno_tsv = tempfile(fileext=".tsv")
        utils::write.table(df, file=fname_pheno_tsv, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        if (verbose) {
            print("Output phenotype file (tsv):")
            print(fname_pheno_tsv)
        }
    }
    if (save_geno_vcf) {
        vcf = fn_G_to_vcf(G=G, min_depth=min_depth, max_depth=max_depth, verbose=verbose)
        if (methods::is(vcf, "gpError")) {
            error = chain(vcf, methods::new("gpError",
                code=231,
                message=paste0(
                    "Error in io::fn_simulate_data(...). ",
                    "Please set n_alleles=2 as we can only convert biallelic loci into VCF format at the moment."
                )))
            return(error)
        }
        # fname_geno_vcf = file.path(getwd(), paste0("simulated_genotype-", date, ".vcf.gz"))
        fname_geno_vcf = tempfile(fileext=".vcf.gz")
        vcfR::write.vcf(vcf, file=fname_geno_vcf)
        if (verbose) {
            print("Output genotype file (sync.gz):")
            print(fname_geno_vcf)
        }
    }
    if (save_geno_tsv)  {
        list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=G, verbose=verbose)
        if (methods::is(list_ids_chr_pos_all, "gpError")) {
            error = chain(list_ids_chr_pos_all, methods::new("gpError",
                code=232,
                message=paste0(
                    "Error in io::fn_simulate_data(...). ",
                    "Error type returned by list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=G, verbose=verbose)."
                )))
            return(error)
        }
        df_geno = data.frame(chr=list_ids_chr_pos_all$vec_chr, pos=list_ids_chr_pos_all$vec_pos, allele=list_ids_chr_pos_all$vec_all, t(G))
        # fname_geno_tsv = file.path(getwd(), paste0("simulated_genotype-", date, ".tsv"))
        fname_geno_tsv = tempfile(fileext=".tsv")
        utils::write.table(df_geno, file=fname_geno_tsv, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        if (verbose) {
            print("Output genotype file (tsv):")
            print(fname_geno_tsv)
        }
    }
    if (save_geno_rds) {
        ### Convert loci names so that the chromosome, position and allele are tab-delimited
        # fname_geno_rds = file.path(getwd(), paste0("simulated_genotype-", date, ".Rds"))
        fname_geno_rds = tempfile(fileext=".Rds")
        ### Revert to numeric if we have more than 52 alleles (limit of the Latin alphabet)
        if (!non_numeric_Rds) {
            saveRDS(G, file=fname_geno_rds)
        } else {
            G_non_numeric = fn_G_numeric_to_non_numeric(G=G, ploidy=ploidy, verbose=verbose)
            if (methods::is(G_non_numeric, "gpError")) {
                error = chain(G_non_numeric, methods::new("gpError",
                    code=233,
                    message=paste0(
                        "Error in io::fn_simulate_data(...). ",
                        "Error type returned by G_non_numeric = fn_G_numeric_to_non_numeric(G=G, ploidy=ploidy, verbose=verbose)."
                    )))
                return(error)
            }
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
#'  This may be in 1 of 3 formats:
#'      - VCF format with AD (allele depth) and/or GT (genotype), and DP (depth) fields,
#'          where multi-allelic loci are split into separate rows, 
#'      - allele frequency table saved as a tab-delimited file with a header line and the first 3 columns refer to the
#'          chromosome (chr), position (pos), and allele (allele),
#'          with subsequent columns referring to the allele frequencies of a sample, entry or pool.
#'          Names of the samples, entries, or pools in the header line can be any unique string of characters.
#'      - Rds file containing a single numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'          Row names can be any string of characters which identify the sample or entry or pool names.
#'          Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'          the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'          subsequent elements are optional which may refer to the allele identifier and other identifiers.
#' @param ploidy ploidy level which will generate genotype classes instead of continuous allele frequencies.
#'  If NULL, then continuous allele frequencies and no binning or classification into genotype classes
#'  will be performed (Default=NULL).
#' @param force_biallelic assume biallelic loci thereby dropping extra rows in vcf files corresponding to multi-allelic loci.
#'  Generate the vcf via `bcftools call -m ...` for multi-allelic and rare-variant calling,
#'  followed by `bcftools -m - ...` to  split mult-iallelic sites into biallelic records (-) (Default=TRUE)
#' @param retain_minus_one_alleles_per_locus omit the alternative or trailing allele per locus? (Default=TRUE)
#' @param min_depth if input is a VCF file: minimum depth per locus beyond which will be set to missing data (Default=0)
#' @param max_depth if input is a VCF file: maximum depth per locus beyond which will be set to missing data (Default=.Machine$integer.max)
#' @param verbose show genotype loading messages? (Default=FALSE)
#' @returns
#'  - Ok: numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(verbose=TRUE, 
#'  save_geno_vcf=TRUE, save_geno_tsv=TRUE, save_geno_rds=TRUE, save_pheno_tsv=TRUE)
#' G_vcf = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf, verbose=TRUE)
#' G_tsv = fn_load_genotype(fname_geno=list_sim$fname_geno_tsv, verbose=TRUE)
#' G_rds = fn_load_genotype(fname_geno=list_sim$fname_geno_rds, verbose=TRUE)
#' @export
fn_load_genotype = function(fname_geno, ploidy=NULL, force_biallelic=TRUE, retain_minus_one_alleles_per_locus=TRUE, 
    min_depth=0, max_depth=.Machine$integer.max, new_format=FALSE, verbose=FALSE) 
{
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(ploidy=8, verbose=TRUE, save_geno_vcf=TRUE, save_geno_tsv=TRUE, save_geno_rds=TRUE, save_pheno_tsv=TRUE)
    # # list_sim = fn_simulate_data(ploidy=8, verbose=TRUE, save_geno_vcf=TRUE, save_geno_tsv=TRUE, save_geno_rds=TRUE, save_pheno_tsv=TRUE, non_numeric_Rds=TRUE)
    # fname_geno = list_sim$fname_geno_vcf
    # fname_geno = list_sim$fname_geno_tsv
    # fname_geno = list_sim$fname_geno_rds
    # ploidy = 4
    # force_biallelic = TRUE
    # retain_minus_one_alleles_per_locus = TRUE
    # min_depth = 10
    # max_depth = 100
    # verbose = TRUE
    ###################################################
    if (verbose) {
        print("##########################")
        print("### Load genotype data ###")
        print("##########################")
    }
    ### Load the genotype matrix (n x p)
    ### TryCatch Rds, vcf, then tsv
    G = tryCatch({
        ################
        ### RDS file ###
        ################
        G = readRDS(fname_geno)
        if (is.numeric(G)==FALSE) {
            G = fn_G_non_numeric_to_numeric(G_non_numeric=G, verbose=verbose)
            if (methods::is(G, "gpError")) {
                error = chain(G, methods::new("gpError",
                    code=234,
                    message=paste0(
                        "Error in io::fn_load_genotype(...). ",
                        "Error type returned by G = fn_G_non_numeric_to_numeric(G_non_numeric=G, verbose=verbose)."
                    )))
                return(error)
            }
        }
        if (verbose) {print("Genotype loaded from an RDS file. No depth information available.")}
        G
    }, 
    error=function(e) {
        tryCatch({
            ################
            ### VCF file ###
            ################
            vcf = vcfR::read.vcfR(fname_geno, verbose=TRUE)
            list_G_D = fn_vcf_to_G(vcf, force_biallelic=force_biallelic, retain_minus_one_alleles_per_locus=retain_minus_one_alleles_per_locus, min_depth=min_depth, max_depth=max_depth, verbose=verbose)
            G = list_G_D$G
            if (methods::is(G, "gpError")) {
                error = chain(G, methods::new("gpError",
                    code=235,
                                        message=paste0(
                        "Error in io::fn_load_genotype(...).", 
                        "Error loading the vcf file: ", fname_geno, ".")))
                return(error)
            } else {
                if (verbose) {print("Genotype loaded from a VCF file. Depth information used.")}
                rm("vcf")
                rm("list_G_D")
                gc()
                return(G)
            }
        }, 
        error=function(e) {
            ########################################
            ### TSV: allele frequency table file ###
            ########################################
            df = utils::read.delim(fname_geno, sep="\t", header=TRUE, check.names=FALSE)
            if (!new_format) {
                if (!((grepl("chr", colnames(df)[1], ignore.case=TRUE)) &
                    (grepl("pos", colnames(df)[2], ignore.case=TRUE)) &
                    (grepl("allele", colnames(df)[3], ignore.case=TRUE)))
                ) {
                    error = methods::new("gpError", 
                        code=236,
                        message=paste0(
                            "Error in io::fn_load_genotype(...). ",
                            "The file: ", fname_geno, " is not in allele frequency table format as described in the README.md. ",
                            "The first 3 columns do not correspond to 'chr', 'pos', and 'allele'."))
                    return(error)
                }
                vec_loci_names = paste(df[,1], df[,2], df[,3], sep="\t")
                vec_entries = colnames(df)[c(-1:-3)]
                G = as.matrix(t(df[, c(-1:-3)]))
                rownames(G) = vec_entries
                colnames(G) = vec_loci_names
            } else {
                if (!((grepl("chr", colnames(df)[1], ignore.case=TRUE)) &
                    (grepl("pos", colnames(df)[2], ignore.case=TRUE)) &
                    (grepl("all_alleles", colnames(df)[3], ignore.case=TRUE)) &
                    (grepl("allele", colnames(df)[4], ignore.case=TRUE)))
                ) {
                    error = methods::new("gpError", 
                        code=236,
                        message=paste0(
                            "Error in io::fn_load_genotype(...). ",
                            "The file: ", fname_geno, " is not in allele frequency table format as described in the README.md. ",
                            "The first 3 columns do not correspond to 'chr', 'pos', and 'allele'."))
                    return(error)
                }
                vec_loci_names = paste(df[,1], df[,2], df[,3], df[,4], sep="\t")
                vec_entries = colnames(df)[c(-1:-4)]
                G = as.matrix(t(df[, c(-1:-4)]))
                rownames(G) = vec_entries
                colnames(G) = vec_loci_names
            }
            if (verbose) {print("Genotype loaded from a tab-delimited allele frequency table file. No depth information available.")}
            rm("df")
            gc()
            return(G)
        })
    })
    if (methods::is(G, "gpError")) {
        return(G)
    }
    ### Retain a-1 allele/s per locus, i.e. remove duplicates assuming all loci are biallelic
    if (retain_minus_one_alleles_per_locus==TRUE) {
        list_G_G_alt = fn_G_split_off_alternative_allele(G=G, verbose=verbose)
        if (methods::is(list_G_G_alt, "gpError")) {
            error = chain(list_G_G_alt, methods::new("gpError",
                code=237,
                message=paste0(
                    "Error in io::fn_load_genotype(...). ",
                    "Error type returned by list_G_G_alt = fn_G_split_off_alternative_allele(G=G, verbose=verbose)."
                )))
            return(error)
        }
        G = list_G_G_alt$G
        rm("list_G_G_alt")
        gc()
    }
    ### Bin allele frequencies
    if (!is.null(ploidy)) {
        G = fn_classify_allele_frequencies(G=G, ploidy=ploidy, verbose=verbose)
        if (methods::is(G, "gpError")) {
            error = chain(G, methods::new("gpError",
                code=238,
                message=paste0(
                    "Error in io::fn_load_genotype(...). ",
                    "Error type returned by G = fn_classify_allele_frequencies(G=G, ploidy=ploidy, verbose=verbose)"
                )))
            return(error)
        }
    }
    ### Show the allele frequency stats
    if (verbose) {
        print(paste0("Genotype data has ", nrow(G), " samples/entries/pools genotyped across ", ncol(G), " loci/alleles/SNPs."))
        print("Distribution of allele frequencies")
        vec_freqs_sample = sample(G, size=min(c(prod(dim(G)), 1e4)))
        vec_freqs_sample = vec_freqs_sample[!is.na(vec_freqs_sample)]
        txtplot::txtdensity(c(vec_freqs_sample, 1-vec_freqs_sample))
        print("Distribution of mean sparsity per locus")
        mat_sparsity = is.na(G)
        vec_sparsity_per_locus = colMeans(mat_sparsity, na.rm=TRUE)
        vec_sparsity_per_sample = rowMeans(mat_sparsity, na.rm=TRUE)
        txtplot::txtdensity(vec_sparsity_per_locus[!is.na(vec_sparsity_per_locus)])
        print("Distribution of mean sparsity per sample")
        txtplot::txtdensity(vec_sparsity_per_sample[!is.na(vec_sparsity_per_sample)])
        rm("mat_sparsity")
        rm("vec_sparsity_per_locus")
        rm("vec_sparsity_per_sample")
        gc()
    }
    gc()
    return(G)
}

#' Filter genotypes by:
#'  - minimum allele frequency,
#'  - minimum variation within locus (in an attempt to prevent duplicating the intercept)
#'  - loci identities using a SNP which specifies the SNP coordinates as well as the reference and alternative alleles
#'  - mean sparsity per locus
#'  - mean sparsity per sample/entry/pool
#'
#' @param G numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#' @param maf minimum allele frequency (Default=0.01)
#' @param sdev_min minimum allele frequency standard deviation (constant allele frequency across samples is just another intercept) (Default=0.0001)
#' @param fname_snp_list name of the file containing the list of expected SNPs including their coordinates and alleles.
#'  This is a tab-delimited file with 3 columns: '#CHROM', 'POS', 'REF,ALT', corresponding to (Default=NULL)
#'  chromosome names (e.g. 'chr1' & 'chrom_1'), 
#'  numeric positions (e.g. 12345 & 100001), and 
#'  reference-alternative allele strings separated by a comma (e.g. 'A,T' & 'C,G') (Default=NULL)
#' @param max_n_alleles maximum number of alleles per locus. Note that at max_n_alleles=1, we are assuming 
#' retain_minus_one_alleles_per_locus=TRUE in fn_load_genotype(...), and not that we want fixed - one allele per locus or site. (Default=NULL)
#' @param max_sparsity_per_locus maximum mean sparsity per locus, e.g. 0.1 or 0.5 (Default=NULL)
#' @param frac_topmost_sparse_loci_to_remove fraction of the top-most sparse loci to remove, e.g. 0.01 or 0.25 (Default=NULL)
#' @param n_topmost_sparse_loci_to_remove number of top-most sparse loci to remove, e.g. 100 or 1000 (Default=NULL)
#' @param max_sparsity_per_sample maximum mean sparsity per sample, e.g. 0.3 or 0.5 (Default=NULL)
#' @param frac_topmost_sparse_samples_to_remove fraction of the top-most sparse samples to remove, e.g. 0.01 or 0.05 (Default=NULL)
#' @param n_topmost_sparse_samples_to_remove number of top-most sparse samples to remove, e.g. 5 or 10 (Default=NULL)
#' @param verbose show genotype filtering messages? (Default=FALSE)
#' @returns
#'  - Ok: numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(verbose=TRUE)
#' G = fn_load_genotype(list_sim$fname_geno_vcf, retain_minus_one_alleles_per_locus=FALSE)
#' ### Rename allele_1 and allele_alt to A and T, respectively to allow filtering using a SNP list
#' colnames(G) = gsub("allele_1", "A", colnames(G))
#' colnames(G) = gsub("allele_alt", "T", colnames(G))
#' maf = 0.05
#' sdev_min = 0.0001
#' verbose = TRUE
#' ### Simulate SNP list for filtering
#' mat_loci = matrix(unlist(strsplit(colnames(G), "\t")), byrow=TRUE, ncol=3)
#' vec_loci = unique(paste0(mat_loci[,1], "\t", mat_loci[,2]))
#' mat_loci = matrix(unlist(strsplit(vec_loci, "\t")), byrow=TRUE, ncol=2)
#' df_snp_list = data.frame(CHROM=mat_loci[,1], 
#'  POS=as.numeric(mat_loci[,2]), 
#'  REF_ALT=paste0("A,T"))
#' df_snp_list$REF_ALT[1:100] = "C,G"
#' colnames(df_snp_list) = c("#CHROM", "POS", "REF,ALT")
#' fname_snp_list = tempfile(fileext=".snplist")
#' utils::write.table(df_snp_list, file=fname_snp_list, sep="\t", 
#'  row.names=FALSE, col.names=TRUE, quote=FALSE)
#' ### Filter
#' G_filtered = fn_filter_genotype(G=G, maf=0.05, fname_snp_list=fname_snp_list, verbose=TRUE)
#' @export
fn_filter_genotype = function(G, maf=0.01, sdev_min=0.0001, 
    fname_snp_list=NULL, max_n_alleles=NULL,
    max_sparsity_per_locus=NULL, frac_topmost_sparse_loci_to_remove=NULL, n_topmost_sparse_loci_to_remove=NULL, 
    max_sparsity_per_sample=NULL, frac_topmost_sparse_samples_to_remove=NULL, n_topmost_sparse_samples_to_remove=NULL, 
    verbose=FALSE)
{
    ###################################################
    ### TEST
    # # list_sim = fn_simulate_data(verbose=TRUE)
    # # G = fn_load_genotype(list_sim$fname_geno_vcf)
    # # maf = 0.05
    # # sdev_min = 0.0001
    # # fname_snp_list = NULL
    # ### Simulate SNP list for filtering
    # list_sim = fn_simulate_data(verbose=TRUE)
    # G = fn_load_genotype(list_sim$fname_geno_vcf, min_depth=42, max_depth=750, retain_minus_one_alleles_per_locus=FALSE)
    # ### Rename allele_1 and allele_alt to A and T, respectively to allow filtering using a SNP list
    # colnames(G) = gsub("allele_1", "A", colnames(G))
    # colnames(G) = gsub("allele_alt", "T", colnames(G))
    # maf = 0.05
    # sdev_min = 0.0001
    # max_n_alleles = 2
    # max_sparsity_per_locus = 0.4
    # frac_topmost_sparse_loci_to_remove = 0.01
    # n_topmost_sparse_loci_to_remove = 100
    # max_sparsity_per_sample = 0.3
    # frac_topmost_sparse_samples_to_remove = 0.01
    # n_topmost_sparse_samples_to_remove = 10
    # mat_loci = matrix(unlist(strsplit(colnames(G), "\t")), byrow=TRUE, ncol=3)
    # vec_loci = unique(paste0(mat_loci[,1], "\t", mat_loci[,2]))
    # mat_loci = matrix(unlist(strsplit(vec_loci, "\t")), byrow=TRUE, ncol=2)
    # df_snp_list = data.frame(CHROM=mat_loci[,1], POS=as.numeric(mat_loci[,2]), REF_ALT=paste0("allele_1,allele_alt"))
    # df_snp_list$REF_ALT[1:100] = "allele_2,allele_4"
    # colnames(df_snp_list) = c("#CHROM", "POS", "REF,ALT")
    # fname_snp_list = tempfile(fileext=".snplist")
    # utils::write.table(df_snp_list, file=fname_snp_list, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    # verbose = TRUE
    ###################################################
    if (verbose) {
        print("############################")
        print("### Filter genotype data ###")
        print("############################")
    }
    ### Input sanity check
    if (methods::is(G, "gpError")) {
        error = chain(G, methods::new("gpError",
            code=239,
            message=paste0(
                "Error in io::fn_filter_genotype(...). ",
                "Input G is an error type."
            )))
        return(error)
    }
    ### Make sure G contains allele frequencies
    if (sum(((G < 0) | (G > 1) | is.infinite(G)), na.rm=TRUE) > 0) {
        error = methods::new("gpError",
            code=240,
            message=paste0(
                "Error in io::fn_filter_genotype(...). ",
                "We are expecting a matrix allele frequencies but ",
                "we are getting negative values and/or values greater than 1 and/or infinite values."))
        return(error)
    }
    ### Make sure the input thresholds are sensible
    if ((maf < 0.0) | (maf > 1.0)) {
        error = methods::new("gpError",
            code=241,
            message=paste0(
                "Error in io::fn_filter_genotype(...). ",
                "Please use a minimum allele frequency (maf) between 0 and 1."
            ))
        return(error)
    }
    if ((sdev_min < 0.0) | (sdev_min > 1.0)) {
        error = methods::new("gpError",
            code=242,
            message=paste0(
                "Error in io::fn_filter_genotype(...). ",
                "Please use a minimum standard deviation in allele frequency (sdev_min) between 0 and 1."
            ))
        return(error)
    }
    if (!is.null(max_n_alleles)) {
        if (max_n_alleles < 1) {
            error = methods::new("gpError",
                code=243,
                message=paste0(
                    "Error in io::fn_filter_genotype(...). ",
                    "Please use a maximum number of alleles per locus (max_n_alleles) of at least 1. ",
                    "Note that at max_n_alleles=1, we are assuming retain_minus_one_alleles_per_locus=TRUE in fn_load_genotype(...), ",
                    "and not that we want fixed - one allele per locus sites."
                ))
            return(error)
        }
    }
    if (!is.null(max_sparsity_per_locus)) {
        if ((max_sparsity_per_locus < 0.0) | (max_sparsity_per_locus > 1.0)) {
            error = methods::new("gpError",
                code=244,
                message=paste0(
                    "Error in io::fn_filter_genotype(...). ",
                    "Please use a maximum sparsity per locus (max_sparsity_per_locus) between 0 and 1."
                ))
            return(error)
        }
    }
    if (!is.null(frac_topmost_sparse_loci_to_remove)) {
        if ((frac_topmost_sparse_loci_to_remove < 0.0) | (frac_topmost_sparse_loci_to_remove > 1.0)) {
            error = methods::new("gpError",
                code=245,
                message=paste0(
                    "Error in io::fn_filter_genotype(...). ",
                    "Please use a fraction of the top-most sparse loci (frac_topmost_sparse_loci_to_remove) between 0 and 1."
                ))
            return(error)
        }
    }
    if (!is.null(n_topmost_sparse_loci_to_remove)) {
        if ((n_topmost_sparse_loci_to_remove < 0) | (n_topmost_sparse_loci_to_remove > ncol(G))) {
            error = methods::new("gpError",
                code=246,
                message=paste0(
                    "Error in io::fn_filter_genotype(...). ",
                    "Please use a number of top-most sparse loci to remove (n_topmost_sparse_loci_to_remove) between 0 and ", ncol(G), "."
                ))
            return(error)
        }
    }
    if (!is.null(max_sparsity_per_sample)) {
        if ((max_sparsity_per_sample < 0.0) | (max_sparsity_per_sample > 1.0)) {
            error = methods::new("gpError",
                code=247,
                message=paste0(
                    "Error in io::fn_filter_genotype(...). ",
                    "Please use a maximum sparsity per sample (max_sparsity_per_sample) between 0 and 1."
                ))
            return(error)
        }
    }
    if (!is.null(frac_topmost_sparse_samples_to_remove)) {
        if ((frac_topmost_sparse_samples_to_remove < 0.0) | (frac_topmost_sparse_samples_to_remove > 1.0)) {
            error = methods::new("gpError",
                code=248,
                message=paste0(
                    "Error in io::fn_filter_samples(...). ",
                    "Please use a fraction of the top-most sparse samples (frac_topmost_sparse_samples_to_remove) between 0 and 1."
                ))
            return(error)
        }
    }
    if (!is.null(n_topmost_sparse_samples_to_remove)) {
        if ((n_topmost_sparse_samples_to_remove < 0) | (n_topmost_sparse_samples_to_remove > nrow(G))) {
            error = methods::new("gpError",
                code=249,
                message=paste0(
                    "Error in io::fn_filter_samples(...). ",
                    "Please use a number of top-most sparse samples to remove (n_topmost_sparse_samples_to_remove) between 0 and ", nrow(G), "."
                ))
            return(error)
        }
    }
    ### Extract the mean and standard deviation in allele frequencies per locus
    vec_freqs = colMeans(G, na.rm=TRUE)
    vec_sdevs = apply(G, MARGIN=2, FUN=stats::sd, na.rm=TRUE)
    if (verbose) {
        print("Distribution of mean allele frequencies per locus:")
        txtplot::txtdensity(vec_freqs[!is.na(vec_freqs)])
        print("Distribution of allele frequency standard deviation per locus:")
        txtplot::txtdensity(vec_sdevs[!is.na(vec_sdevs)])
    }
    ### Filter by minimum allele frequency and minimum allele frequency variance (constant allele frequency across samples is just another intercept)
    vec_idx = which(
        (vec_freqs >= maf) &
        (vec_freqs <= (1-maf)) &
        (vec_sdevs >= sdev_min))
    rm("vec_freqs")
    rm("vec_sdevs")
    gc()
    if (length(vec_idx) == 0) {
        error = methods::new("gpError",
            code=250,
            message=paste0(
                "Error in io::fn_filter_genotype(...). ",
                "All loci did not pass the minimum allele frequency (", maf, ") and minimum allele frequency standard deviation (", sdev_min, ")."
            ))
        return(error)
    } else if (length(vec_idx) < ncol(G)) {
        if (verbose) {
            print(paste0("Filtering by minimum allele frequency (", maf, ") and allele frequency standard deviation (", sdev_min, "):"))
            print(paste0("Retaining ", length(vec_idx), " loci, i.e. ", length(vec_idx), "/", ncol(G), " (", round(length(vec_idx)*100/ncol(G)), "% retained)"))
        }
        G = G[, vec_idx, drop=FALSE]
    } else {
        if (verbose) {print("All loci passed the minimum allele frequency and standard deviation thresholds.")}
    }
    ### Filter using a SNP list
    if (!is.null(fname_snp_list)) {
        list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=G, verbose=verbose)
        if (methods::is(list_ids_chr_pos_all, "gpError")) {
            error = chain(list_ids_chr_pos_all, methods::new("gpError",
                code=251,
                message=paste0(
                    "Error in io::fn_filter_genotype(...). ",
                    "Error in list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=G, verbose=verbose)."
                )))
            return(error)
        }
        if (length(unique(list_ids_chr_pos_all$vec_all)) < 1) {
            error = methods::new("gpError",
                code=252,
                message=paste0(
                    "Error in io::fn_filter_genotype(...). ",
                    "Please make sure all loci have an associated allele, if you wish to filter using a SNP list."
                ))
            return(error)
        }
        vec_observed_snps = paste(
            list_ids_chr_pos_all$vec_chr,
            list_ids_chr_pos_all$vec_pos,
            list_ids_chr_pos_all$vec_all, sep="\t")
        df = utils::read.delim(fname_snp_list, sep="\t", header=TRUE, check.names=FALSE)
        if (sum(colnames(df) != c("#CHROM", "POS", "REF,ALT")) > 0) {
            error = methods::new("gpError",
                code=253,
                message=paste0(
                    "Error in io::fn_filter_genotype(...). ",
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
        if (verbose) {print("Filtering using a SNP list.")}
        vec_idx = which(vec_observed_snps %in% vec_expected_snps)
        if (length(vec_idx) == 0) {
            error = methods::new("gpError",
                code=254,
                message=paste0(
                    "Error in io::fn_filter_genotype(...). ",
                    "None of the loci are in the SNP list."
                ))
            return(error)
        } else if (length(vec_idx) < ncol(G)) {
            if (verbose) {print(paste0("Filtered out ", ncol(G)-length(vec_idx), " loci not included in the SNP list: ", fname_snp_list, "."))}
            G = G[, vec_idx, drop=FALSE]
        } else {
            if (verbose) {print("All loci were found in the SNP list.")}
        }
    }
    ### Filter by maximum number of alleles per locus
    if (!is.null(max_n_alleles)) {
        list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=G, verbose=verbose)
        if (methods::is(list_ids_chr_pos_all, "gpError")) {
            error = chain(list_ids_chr_pos_all, methods::new("gpError",
                code=255,
                message=paste0(
                    "Error in io::fn_filter_genotype(...). ",
                    "Error in list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=G, verbose=verbose)."
                )))
            return(error)
        }
        vec_loci_names = paste0(list_ids_chr_pos_all$vec_chr, "\t", list_ids_chr_pos_all$vec_pos)
        vec_allele_counts = table(vec_loci_names)
        if (verbose) {
            print("Distribution of the number of alleles per locus:")
            txtplot::txtdensity(vec_allele_counts[!is.na(vec_allele_counts)])
        }
        vec_loci_names_passed = names(vec_allele_counts)[vec_allele_counts <= max_n_alleles]
        vec_idx = which(vec_loci_names %in% vec_loci_names_passed)
        if (length(vec_idx) == 0) {
            error = methods::new("gpError",
                code=256,
                message=paste0(
                    "Error in io::fn_filter_genotype(...). ",
                    "All loci were filtered out. ",
                    "All loci have more than ", max_n_alleles, " alleles per locus. ",
                    "Are you including all alleles per locus and using max_n_alleles=1? ",
                    "If so, then please use retain_minus_one_alleles_per_locus=TRUE in fn_load_genotype(...)."
                ))
            return(error)
        } else if (length(vec_idx) < ncol(G)) {
            if (verbose) {print(paste0("Removing ", ncol(G)-length(vec_idx), " loci which have more than ", max_n_alleles, " allele/s per locus."))}
            G = G[, vec_idx, drop=FALSE]
        } else {
            if (verbose) {print(paste0("All loci have a maximum number of alleles per locus of ", max_n_alleles, "."))}
        }
        rm("list_ids_chr_pos_all")
        gc()
    }
    ### Filtering by sparsity below only works if the input genotype data have NA
    ### Define the sparsity matrix
    mat_sparsity = is.na(G)
    ### Filter by mean sparsity per locus
    vec_sparsity_per_locus = colMeans(mat_sparsity, na.rm=TRUE)
    if (verbose) {
        print("Distribution of mean sparsity per locus")
        txtplot::txtdensity(vec_sparsity_per_locus[!is.na(vec_sparsity_per_locus)])
    }
    vec_idx_loci_to_remove = c()
    if (!is.null(max_sparsity_per_locus)) {
        if (verbose) {print(paste0("Filtering by maximum mean sparsity per locus of ", max_sparsity_per_locus))}
        vec_idx_loci_to_remove = which(vec_sparsity_per_locus > max_sparsity_per_locus)
        if (verbose & (length(vec_idx_loci_to_remove)==ncol(G))) {
            error = methods::new("gpError",
                code=257,
                message=paste0("All loci were filtered out. Please consider increasing the max_sparsity_per_locus from ", max_sparsity_per_locus, 
                ", given that the mean sparsity per locus ranges from ", min(vec_sparsity_per_locus, na.rm=TRUE), 
                " to ", max(vec_sparsity_per_locus, na.rm=TRUE), 
                " with a mean of ", mean(vec_sparsity_per_locus, na.rm=TRUE), 
                " and a median of ", stats::median(vec_sparsity_per_locus, na.rm=TRUE)))
            return(error)
        }
    }
    if (!is.null(frac_topmost_sparse_loci_to_remove)) {
        if (verbose) {print(paste0("Filtering out the top ", round(frac_topmost_sparse_loci_to_remove*100), "% most sparse loci."))}
        n_loci_to_remove = round(frac_topmost_sparse_loci_to_remove * ncol(G))
        if (n_loci_to_remove > length(vec_idx_loci_to_remove)) {
            vec_idx_sort_decreasing_sparsity = order(vec_sparsity_per_locus, decreasing=TRUE)
            vec_idx_loci_to_remove = sort(unique(c(vec_idx_loci_to_remove, vec_idx_sort_decreasing_sparsity[1:n_loci_to_remove])))
        }
        if (verbose & (length(vec_idx_loci_to_remove)==ncol(G))) {
            error = methods::new("gpError",
                code=258,
                message=paste0("All loci were filtered out. Please consider decreasing the frac_topmost_sparse_loci_to_remove from ", frac_topmost_sparse_loci_to_remove, 
                " to something more reasonable."))
            return(error)
        }
    }
    if (!is.null(n_topmost_sparse_loci_to_remove)) {
        if (verbose) {print(paste0("Filtering out the top ", n_topmost_sparse_loci_to_remove, " most sparse loci."))}
        if (n_topmost_sparse_loci_to_remove > length(vec_idx_loci_to_remove)) {
            vec_idx_sort_decreasing_sparsity = order(vec_sparsity_per_locus, decreasing=TRUE)
            vec_idx_loci_to_remove = sort(unique(c(vec_idx_loci_to_remove, vec_idx_sort_decreasing_sparsity[1:n_topmost_sparse_loci_to_remove])))
        }
        if (verbose & (length(vec_idx_loci_to_remove)==ncol(G))) {
            error = methods::new("gpError",
                code=259,
                message=paste0("All loci were filtered out. Please consider decreasing the n_topmost_sparse_loci_to_remove from ", n_topmost_sparse_loci_to_remove, 
                " to something more reasonable, if it please you m'lady/m'lord."))
            return(error)
        }
    }
    if (length(vec_idx_loci_to_remove) > 0) {
        if (verbose) {print(paste0("Filtered out ", length(vec_idx_loci_to_remove), " most sparse loci."))}
        vec_idx = which(!(c(1:ncol(G)) %in% vec_idx_loci_to_remove))
        G = G[, vec_idx]
    } else {
        if (verbose) {print("All loci passed the filtering by mean sparsity per locus.")}
    }
    ### Filter by mean sparsity per sample
    vec_sparsity_per_sample = rowMeans(mat_sparsity, na.rm=TRUE)
    if (verbose) {
        print("Distribution of mean sparsity per sample")
        txtplot::txtdensity(vec_sparsity_per_sample[!is.na(vec_sparsity_per_sample)])
    }
    vec_idx_samples_to_remove = c()
    if (!is.null(max_sparsity_per_sample)) {
        if (verbose) {print(paste0("Filtering by maximum mean sparsity per sample of ", max_sparsity_per_sample))}
        vec_idx_samples_to_remove = which(vec_sparsity_per_sample > max_sparsity_per_sample)
        if (verbose & (length(vec_idx_samples_to_remove)==nrow(G))) {
            error = methods::new("gpError", 
                code=260,
                message=paste0("All samples were filtered out. Please consider increasing the max_sparsity_per_sample from ", max_sparsity_per_sample, 
                ", given that the mean sparsity per sample ranges from ", min(vec_sparsity_per_sample, na.rm=TRUE), 
                " to ", max(vec_sparsity_per_sample, na.rm=TRUE), 
                " with a mean of ", mean(vec_sparsity_per_sample, na.rm=TRUE), 
                " and a median of ", stats::median(vec_sparsity_per_sample, na.rm=TRUE)))
            return(error)
        }
    }
    if (!is.null(frac_topmost_sparse_samples_to_remove)) {
        if (verbose) {print(paste0("Filtering out the top ", round(frac_topmost_sparse_samples_to_remove*100), "% most sparse samples."))}
        n_samples_to_remove = round(frac_topmost_sparse_samples_to_remove * nrow(G))
        if (n_samples_to_remove > length(vec_idx_samples_to_remove)) {
            vec_idx_sort_decreasing_sparsity = order(vec_sparsity_per_sample, decreasing=TRUE)
            vec_idx_samples_to_remove = sort(unique(c(vec_idx_samples_to_remove, vec_idx_sort_decreasing_sparsity[1:n_samples_to_remove])))
        }
        if (verbose & (length(vec_idx_samples_to_remove)==nrow(G))) {
            error = methods::new("gpError", 
                code=261,
                message=paste0("All samples were filtered out. Please consider decreasing the frac_topmost_sparse_samples_to_remove from ", frac_topmost_sparse_samples_to_remove, 
                " to something more reasonable."))
            return(error)
        }
    }
    if (!is.null(n_topmost_sparse_samples_to_remove)) {
        if (verbose) {print(paste0("Filtering out the top ", n_topmost_sparse_samples_to_remove, " most sparse samples."))}
        if (n_topmost_sparse_samples_to_remove > length(vec_idx_samples_to_remove)) {
            vec_idx_sort_decreasing_sparsity = order(vec_sparsity_per_sample, decreasing=TRUE)
            vec_idx_samples_to_remove = sort(unique(c(vec_idx_samples_to_remove, vec_idx_sort_decreasing_sparsity[1:n_topmost_sparse_samples_to_remove])))
        }
        if (verbose & (length(vec_idx_samples_to_remove)==nrow(G))) {
            error = methods::new("gpError", 
                code=262,
                message=paste0("All samples were filtered out. Please consider decreasing the n_topmost_sparse_samples_to_remove from ", n_topmost_sparse_samples_to_remove, 
                " to something more reasonable."))
            return(error)
        }
    }
    if (length(vec_idx_samples_to_remove) > 0) {
        if (verbose) {print(paste0("Filtered out ", length(vec_idx_samples_to_remove), " most sparse samples."))}
        vec_idx = which(!(c(1:nrow(G)) %in% vec_idx_samples_to_remove))
        G = G[vec_idx, ]
    } else {
        if (verbose) {print("All samples passed the filtering by mean sparsity per sample.")}
    }
    ### We need to repeat filtering by minimum allele frequency and minimum allele frequency variance 
    ### because the mean allele frequencies would have changed significantly after the multiple filtering steps above
    vec_freqs = colMeans(G, na.rm=TRUE)
    vec_sdevs = apply(G, MARGIN=2, FUN=stats::sd, na.rm=TRUE)
    vec_idx = which(
        (vec_freqs >= maf) &
        (vec_freqs <= (1-maf)) &
        (vec_sdevs >= sdev_min))
    if (length(vec_idx) == 0) {
        error = methods::new("gpError",
            code=263,
            message=paste0(
                "Error in io::fn_filter_genotype(...). ",
                "All loci did not pass the minimum allele frequency (", maf, ") and minimum allele frequency standard deviation (", sdev_min, ")."
            ))
        return(error)
    } else if (length(vec_idx) < ncol(G)) {
        if (verbose) {
            print(paste0("[Repeat] Filtering by minimum allele frequency (", maf, ") and allele frequency standard deviation (", sdev_min, "):"))
            print(paste0("[Repeat] Retaining ", length(vec_idx), " loci, i.e. ", length(vec_idx), "/", ncol(G), " (", round(length(vec_idx)*100/ncol(G)), "% retained)"))
        }
        G = G[, vec_idx, drop=FALSE]
    } else {
        if (verbose) {print("[Repeat] All loci passed the minimum allele frequency and standard deviation thresholds.")}
    }
    if (verbose) {
        print(paste0("After genotype filtering we retain n=", nrow(G), " with p=", ncol(G)))
        mat_idx_missing = is.na(G)
        print(paste0("Mean sparsity = ", round(100*mean(mat_idx_missing), 2), "%"))
        print("Allele frequency distribution: ")
        vec_freqs_sample = sample(unlist(G[!mat_idx_missing]), size=min(c(sum(!mat_idx_missing), 1e4)))
        txtplot::txtdensity(c(vec_freqs_sample, 1-vec_freqs_sample))
        vec_freqs_per_locus = colMeans(G, na.rm=TRUE)
        vec_freqs_per_locus = vec_freqs_per_locus[!is.na(vec_freqs_per_locus)]
        print(paste0("Mean allele frequencies per locus range from ", min(vec_freqs_per_locus), " to ", max(vec_freqs_per_locus)))
    }
    ### Clean-up
    gc()
    ### Return filtered allele frequency matrix
    return(G)
}

#' Merge two genotypes matrices where if there are conflicts:
#'  - data on the first matrix will be used, or
#'  - data on the second matrix will be used, or
#'  - arithmetic mean between the two matrices will be used.
#'
#' @param G1 numeric n1 samples x p1 loci-alleles matrix of allele frequencies with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#' @param G2 numeric n2 samples x p2 loci-alleles matrix of allele frequencies with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#' @param str_conflict_resolution conflict resolution mode. Use "1" to always use the genotype data from the 
#'  first matrix, "2" to to always use the data from the second matrix, and "3" to compute the arithmetic mean
#'  between the two matrices.
#' @param keep_common_loci_only keep only the common loci across the 2 genotype matrices? (Default=FALSE)
#' @param verbose show genotype merging messages? (Default=FALSE)
#' @returns
#'  - Ok: numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'  - Err: gpError
#' @examples
#' list_sim = gp::fn_simulate_data(verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' G1 = G[1:ceiling(0.5*nrow(G)), ]
#' G2 = G[(ceiling(0.5*nrow(G))+1):nrow(G), ]
#' G_merged = fn_merge_genotype_genotype(G1=G1, G2=G2, verbose=TRUE)
#' @export
fn_merge_genotype_genotype = function(G1, G2, str_conflict_resolution=c("1-use-G1", "2-use-G2", "3-use_mean")[3], keep_common_loci_only=FALSE, verbose=FALSE) {
    ###################################################
    ### TEST
	# list_sim = gp::fn_simulate_data(verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
	# G1 = G[sample.int(nrow(G), size=50), sample.int(ncol(G), size=500), drop=FALSE]
	# G2 = G[sample.int(nrow(G), size=50), sample.int(ncol(G), size=500), drop=FALSE]
	# str_conflict_resolution = c("1-use-G1", "2-use-G2", "3-use_mean")[3]
    # keep_common_loci_only = TRUE
    # verbose = TRUE
    ###################################################
    if (verbose) {
        print("###########################")
        print("### Merge genotype data ###")
        print("###########################")
    }
    ### Input sanity check
    if (methods::is(G1, "gpError") | methods::is(G2, "gpError")) {
        if (methods::is(G1, "gpError")) {
            error = chain(G1, methods::new("gpError",
                code=281,
                message=paste0(
                    "Error in io::fn_merge_genotype_genotype(...). ",
                    "Input G1 is an error type."
                )))
        } else if (methods::is(G2, "gpError")) {
            error = chain(G2, methods::new("gpError",
                code=282,
                message=paste0(
                    "Error in io::fn_merge_genotype_genotype(...). ",
                    "Input G2 is an error type."
                )))
        } else {
            error = chain(G1, chain(G2, methods::new("gpError",
                code=283,
                message=paste0(
                    "Error in io::fn_merge_genotype_genotype(...). ",
                    "Inputs G1 and G2 are error types."
                ))))
        }
        return(error)
    }
    ### Do we want to keep only the common loci
    if (keep_common_loci_only) {
        vec_common_loci_names = intersect(colnames(G1), colnames(G2))
        if (length(vec_common_loci_names) == 0) {
            error = methods::new("gpError",
                code=281,
                message=paste0(
                    "Error in io::fn_merge_genotype_genotype(...). ",
                    "No common loci across the two genotype matrices."
                )
            )
            return(error)
        }
        G1 = G1[, (colnames(G1) %in% vec_common_loci_names), drop=FALSE]
        G2 = G2[, (colnames(G2) %in% vec_common_loci_names), drop=FALSE]
    }
    ### Extract row and column names, i.e. sample and loci names
	vec_G1_row_names = rownames(G1); vec_G1_column_names = colnames(G1)
	vec_G2_row_names = rownames(G2); vec_G2_column_names = colnames(G2)
	### Merge the 2 matrices where G1 takes precedence over G2, where we simply add the unique columns in G2.
    ### This means that in the merged matrix, data are missing at common loci in the G2 samples.
	vec_G2_bool_unique_loci = !(vec_G2_column_names %in% vec_G1_column_names)
    df_G_merged = merge(
		data.frame(ID=vec_G1_row_names, G1, check.names=FALSE),
		data.frame(ID=vec_G2_row_names, G2[, vec_G2_bool_unique_loci, drop=FALSE], check.names=FALSE),
		by="ID", all=TRUE)
	### Convert the merged genotype data frames into a matrix
	G_merged = as.matrix(df_G_merged[, -1, drop=FALSE]); rownames(G_merged) = df_G_merged$ID
	vec_G_merged_row_names = rownames(G_merged); vec_G_merged_column_names = colnames(G_merged)
	### Define the intersections
	vec_common_row_names = intersect(vec_G1_row_names, vec_G2_row_names)
	vec_common_column_names = intersect(vec_G1_column_names, vec_G2_column_names)
	### Insert G2 data into the intersecting columns (loci)
	if (sum(!vec_G2_bool_unique_loci) > 0) {
		pb = utils::txtProgressBar(min=0, max=nrow(G2), style=3)
		for (i in 1:nrow(G2)) {
			# i = 1
			row_name = vec_G2_row_names[i]
			idx_G_merged_row = which(vec_G_merged_row_names == row_name)
			idx_G2_row = which(vec_G2_row_names == row_name)
			vec_G_merged_idx_column_sort = which(vec_G_merged_column_names %in% vec_G2_column_names)
			vec_G_merged_idx_column_sort = vec_G_merged_idx_column_sort[order(vec_G_merged_column_names[vec_G_merged_idx_column_sort])]
			vec_G2_idx_column_sort = which(vec_G2_column_names %in% vec_G2_column_names)
			vec_G2_idx_column_sort = vec_G2_idx_column_sort[order(vec_G2_column_names[vec_G2_idx_column_sort])]
			if (sum(vec_G_merged_column_names[vec_G_merged_idx_column_sort] == vec_G2_column_names[vec_G2_idx_column_sort]) != length(vec_G_merged_idx_column_sort)) {
				error = methods::new("dbError",
					code=000,
					message=paste0("Error in fn_merge_genotype_genotype(...): the merged sample names do not match with G2 sample names."))
				return(error)
			}
			G_merged[idx_G_merged_row, vec_G_merged_idx_column_sort] = G2[idx_G2_row, vec_G2_idx_column_sort]
			utils::setTxtProgressBar(pb, i)
		}
		close(pb)
	}
	### Fix conflicts on common rows and columns
	if ((length(vec_common_row_names) > 0) && (length(vec_common_column_names) > 0)) {
		pb = utils::txtProgressBar(min=0, max=(length(vec_common_row_names)*length(vec_common_column_names)), style=3); counter = 1
		for (row_name in vec_common_row_names) {
			for (column_name in vec_common_column_names) {
				# row_name = vec_common_row_names[1]; column_name = vec_common_column_names[1]
				idx_G1_row = which(vec_G1_row_names == row_name)
				idx_G1_column = which(vec_G1_column_names == column_name)
				idx_G2_row = which(vec_G2_row_names == row_name)
				idx_G2_column = which(vec_G2_column_names == column_name)
				g1 = G1[idx_G1_row, idx_G1_column]
				g2 = G2[idx_G2_row, idx_G2_column]
				if (grepl("^1", str_conflict_resolution)) {
					g = g1
				} else if (grepl("^2", str_conflict_resolution)) {
					g = g2
				} else {
					if (is.na(g1) & is.na(g2)) {
						g = NA
					} else {
						g = mean(c(g1, g2), na.rm=TRUE)
					}
				}
				G_merged[vec_G_merged_row_names == row_name, vec_G_merged_column_names == column_name] = g
				utils::setTxtProgressBar(pb, counter); counter = counter + 1
			}
		}
		close(pb)
	} else if (verbose) {
		if ((length(vec_common_row_names) == 0) && (length(vec_common_column_names) == 0)) {
			print("There are no common samples and loci between the 2 genotype matrices.")
		} else if (length(vec_common_row_names) == 0) {
			print("There are no common samples between the 2 genotype matrices.")
		} else {
			print("There are no common loci between the 2 genotype matrices.")
		}
	}
	### Output
	return(G_merged)
}

#' Save numeric genotype matrix into a file
#'
#' @param G numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#' @param fname file name of the output file
#' @param file_type output file format (Default="RDS"). Choose from:
#'  - "RDS": Rds file containing a single numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'      Row names can be any string of characters which identify the sample or entry or pool names.
#'      Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'      the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'      subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'  - "TSV": allele frequency table saved as a tab-delimited file with a header line and the first 3 columns refer to the
#'      chromosome (chr), position (pos), and allele (allele),
#'      with subsequent columns referring to the allele frequencies of a sample, entry or pool.
#'      Names of the samples, entries, or pools in the header line can be any unique string of characters.
#' @param verbose show genotype saving messages? (Default=FALSE)
#' @returns
#'  - Ok: file name of the output file
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(verbose=TRUE)
#' G = fn_load_genotype(list_sim$fname_geno_vcf, min_depth=42, max_depth=750)
#' fn_save_genotype(G=G, fname=tempfile(fileext=".Rds"), file_type="Rds", verbose=TRUE)
#' fn_save_genotype(G=G, fname=tempfile(fileext=".tsv"), file_type="tsv", verbose=TRUE)
#' @export
fn_save_genotype = function(G, fname, file_type=c("RDS", "TSV")[1], verbose=FALSE) {
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(verbose=TRUE)
    # G = fn_load_genotype(list_sim$fname_geno_vcf, min_depth=42, max_depth=750, retain_minus_one_alleles_per_locus=FALSE)
    # fname = "test.tsv"
    # file_type=c("RDS", "TSV")[2]
    # verbose = TRUE
    ###################################################### Input sanity check
    if (methods::is(G, "gpError")) {
        error = chain(G, methods::new("gpError",
            code=264,
            message=paste0(
                "Error in io::fn_save_genotype(...). ",
                "Input G is an error type."
            )))
        return(error)
    }
    ### Make sure G contains allele frequencies
    if (sum(((G < 0) | (G > 1) | is.infinite(G)), na.rm=TRUE) > 0) {
        error = methods::new("gpError",
            code=265,
            message=paste0(
                "Error in io::fn_save_genotype(...). ",
                "We are expecting a matrix allele frequencies but ",
                "we are getting negative values and/or values greater than 1 and/or infinite values."))
        return(error)
    }
    ### Make sure the SNP names are tab-delimited by extracting the names
    list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=G, verbose=verbose)
    if (methods::is(list_ids_chr_pos_all, "gpError")) {
        error = chain(list_ids_chr_pos_all, methods::new("gpError",
            code=266,
            message=paste0(
                "Error in io::fn_save_genotype(...). ",
                "Error in list_ids_chr_pos_all = fn_G_extract_names(mat_genotypes=G, verbose=verbose)."
            )))
        return(error)
    }
    ### Save
    if (grepl("RDS", file_type, ignore.case=TRUE)) {
        if (verbose) {print(paste0("Saving the genotype matrix as an R object (RDS): ", fname))}
        saveRDS(G, file=fname)
    } else if (grepl("TSV", file_type, ignore.case=TRUE)) {
        if (verbose) {print(paste0("Saving the genotype matrix as tab-delimited allele frequency table (TSV): ", fname))}
        df_allele_freq = data.frame(chr=list_ids_chr_pos_all$vec_chr, pos=list_ids_chr_pos_all$vec_pos, allele=list_ids_chr_pos_all$vec_all, t(G))
        colnames(df_allele_freq) = c("chr", "pos", "allele", rownames(G))
        utils::write.table(df_allele_freq, file=fname, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    } else {
        error = methods::new("gpError",
            code=267,
            message=paste0(
                "Error in io::fn_save_genotype(...). ",
                "Unrecognised file type, please use 'RDS' or 'TSV'"
            ))
        return(error)
    }
    ### Output the filename
    return(fname)
}

#' Load phenotype data from a text-delimited file
#'
#' @param fname_pheno filename of the text-delimited phenotype data
#' @param sep column-delimiter in the phenotype file (Default="\\t")
#' @param header does the phenotype file have a header line? (Default=TRUE)
#' @param idx_col_id which column correspond to the sample/entry/pool/genotype names? (Default=1)
#' @param idx_col_pop which column correspond to the population groupings? (Default=2)
#' @param idx_col_y which column correspond to the phenotype data? (Default=3)
#' @param na_strings vector of string corresponding to missing data in the phenotype column (Default=c("", "-", "NA", "na", "NaN", "missing", "MISSING"))
#' @param verbose show phenotype loading messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'      + $y: named vector of numeric phenotype data
#'      + $pop: population or groupings corresponding to each element of y
#'      + $trait_name: name of the trait
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' @export
fn_load_phenotype = function(fname_pheno, sep="\t", header=TRUE, 
    idx_col_id=1, idx_col_pop=2, idx_col_y=3, 
    na_strings=c("", "-", "NA", "na", "NaN", "missing", "MISSING"), verbose=FALSE)
{
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    # fname_pheno = list_sim$fname_pheno_tsv
    # sep="\t"
    # header=TRUE
    # idx_col_id=1
    # idx_col_pop=2
    # idx_col_y=3
    # na_strings=c("", "-", "NA", "na", "NaN", "missing", "MISSING")
    # verbose = TRUE
    ###################################################
    if (verbose) {
        print("###########################")
        print("### Load phenotype data ###")
        print("###########################")
    }
    df = utils::read.table(fname_pheno, sep=sep, header=header, na.strings=na_strings, check.names=FALSE)
    ### Convert special symbols into underscores
    vec_column_names = colnames(df)
    vec_symbols = c(" ", "\t", "`", "~", "!", "@", "#", "\\$", "%", "\\^", "&", "\\*", "\\(", "\\)", "=", "\\+", "\\[", "\\]", "\\{", "\\}", "<", ">", "\\\\", "\\|", "'", '"', ";", "/", "\\?", "__", "__")
    for (symbol in vec_symbols) {
        vec_column_names = stringr::str_replace_all(string=vec_column_names, pattern=symbol, replacement="_")
    }
    colnames(df) = vec_column_names
    if (max(c(idx_col_y, idx_col_id, idx_col_pop)) > ncol(df)) {
        error = methods::new("gpError",
            code=268,
            message=paste0(
                "Error in io::fn_load_phenotype(...). ",
                "The requested columns: ", idx_col_id, ", ", idx_col_pop, " and ", idx_col_y, " are incompatible with ",
                "the dimensions of the loaded phenotype file: ", fname_pheno,
                ", which has ", ncol(df), " columns and ", nrow(df), " rows. ",
                "Are you certain that your file is separated by: '", sep, "'? ",
                "Are the sample IDs really at column '", idx_col_id, "'? ",
                "Are the population IDs really at column '", idx_col_pop, "'? ",
                "Are the phenotype data really at column '", idx_col_y, "'? ",
                "Are your missing data encoded as any of these: ", paste(paste0("'", na_strings, "'"), collapse=", "), "?"))
        return(error)
    }
    entry = as.character(df[, idx_col_id])
    pop = as.character(df[, idx_col_pop])
    y = df[, idx_col_y]
    if (sum(duplicated(entry)) > 0) {
        error = methods::new("gpError",
            code=281,
            message=paste0(
                "Error in io::fn_load_phenotype(...). ",
                "The sample names have duplicates. ",
                "We expect unique samples in the phenotype file. ",
                "Please remove duplicated or extract BLUEs/BLUPs using an appropriate linear model."))
        return(error)
    }
    if (!is.numeric(y)) {
        error = methods::new("gpError",
            code=269,
            message=paste0(
                "Error in io::fn_load_phenotype(...). ",
                "Phenotype file: ", fname_pheno, ", contains non-numeric data at column ", idx_col_y, ". ",
                "Are you certain that your file is separated by: '", sep, "'? ",
                "Are the phenotype data really at column '", idx_col_y, "'? ",
                "Are your missing data encoded as any of these: ", paste(paste0("'", na_strings, "'"), collapse=", "), "? ",
                "Do these look like numbers to you: ", paste(utils::head(y), collapse=", "), "? ",
                "These too: ", paste(utils::tail(y), collapse=", "), "?"))
        return(error)
    }
    ### Emit an error if there is no phenotypic variance
    if (stats::var(y, na.rm=TRUE) < .Machine$double.eps) {
        error = methods::new("gpError",
            code=270,
            message=paste0(
                "Error in io::fn_load_phenotype(...). ",
                "No variance in phenotype data. ",
                "We require variance because without it, we are lost in the dark without even a match to guide us out."))
        return(error)
    }
    names(y) = entry
    if (verbose) {
        print("Distribution of phenotype data:")
        txtplot::txtdensity(y[!is.na(y)])
    }
    if (header) {
        trait_name = colnames(df)[idx_col_y]
    } else {
        trait_name = paste0("trait_", idx_col_y)
    }
    return(list(y=y, pop=pop, trait_name=trait_name))
}

#' Filter phenotype data by converting outliers into missing data and optionally excluding missing data
#'
#' @param list_pheno list with 3 elements: 
#'  (1) $y: a named vector of numeric phenotype data, 
#'  (2) $pop: population or groupings corresponding to each element of y, and
#'  (3) $trait_name: name of the trait.
#' @param remove_outliers remove missing data in the y vector? 
#'  If true, this removes the corresponding element/s in the pop vector. (Default=TRUE)
#' @param remove_NA remove missing data in the y vector? 
#'  If true, this removes the corresponding element/s in the pop vector. (Default=FALSE)
#' @param verbose show phenotype filtering messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'      + $y: named vector of numeric phenotype data
#'      + $pop: population or groupings corresponding to each element of y
#'      + $trait_name: name of the trait
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' list_pheno$y[1] = Inf
#' list_pheno$y[2] = NA
#' list_pheno_filtered = fn_filter_phenotype(list_pheno, remove_NA=TRUE, verbose=TRUE)
#' @export
fn_filter_phenotype = function(list_pheno, remove_outliers=TRUE, remove_NA=FALSE, verbose=FALSE) {
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # verbose = TRUE
    ###################################################
    if (verbose) {
        print("#############################")
        print("### Filter phenotype data ###")
        print("#############################")
    }
    if (methods::is(list_pheno, "gpError")) {
        error = chain(list_pheno,
            methods::new("gpError", 
                code=271,
                message=paste0(
                    "Error in io::fn_filter_phenotype(...). ",
                    "The loaded phenotype data returned an error."
                ))
        )
        return(error)
    }
    ### Check lengths of the phenotype data and population/grouping vector
    n = length(list_pheno$y)
    if (n != length(list_pheno$pop)) {
        error = methods::new("gpError",
            code=272,
            message=paste0(
                "Error in io::fn_filter_phenotype(...). ",
                "The length of the phenotype vector (n=", n, ") is not equal to ",
                "the length of population/grouping vector (n=",length(list_pheno$pop) , "). ",
                "This is invalid. Please make sure the phenotype vector (y) correspond element-wise to the ",
                "population/grouping vector (pop)."))
        return(error)
    }
    ### Identify outliers with graphics::boxplot, i.e. values beyond -2.698 standard deviations (definition of R::graphics::boxplot whiskers)
    if (remove_outliers) {
        b = graphics::boxplot(list_pheno$y, plot=FALSE)
        vec_idx = which(list_pheno$y %in% b$out)
        if (length(vec_idx) == 0) {
            if (verbose) {"The phenotype data do not have any outliers."}
        } else {
            if (verbose) {
                print("Before removing outlier/s:")
                print(paste0("n=", n))
                txtplot::txtdensity(list_pheno$y[!is.na(list_pheno$y) & !is.infinite(list_pheno$y)])
            }
            list_pheno$y[vec_idx] = NA
            if (verbose) {
                print("After removing outlier/s:")
                print(paste0("n=", length(list_pheno$y)))
                txtplot::txtdensity(list_pheno$y[!is.na(list_pheno$y) & !is.infinite(list_pheno$y)])
            }
        }
    }
    ### Removing missing data
    if (remove_NA) {
        vec_idx = !is.na(list_pheno$y)
        if (length(vec_idx) == n) {
            if (verbose) {"The phenotype data do not have missing values."}
        } else {
            if (verbose) {
                print("Before removing missing values:")
                print(paste0("n=", n))
                txtplot::txtdensity(list_pheno$y[!is.na(list_pheno$y) & !is.infinite(list_pheno$y)])
            }
            list_pheno$y = list_pheno$y[vec_idx]
            list_pheno$pop = list_pheno$pop[vec_idx]
            if (verbose) {
                print("After removing missing values:")
                print(paste0("n=", length(list_pheno$y)))
                txtplot::txtdensity(list_pheno$y[!is.na(list_pheno$y) & !is.infinite(list_pheno$y)])
            }
        }
    }
    ### Emit an error if there is no phenotypic variance after filtering
    if (stats::var(list_pheno$y, na.rm=TRUE) < .Machine$double.eps) {
        error = methods::new("gpError",
            code=273,
            message=paste0(
                "Error in io::fn_filter_phenotype(...). ",
                "No variance in phenotype data after filtering. ",
                "Consider transforming your phenotype data."))
        return(error)
    }
    return(list_pheno)
}

#' Save phenotype data into a delimited file
#'  with 3 columns referring to the sample ID, population ID, and trait values
#'
#' @param list_pheno list with 3 elements: 
#'  (1) $y: a named vector of numeric phenotype data, 
#'  (2) $pop: population or groupings corresponding to each element of y, and
#'  (3) $trait_name: name of the trait.
#' @param fname file name of the output file
#' @param sep delimited of the output file (Default="\\t")
#' @param verbose show phenotype filtering messages? (Default=FALSE)
#' @returns
#'  - Ok: file name of the output file
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' fn_save_phenotype(list_pheno, fname=tempfile(fileext=".tsv"), sep="\t", verbose=TRUE)
#' fn_save_phenotype(list_pheno, fname=tempfile(fileext=".csv"), sep=",", verbose=TRUE)
#' @export
fn_save_phenotype = function(list_pheno, fname, sep="\t", verbose=FALSE) {
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # fname = "test.tsv"
    # sep = "\t"
    # verbose = TRUE
    ###################################################
    if (methods::is(list_pheno, "gpError")) {
        error = chain(list_pheno,
            methods::new("gpError", 
                code=274,
                message=paste0(
                    "Error in io::fn_save_phenotype(...). ",
                    "The loaded phenotype data returned an error."
                ))
        )
        return(error)
    }
    ### Check lengths of the phenotype data and population/grouping vector
    n = length(list_pheno$y)
    if (n != length(list_pheno$pop)) {
        error = methods::new("gpError",
            code=275,
            message=paste0(
                "Error in io::fn_save_phenotype(...). ",
                "The length of the phenotype vector (n=", n, ") is not equal to ",
                "the length of population/grouping vector (n=",length(list_pheno$pop) , "). ",
                "This is invalid. Please make sure the phenotype vector (y) correspond element-wise to the ",
                "population/grouping vector (pop)."))
        return(error)
    }
    ### Save
    df = data.frame(id=names(list_pheno$y), pop=list_pheno$pop, trait=list_pheno$y)
    utils::write.table(df, file=fname, sep=sep, row.names=FALSE, col.names=TRUE, quote=FALSE)
    if (verbose) {print(paste0("Phenotype data saved into: ", fname))}
    ### Return output file name
    return(fname)
}

#' Merge genotype and phenotype data by their entry names, i.e. rownames for G and names for y.
#' All samples with genotype data will be included and samples without phenotype data 
#' will be set to NA (all.x=TRUE). Samples with phenotype but without genotype data are omitted.
#'
#' @param G numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#' @param list_pheno list with 3 elements: 
#'  (1) $y: a named vector of numeric phenotype data, 
#'  (2) $pop: population or groupings corresponding to each element of y, and
#'  (3) $trait_name: name of the trait.
#' @param COVAR numeric n samples x k covariates matrix with non-null row and column names.
#'  The row names correspond to the sample or entry or pool names,
#'  which fully or partially match with the row names in G and the names in list_pheno$y. (Default=NULL)
#' @param verbose show genotype, phenotype, and covariate merging messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'      + $G: numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'      Row names can be any string of characters which identify the sample or entry or pool names.
#'      Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'      the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'      subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'      + $list_pheno:
#'      $y: named vector of numeric phenotype data
#'      $pop: population or groupings corresponding to each element of y
#'      $trait_name: name of the trait
#'      + $COVAR: numeric n samples x k covariates matrix with non-null row and column names.
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' rownames(G)[1] = "entry_exclude_me"
#' rownames(G)[2] = "entry_exclude_me_too"
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' COVAR = matrix(stats::rnorm(n=(10*nrow(G))), nrow=nrow(G))
#' rownames(COVAR) = rownames(G); colnames(COVAR) = paste0("covariate_", 1:ncol(COVAR))
#' list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
#' @export
fn_merge_genotype_and_phenotype = function(G, list_pheno, COVAR=NULL, verbose=FALSE) {
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    # rownames(G)[1] = "entry_exclude_me"
    # rownames(G)[2] = "entry_exclude_me_too"
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # COVAR = matrix(stats::rnorm(n=(10*nrow(G))), nrow=nrow(G))
    # rownames(COVAR) = rownames(G); colnames(COVAR) = paste0("covariate_", 1:ncol(COVAR))
    # verbose = TRUE
    ###################################################
    if (verbose) {
        print("###########################################")
        print("### Merging genotype and phenotype data ###")
        print("###########################################")
    }
    ### Check list_pheno
    if (sum(is.na(list_pheno$y))==length(list_pheno$y)) {
        error = methods::new("gpError",
            code=284,
            message=paste0(
                "Error in io::fn_merge_genotype_and_phenotype(...). ",
                "All phenotype data are missing."))
        return(error)
    }
    if (length(names(list_pheno$y)) < length(list_pheno$y)) {
        error = methods::new("gpError",
            code=285,
            message=paste0(
                "Error in io::fn_merge_genotype_and_phenotype(...). ",
                "Phenotype data are missing names."))
        return(error)
    }
    ### All samples with genotype data will be included and samples without phenotype data will be set to NA (all.x=TRUE)
    ### Samples with phenotype but without genotype data are omitted.
    M = merge(
        data.frame(id=rownames(G), G, check.names=FALSE), 
        data.frame(id=names(list_pheno$y), pop=list_pheno$pop, y=list_pheno$y), 
        by="id", all.x=TRUE)
    if (sum(!is.na(M$y)) == 0) {
        error = methods::new("gpError",
            code=276,
            message=paste0(
                "Error in io::fn_merge_genotype_and_phenotype(...). ",
                "None of the samples/entries/pools with genotype data have phenotype data. ",
                "IDs in the genotype data: ", paste(c(utils::head(rownames(G)), "... "), collapse=", "),
                "IDs in the phenotype data: ", paste(c(utils::head(names(list_pheno$y)), "... "), collapse=", ")))
        return(error)
    }
    if (!is.null(COVAR)) {
        vec_rownames = rownames(COVAR)
        vec_colnames = colnames(COVAR)
        if (is.null(vec_rownames) | is.null(vec_colnames)) {
            error = methods::new("gpError",
                code=277,
                message=paste0(
                    "Error in io::fn_merge_genotype_and_phenotype(...). ",
                    "The covariance matrix (COVAR) need to have both row and column names."))
            return(error)
        }
        ### Append 'covariates_' to the column names of COVAR if some names intersect with the column names of M
        if (sum(colnames(M) %in% vec_colnames) > 0) {
            colnames(COVAR) = paste0("covariates_", vec_colnames)
        }
        M = merge(M, data.frame(id=vec_rownames, COVAR), by="id", all.x=TRUE)
        vec_idx_G = which(!(colnames(M) %in% c("id", "pop", "y", gsub("-", ".", colnames(COVAR))))) ### merging converts all the dashes into dots in the column names
        vec_idx_C = which(colnames(M) %in% gsub("-", ".", colnames(COVAR))) ### merging converts all the dashes into dots in the column names
        COVAR = as.matrix(M[, vec_idx_C])
        rownames(COVAR) = M$id
        colnames(COVAR) = vec_colnames ### Revert to original column names
    } else {
        vec_idx_G = which(!(colnames(M) %in% c("id", "pop", "y")))
        COVAR = NULL
    }
    G = as.matrix(M[, vec_idx_G])
    rownames(G) = M$id
    y = M$y
    names(y) = M$id
    pop = M$pop
    pop[is.na(pop)] = "Unknown"
    trait_name = list_pheno$trait_name
    ### Clean-up
    rm("M")
    gc()
    if (verbose) {
        print(paste0("n=", nrow(G), "; p=", ncol(G)))
        print("Allele frequency distribution: ")
        vec_freqs_sample = sample(G, size=min(c(prod(dim(G)), 1e4)))
        vec_freqs_sample = vec_freqs_sample[!is.na(vec_freqs_sample)]
        txtplot::txtdensity(c(vec_freqs_sample, 1-vec_freqs_sample))
        print("Phenotype distribution:")
        txtplot::txtdensity(y[!is.na(y)])
        if (!is.null(COVAR)) {
            print("Covariate distribution:")
            print(paste0("m=", ncol(COVAR)))
            txtplot::txtdensity(COVAR[!is.na(COVAR)])
            print("Clustering of the samples/entries/pools based on the covariate matrix (COVAR):")
            if (ncol(COVAR) == 2) {
                txtplot::txtplot(COVAR[,1], COVAR[,2])
            } else {
                print("Using the first 2 PCs of COVAR:")
                list_PCs = stats::prcomp(COVAR)
                print(summary(list_PCs)$importance[, 1:2])
                txtplot::txtplot(list_PCs$x[,1], list_PCs$x[,2])
                list_PCs = NULL
            }
        } else {
            print("Covariate is null.")
        }
    }
    ### Output
    return(list(G=G, list_pheno=list(y=y, pop=pop, trait_name=trait_name), COVAR=COVAR))
}

#' Subset the list of merged genotype and phenotype data using a vector of sample/entry/pool indexes.
#'
#' @param list_merged list of merged genotype matrix, and phenotype vector, as well as an optional covariate matrix
#'      + $G: numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'      Row names can be any string of characters which identify the sample or entry or pool names.
#'      Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'      the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'      subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'      + $list_pheno:
#'      $y: named vector of numeric phenotype data
#'      $pop: population or groupings corresponding to each element of y
#'      $trait_name: name of the trait
#'      + $COVAR: numeric n samples x k covariates matrix with non-null row and column names.
#' @param vec_idx numeric vector of sample/entry/pool indexes to extract from list_merged
#' @param verbose show genotype, phenotype, and covariate subsetting messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'      + $G: numeric n samples x p loci-alleles matrix of allele frequencies with non-null row and column names.
#'      Row names can be any string of characters which identify the sample or entry or pool names.
#'      Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'      the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'      subsequent elements are optional which may refer to the allele identifier and other identifiers.
#'      + $list_pheno:
#'      $y: named vector of numeric phenotype data
#'      $pop: population or groupings corresponding to each element of y
#'      $trait_name: name of the trait
#'      + $COVAR: numeric n samples x k covariates matrix with non-null row and column names.
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' rownames(G)[1] = "entry_exclude_me"
#' rownames(G)[2] = "entry_exclude_me_too"
#' list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
#' COVAR = matrix(stats::rnorm(n=(10*nrow(G))), nrow=nrow(G))
#' rownames(COVAR) = rownames(G); colnames(COVAR) = paste0("covariate_", 1:ncol(COVAR))
#' list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
#' vec_idx = which(list_merged$list_pheno$pop == list_merged$list_pheno$pop[1])
#' list_merged_subset = fn_subset_merged_genotype_and_phenotype(list_merged=list_merged,
#'  vec_idx=vec_idx, verbose=TRUE)
#' @export
fn_subset_merged_genotype_and_phenotype = function(list_merged, vec_idx, verbose=FALSE) {
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(n_pop=3, verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # COVAR = G %*% t(G)
    # list_merged = fn_merge_genotype_and_phenotype(G=G, list_pheno=list_pheno, COVAR=COVAR, verbose=TRUE)
    # vec_idx = which(list_merged$list_pheno$pop == list_merged$list_pheno$pop[1])
    # verbose = TRUE
    ###################################################
    ### Input sanity check
    if (methods::is(list_merged, "gpError")) {
        error = chain(list_merged, 
            methods::new("gpError",
                code=278,
                message=paste0(
                    "Error in cross_validation::fn_subset_merged_genotype_and_phenotype(...). ",
                    "Input data (list_merged) is an error type."
                )))
        return(error)
    }
    if (sum(c(1:nrow(list_merged$G)) %in% vec_idx) != length(vec_idx)) {
        error = methods::new("gpError",
            code=279,
            message=paste0(
                "Error in cross_validation::fn_subset_merged_genotype_and_phenotype(...). ",
                "The indexes of samples/entries/pools do not match the indexes in the data set. ",
                "The indexes asked for ranges from ", min(vec_idx), " to ", max(vec_idx), " while ",
                "the indexes in the data set ranges from 1 to ", nrow(list_merged$G), "."
            ))
        return(error)
    }
    ## Subset
    G = list_merged$G[vec_idx, , drop=FALSE]
    list_pheno = list(
        y=list_merged$list_pheno$y[vec_idx],
        pop=list_merged$list_pheno$pop[vec_idx],
        trait_name=list_merged$list_pheno$trait_name
    )
    if (!is.null(list_merged$COVAR)) {
        COVAR = list_merged$COVAR[vec_idx, , drop=FALSE]
    } else {
        COVAR = NULL
    }
    ### Clean-up
    rm("list_merged")
    gc()
    ### Output
    return(list(
        G=G,
        list_pheno=list_pheno,
        COVAR=COVAR
    ))
}

#' Estimate memory usage for parallel replicated k-fold cross validation of multiple genomic prediction models
#'
#' @param X any R object but the main intended object is a list containing the 
#'  genotype matrix, phenotype list and covariate matrix for genomic prediction
#' @param n_models number of genomic prediction models to fit. Note that Bayesian and 
#'  gBLUP models are more memory-intensive than penalised regression ones. (Default=7)
#' @param n_folds number of cross-validation folds (Default=10)
#' @param n_reps number of cross-validation replication (Default=10)
#' @param memory_requested_Gb memory requested or available for use (Default=400)
#' @param memory_multiplier estimated memory usage multiplier (Default=40)
#' @param verbose show memory usage estimation messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'      + $size_X: memory used for a single genomic prediction instance
#'      + $size_total: total memory required for parallel computations across models, folds and replications
#'      + $n_threads: recommended and estimated maximum number of threads to use to prevent out-of-memory (OOM) error
#'  - Err: gpError
#' @examples
#' list_mem = fn_estimate_memory_footprint(X=rnorm(10000), verbose=TRUE)
#' @export
fn_estimate_memory_footprint = function(X, n_models=7, n_folds=10, n_reps=10, 
    memory_requested_Gb=400, memory_multiplier=40, verbose=FALSE)
{
    ###################################################
    ### TEST
    # X = matrix(0.0, nrow=492, ncol=455255)
    # n_models = 7
    # n_reps = 10
    # n_folds = 10
    # memory_requested_Gb = 400
    # memory_multiplier = 50
    # verbose = TRUE
    ###################################################
    if ((as.numeric(utils::object.size(X)) == 0) | (n_models <= 0) | (n_folds <= 0) | (n_reps <= 0) | (memory_requested_Gb <= 0) | (memory_multiplier <= 0)) {
        error = methods::new("gpError",
            code=280,
            message=paste0(
                "Error in io::fn_estimate_memory_footprint(...). ",
                "The size of the input data, number of models, folds and replications, ",
                "as well as the memory requested or available and memory usage multiplier cannot be zero. ",
                "X (size=", utils::object.size(X), " bytes); ", 
                "n_models=", n_models, "; n_folds=", n_folds, "; n_reps=", n_reps, 
                "; memory_requested_Gb=", memory_requested_Gb, " Gb; ",
                "memory_multiplier=", memory_multiplier))
        return(error)
    }
    size_RAM = memory_requested_Gb * 2^30 * (utils::object.size(1) / utils::object.size(1))
    size_X = utils::object.size(X)
    size_total = size_X * n_models * n_folds * n_reps * memory_multiplier
    n_threads = floor(as.numeric(gsub(" bytes", "",  size_RAM / (size_X * memory_multiplier))))
    if (verbose) {
        print(paste0("Size of X: = ", format(size_X, units="b")))
        print(paste0("n_models = ", n_models, "; n_folds = ", n_folds, "; n_reps = ", n_reps))
        print(paste0("Total memory requested = ", 
            format(size_RAM, units="Gb"), " (", 
            format(size_RAM, units="Mb"), "; ", 
            format(size_RAM, units="Kb"), "; ", 
            format(size_RAM, units="b"), ")"))
        print(paste0("Size of each dataset (merged genotype, phenotype and covariate data) = ", 
            format(size_X, units="Gb"), " (", 
            format(size_X, units="Mb"), "; ", 
            format(size_X, units="Kb"), "; ", 
            format(size_X, units="b"), ")"))
        print(paste0("Estimated total memory required for simultaneous parallel computations = ", 
            format(size_total, units="Gb"), " (", 
            format(size_RAM, units="Mb"), "; ", 
            format(size_RAM, units="Kb"), "; ", 
            format(size_RAM, units="b"), ")"))
        print(paste0("To prevent out-of-memory (OOM) errors, we recommend a maximum number of threads = ", n_threads))
    }
    return(list(
        size_X=size_X,
        size_total=size_total,
        n_threads=n_threads))
}

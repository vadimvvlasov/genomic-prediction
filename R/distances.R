# source("R/io.R")

#' Genetic relationship matrices and their inverses
#'
#' @param G n samples x p loci matrix of allele frequencies (numeric ranging from 0 to 1) with non-null row and column names.
#'  Row names can be any string of characters which identify the sample or entry or pool names.
#'  Column names need to be tab-delimited, where first element refers to the chromosome or scaffold name, 
#'  the second should be numeric which refers to the position in the chromosome/scaffold, and 
#'  subsequent elements are optional which may refer to the allele identifier and other identifiers.
#' @param maf minimum allele frequency (Default=0.01)
#' @param ploidy (for the GRM of Bell et al, 2017 & VanRaden et al, 2008) ploidy level which can refer to the number of haploid genomes to simulate pools (Default=2)
#' @param diagonal_load a small positive value to add to the diagonals of the genetic relationship matrix to ensure invertibility (Default=0.001)
#' @param verbose show messages? (Default=FALSE)
#' @returns
#'  - Ok (returns several genetic relationship and distances matrices):
#'      + $grm: simple X%*%t(X)/ncol(X)
#'      + $grm_VanRaden: ploidy-awaren (see VanRaden et al, 2008 and Bell et al, 2017)
#'      + $dist_binary: Jaccard's distance, where non-zero values are converted to ones.
#'      + $dist_euclidean: Euclidean's distance (L2 norm)
#'      + $dist_taxicab: Taxicab or Manhattan distance (L1 norm)
#'      + $inverse_grm: inverse of $grm
#'      + $inverse_grm_VanRaden: inverse of $grm_VanRaden
#'      + $inverse_one_minus_dist_binary: inverse of 1.00 - $dist_binary
#'      + $inverse_one_minus_dist_euclidean: inverse of 1.00 - $dist_euclidean
#'      + $inverse_one_minus_dist_taxicab: inverse of 1.00 - $dist_taxicab
#'  - Err: gpError
#' @examples
#' list_sim = fn_simulate_data(verbose=TRUE)
#' G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
#' list_grm_and_dist = fn_grm_and_dist(G=G, verbose=TRUE)
#' @export
fn_grm_and_dist = function(G, maf=0.01, ploidy=2, diagonal_load=0.001, verbose=FALSE) {
    ###################################################
    ### TEST
    # list_sim = fn_simulate_data(verbose=TRUE)
    # G = fn_load_genotype(fname_geno=list_sim$fname_geno_vcf)
    # # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # maf = 0.01
    # ploidy = 2
    # diagonal_load = 0.001
    # verbose = TRUE
    ###################################################
    if (methods::is(G, "gpError")) {
        error = chain(G, 
            methods::new("gpError",
                code=600,
                message=paste0(
                    "Error in distances::fn_grm_and_dist(...). ",
                    "Input data (G) is an error type."
                )))
        return(error)
    }
    ### Retain only the loci (columns) which passed the minimum allele frequency (maf)
    n = nrow(G) ### number of samples
    p = ncol(G) ### number of markers
    q = colMeans(G)
    if (sum(is.na(q)) > 0) {
        vec_loci_names_with_missing_data = colnames(G)[is.na(q)]
        error = methods::new("gpError",
                code=601,
                message=paste0(
                    "We expect no missing data in the genotype matrix. ",
                    "See the following loci: ", 
                    paste(vec_loci_names_with_missing_data, collapse=", ")
                ))
        return(error)
    }
    vec_idx = which((q >= maf) & (q <= (1.00-maf)))
    if (length(vec_idx) < p) {
        G = G[, vec_idx, drop=FALSE]
        q = q[vec_idx]
    }
    ### (1) Simple GRM
    if (verbose) {print("Calculating GRM via G%*%t(G) / p:")}
    grm = G%*%t(G) / p
    inverse_grm = tryCatch(solve(grm), error=function(e) {NA})
    if (is.na(inverse_grm[1])) {
        diag(grm) = diag(grm) + diagonal_load
        inverse_grm = solve(grm)
    }
    if (verbose) {
        if (abs(sum(grm %*% inverse_grm) - n) > 1e-5) {
            print("Warning: the inverse of grm may be unreliable due to the addition of the diagonal_load to enable inversion.")
        }
    }
    ### (2) Ploidy-aware GRM of Bell et al (2017) and VanRaden et al (2008)
    if (verbose) {print("Calculating GRM via Bell et al (2017) and VanRaden et al (2008):")}
    G_star = ploidy*(G-0.5)
    q_star = ploidy*(q-0.5)
    Z = G_star - matrix(rep(q_star, each=n), nrow=n, byrow=FALSE)
    grm_VanRaden = Z%*%t(Z) / (ploidy*sum(q*(1-q)))
    inverse_grm_VanRaden = tryCatch(solve(grm_VanRaden), error=function(e) {NA})
    if (is.na(inverse_grm_VanRaden[1])) {
        diag(grm_VanRaden) = diag(grm_VanRaden) + diagonal_load
        inverse_grm_VanRaden = solve(grm_VanRaden)
    }
    if (verbose) {
        if (abs(sum(grm_VanRaden %*% inverse_grm_VanRaden) - n) > 1e-5) {
            print("Warning: the inverse of grm_VanRaden may be unreliable due to the addition of the diagonal_load to enable inversion.")
        }
    }
    ### (3) Jaccard's distance
    if (verbose) {print("Calculating Jaccard's distances:")}
    dist_binary = as.matrix(stats::dist(G, method="binary", diag=TRUE, upper=TRUE))
    one_minus_dist_binary =  1.00 - dist_binary
    inverse_one_minus_dist_binary = tryCatch(solve(one_minus_dist_binary), error=function(e) {NA})
    if (is.na(inverse_one_minus_dist_binary[1])) {
        diag(one_minus_dist_binary) = diag(one_minus_dist_binary) + diagonal_load
        inverse_one_minus_dist_binary = solve(one_minus_dist_binary)
    }
    if (verbose) {
        if (abs(sum(one_minus_dist_binary %*% inverse_one_minus_dist_binary) - n) > 1e-5) {
            print("Warning: the inverse of grm may be unreliable due to the addition of the diagonal_load to enable inversion.")
        }
    }
    ### (4) Euclidean's distance
    if (verbose) {print("Calculating Euclidean distances:")}
    dist_euclidean = as.matrix(stats::dist(G, method="euclidean", diag=TRUE, upper=TRUE))
    one_minus_dist_euclidean =  1.00 - dist_euclidean
    inverse_one_minus_dist_euclidean = tryCatch(solve(one_minus_dist_euclidean), error=function(e) {NA})
    if (is.na(inverse_one_minus_dist_euclidean[1])) {
        diag(one_minus_dist_euclidean) = diag(one_minus_dist_euclidean) + diagonal_load
        inverse_one_minus_dist_euclidean = solve(one_minus_dist_euclidean)
    }
    if (verbose) {
        if (abs(sum(one_minus_dist_euclidean %*% inverse_one_minus_dist_euclidean) - n) > 1e-5) {
            print("Warning: the inverse of grm may be unreliable due to the addition of the diagonal_load to enable inversion.")
        }
    }
    ### (5) Taxicab distance
    if (verbose) {print("Calculating Taxicab distances:")}
    dist_taxicab = as.matrix(stats::dist(G, method="manhattan", diag=TRUE, upper=TRUE))
    one_minus_dist_taxicab =  1.00 - dist_taxicab
    inverse_one_minus_dist_taxicab = tryCatch(solve(one_minus_dist_taxicab), error=function(e) {NA})
    if (is.na(inverse_one_minus_dist_taxicab[1])) {
        diag(one_minus_dist_taxicab) = diag(one_minus_dist_taxicab) + diagonal_load
        inverse_one_minus_dist_taxicab = solve(one_minus_dist_taxicab)
    }
    if (verbose) {
        if (abs(sum(one_minus_dist_taxicab %*% inverse_one_minus_dist_taxicab) - n) > 1e-5) {
            print("Warning: the inverse of grm may be unreliable due to the addition of the diagonal_load to enable inversion.")
        }
    }
    # ### TESTS WITH ASREML
    # system(command="module load ASReml-R", intern=TRUE)
    # library(asreml)
    # GRM = inverse_one_minus_dist_taxicab
    # list_pheno = fn_load_phenotype(fname_pheno=list_sim$fname_pheno_tsv)
    # df = data.frame(y=list_pheno$y, gen=names(list_pheno$y))
    # rownames(df) = NULL
    # df$gen = as.factor(df$gen)
    # str(df)
    # ### With spherical error (e~N(0,se))
    # mod_6 = asreml::asreml(
    #     fixed = y ~ 1,
    #     random = ~giv(gen),
    #     ginverse = list(gen=GRM),
    #     rcov = ~idv(units),
    #     data = df,
    #     workspace = 3e9
    # )
    # coef_6 = coef(mod_6)$random
    # pred_6 = predict(mod_6, classify="gen")$pred$pvals$predicted.value
    # stats::cor(coef_6, pred_6)
    # summary(mod_6)
    # ### With standard normal error (e~N(0,1))
    # mod_7 = asreml::asreml(
    #     fixed = y ~ 1,
    #     random = ~giv(gen),
    #     ginverse = list(gen=GRM),
    #     rcov = ~id(units),
    #     data = df,
    #     workspace = 3e9
    # )
    # coef_7 = coef(mod_7)$random
    # pred_7 = predict(mod_7, classify="gen")$pred$pvals$predicted.value
    # stats::cor(coef_7, pred_7)
    # summary(mod_7)
    # ### Correlate the BLUPs of the 2 models
    # txtplot::txtdensity(coef_6)
    # txtplot::txtdensity(coef_7)
    # txtplot::txtplot(coef_6, coef_7)
    # stats::cor(coef_6, coef_7)
    # mean(abs(coef_6 - coef_7))
    # sqrt(mean((coef_6 - coef_7)^2))
    if (verbose) {
        vec_grm = as.vector(grm)
        vec_grm_VanRaden = as.vector(grm_VanRaden)
        vec_dist_binary = as.vector(dist_binary)
        vec_dist_euclidean = as.vector(dist_euclidean)
        vec_dist_taxicab = as.vector(dist_taxicab)
        print("grm vs grm_VanRaden:")
        txtplot::txtplot(vec_grm, vec_grm_VanRaden)
        print("grm vs (binary, euclidean, taxicab):")
        txtplot::txtplot(vec_grm, vec_dist_binary)
        txtplot::txtplot(vec_grm, vec_dist_euclidean)
        txtplot::txtplot(vec_grm, vec_dist_taxicab)
        print("grm_VanRaden vs (binary, euclidean, taxicab):")
        txtplot::txtplot(vec_grm_VanRaden, vec_dist_binary)
        txtplot::txtplot(vec_grm_VanRaden, vec_dist_euclidean)
        txtplot::txtplot(vec_grm_VanRaden, vec_dist_taxicab)
        print("Correlations:")
        mat_grms_and_dists = cbind(vec_grm, vec_grm_VanRaden, vec_dist_binary, vec_dist_euclidean, vec_dist_taxicab)
        colnames(mat_grms_and_dists) = c("grm", "grm_VanRaden", "dist_binary", "dist_euclidean", "dist_taxicab")
        print(stats::cor(mat_grms_and_dists))
    }
    return(list(
        grm=grm,
        grm_VanRaden=grm_VanRaden,
        dist_binary=dist_binary,
        dist_euclidean=dist_euclidean,
        dist_taxicab=dist_taxicab,
        inverse_grm=inverse_grm,
        inverse_grm_VanRaden=inverse_grm_VanRaden,
        inverse_one_minus_dist_binary=inverse_one_minus_dist_binary,
        inverse_one_minus_dist_euclidean=inverse_one_minus_dist_euclidean,
        inverse_one_minus_dist_taxicab=inverse_one_minus_dist_taxicab
    ))
}

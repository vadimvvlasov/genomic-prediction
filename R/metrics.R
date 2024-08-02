# source("R/io.R")

#' Genomic prediction metrics
#' NOTES:
#'  - Add more metrics which we may be interested in. while keeping the function signatures, i.e. input and output consistent
#'  - The number of non-zero estimated marker effects are extracted within each model, i.e. in `models.R`
#'  - Run time of each genomic prediction and validation on a single set is measured during each cross-validation replicate, i.e in `cross-validation.R`
#'
#' @param y_true numeric vector observed phenotype values
#' @param y_pred numeric vector predicted phenotype values
#' @param verbose show genomic prediction performance metric calculation messages? (Default=FALSE)
#' @returns
#'  - Ok:
#'      + $mbe: mean bias error
#'      + $mae: mean absolute error
#'      + $rmse: root mean squared error
#'      + $r2: coefficient of determination
#'      + $corr: Pearson's product moment correlation
#'      + $power_t10: fraction of observed top 10 phenotype values correctly predicted
#'      + $power_b10: fraction of observed bottom 10 phenotype values correctly predicted
#'      + $var_additive: variance of predicted phenotype values (estimator of additive genetic variance)
#'      + $var_residual: variance of differnce between observed and predicted phenotype values (estimator of residual variance)
#'      + $h2: estimated narrow-sense heritability weighted by Pearson's correlation between observed 
#'          and predicted phenotype values.
#'  - Err: gpError
#' @examples
#' y_pred = stats::rnorm(100)
#' y_true = y_pred + stats::rnorm(100)
#' list_metrics = fn_prediction_performance_metrics(y_true=y_true, y_pred=y_pred, verbose=TRUE)
#' @export
fn_prediction_performance_metrics = function(y_true, y_pred, verbose=FALSE) {
    ###################################################
    ### TEST
    # y_true = stats::rnorm(100)
    # y_pred = stats::rnorm(100); # y_pred = y_true + stats::rnorm(n=100, mean=0, stats::sd=0.1)
    # # y_pred = rep(NA, length(y_true))
    # verbose = TRUE
    ###################################################
    y_true = as.vector(y_true)
    y_pred = as.vector(y_pred)
    n = length(y_true)
    if (n != length(y_pred)) {
        error = methods::new("gpError",
            code=500,
            message=paste0(
                "Error in metrics::fn_prediction_performance_metrics(...). ",
                "The vector of observed phenotypes has a length of ", n,
                " which does not match the length of the vector of predicted phenotypes with ",
                length(y_pred), " elements."
            ))
        return(error)
    }
    error = y_true - y_pred
    mbe = mean(error, na.rm=TRUE)
    mae = mean(abs(error), na.rm=TRUE)
    rmse = sqrt(mean(error^2, na.rm=TRUE))
    if (sum(is.na(error)) == n) {
        r2 = NA
    } else {
        r2 = 1 - (sum(error^2, na.rm=TRUE) / sum((y_true-mean(y_true, na.rm=TRUE))^2, na.rm=TRUE))
    }
    corr = suppressWarnings(stats::cor(y_true, y_pred, method="pearson", use="na.or.complete"))
    ### Power to select true top and bottom 10%
    n_top_or_bottom_10 = max(c(1, round(0.1*n)))
    top10_dec_true = order(y_true, decreasing=TRUE)[1:n_top_or_bottom_10]
    bottom10_dec_true = order(y_true, decreasing=FALSE)[1:n_top_or_bottom_10]
    if (sum(is.na(error)) == n) {
        power_t10 = NA
        power_b10 = NA
    } else {
        top10_dec_pred = order(y_pred, decreasing=TRUE)[1:n_top_or_bottom_10]
        bottom10_dec_pred = order(y_pred, decreasing=FALSE)[1:n_top_or_bottom_10]
        power_t10 = sum((top10_dec_true %in% top10_dec_pred)) / n_top_or_bottom_10
        power_b10 = sum((bottom10_dec_true %in% bottom10_dec_pred)) / n_top_or_bottom_10
    }
    ### Narrow-sense heritability scaled by prediction accuracy in terms of Pearson's correlation
    var_additive = stats::var(y_pred, na.rm=TRUE)
    var_residual = stats::var(y_true-y_pred, na.rm=TRUE)
    h2 = round(var_additive/(var_additive+var_residual), 10)
    if (!is.na(h2)) {
        if ((h2 < 0.0) | (h2 > 1.0)) {
            h2 = NA
        }
    }
    if (verbose) {
        vec_idx = which(!is.na(y_true) & !is.na(y_pred))
        if (length(vec_idx) > 0) {
            print("Scatter plot of the observed and predicted phenotypes")
            txtplot::txtplot(x=y_true, y=y_pred)
        } else {
            print("All pairs of phenotype values have missing data.")
        }
        print(paste0("Mean bias error (mbe) = ", mbe))
        print(paste0("Mean absolute error (mae) = ", mae))
        print(paste0("Root mean square error (rmse) = ", rmse))
        print(paste0("Coefficient of determination (r2) = ", r2))
        print(paste0("Power to identify top 10 (power_t10) = ", power_t10))
        print(paste0("Power to identify bottom 10 (power_b10) = ", power_b10))
        print(paste0("Pearson's correlation (corr) = ", corr))
        print(paste0("Variance of predicted phenotypes (var_additive) = ", var_additive, " (estimator of additive genetic variance)"))
        print(paste0("Variance of the difference between observed and predicted phenotypes (var_residual) = ", var_residual, " (estimator of residual variance)"))
        print(paste0("Narrow-sense heritability estimate (h2) = ", h2))
    }
    ### Output
    return(list(
        mbe=mbe,
        mae=mae,
        rmse=rmse,
        r2=r2,
        corr=corr,
        power_t10=power_t10,
        power_b10=power_b10,
        var_additive=var_additive,
        var_residual=var_residual,
        h2=h2))
}

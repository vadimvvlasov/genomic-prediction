### Genomic prediction metrics
### Add more metrics which we may be interested in.
### NOTES:
###     - the number of non-zero estimated marker effects are extracted within each model, i.e. in `models.R`
###     - the run-time of each genomic predicition and validation on a single set is measured during each cross-validation replicate, i.e in `cross-validation.R`
library(testthat)

fn_prediction_performance_metrics = function(y_true, y_pred) {
    # y_true = rnorm(100)
    # y_pred = rnorm(100); # y_pred = y_true + rnorm(n=100, mean=0, sd=0.1)
    y_true = as.vector(y_true)
    y_pred = as.vector(y_pred)
    if (length(y_true) != length(y_pred)) {
        print("Please make sure the observed and predicted phenotype vectors are the same length.")
        return(list(mbe=NA, mae=NA, rmse=NA, r2=NA, corr=NA, , power_t10=NA, power_b10=NA))
    }
    error = y_true - y_pred
    mbe = mean(error)
    mae = mean(abs(error))
    rmse = sqrt(mean(error^2))
    r2 = 1 - (sum(error^2) / sum((y_true-mean(y_true))^2))
    corr = suppressWarnings(cor(y_true, y_pred, method="pearson"))
    if (is.na(corr)) {
        corr = 0.0
    }
    ### Power to select true top and bottom 10%
    n = length(y_true)
    n_top_or_bottom_10 = max(c(1, round(0.1*n)))
    top10_dec_true = order(y_true, decreasing=TRUE)[1:n_top_or_bottom_10]
    top10_dec_pred = order(y_pred, decreasing=TRUE)[1:n_top_or_bottom_10]
    power_t10 = sum((top10_dec_true %in% top10_dec_pred)) / n_top_or_bottom_10
    bottom10_dec_true = order(y_true, decreasing=FALSE)[1:n_top_or_bottom_10]
    bottom10_dec_pred = order(y_pred, decreasing=FALSE)[1:n_top_or_bottom_10]
    power_b10 = sum((bottom10_dec_true %in% bottom10_dec_pred)) / n_top_or_bottom_10
    ### Narrow-sense heritability scaled by prediction accuracy in terms of Pearson's correlation
    ### Assumptions:
    ###     - the predicted trait values are due to purely additive effects,
    ###     - the variance in the predicted trait values is proportional to the additive genetic variance (~Va)
    ###     - the variance in the true trait values represent the total phenotypic variance (Vp)
    ###     - the correlation between predicted and true trait values is the factor that relates the variance of the predicted trait values with the true heritability, i.e. h2_predicted = corr(true,pred) * h2_true
    h2 = round((var(y_pred) / var(y_true)) / corr, 10)
    if (is.na(h2) == FALSE) {
        if ((h2 < 0.0) | (h2 > 1.0)) {
            h2 = NA
        }
    } else {
        h2 = 0.0
    }
    ### Output
    out = list(mbe=mbe, mae=mae, rmse=rmse, r2=r2, corr=corr, power_t10=power_t10, power_b10=power_b10, h2=h2)
    return(out)
}

#####################################################################
############################### Tests ###############################
#####################################################################
tests_metrics = function() {
    test_that(
        "fn_prediction_performance_metrics", {
            print("fn_prediction_performance_metrics:")
            n = 100
            y1 = 1:n
            y2 = y1 + pi
            m0 = fn_prediction_performance_metrics(y1, y1)
            m1 = fn_prediction_performance_metrics(y1, y2)
            expect_equal(m0, list(mbe=0, mae=0, rmse=0, r2=1, corr=1, power_t10=1, power_b10=1, h2=1))
            expect_equal(m1, list(mbe=-pi, mae=pi, rmse=pi, r2=(1-(n*(pi^2)/sum((y1-mean(y1))^2))), corr=1, power_t10=1, power_b10=1, h2=1))
        }
    )
}

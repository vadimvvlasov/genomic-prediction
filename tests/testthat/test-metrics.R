# library(testthat)
# source("R/io.R")
# source("R/metrics.R")

test_that("fn_prediction_performance_metrics", {
    n = 100
    y_true = 1:n
    y_pred = y_true + pi
    list_metrics_1 = fn_prediction_performance_metrics(y_true, y_true)
    list_metrics_2 = fn_prediction_performance_metrics(y_true, y_pred)
    list_metrics_3 = fn_prediction_performance_metrics(y_true, NA)
    list_metrics_4 = fn_prediction_performance_metrics(y_true, rep(NA, length(y_true)))
    expect_equal(list_metrics_1, list(mbe=0, mae=0, rmse=0, r2=1, corr=1, power_t10=1, power_b10=1, var_additive=841.666666666666, var_residual=0.0, h2=1))
    expect_equal(list_metrics_2, list(mbe=-pi, mae=pi, rmse=pi, r2=(1-(n*(pi^2)/sum((y_true-mean(y_true))^2))), corr=1, power_t10=1, power_b10=1, var_additive=841.666666666666, var_residual=0.0, h2=1))
    expect_equal(methods::is(list_metrics_3, "gpError"), TRUE)
    expect_equal(is.na(list_metrics_4$mbe), TRUE)
    expect_equal(is.na(list_metrics_4$mae), TRUE)
    expect_equal(is.na(list_metrics_4$rmse), TRUE)
    expect_equal(is.na(list_metrics_4$r2), TRUE)
    expect_equal(is.na(list_metrics_4$corr), TRUE)
    expect_equal(is.na(list_metrics_4$power_t10), TRUE)
    expect_equal(is.na(list_metrics_4$power_b10), TRUE)
    expect_equal(is.na(list_metrics_4$var_additive), TRUE)
    expect_equal(is.na(list_metrics_4$var_residual), TRUE)
    expect_equal(is.na(list_metrics_4$h2), TRUE)
})

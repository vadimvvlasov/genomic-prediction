from __future__ import annotations

import numpy as np
from scipy import stats as sp_stats


def fn_prediction_performance_metrics(
    y_true, y_pred, *, verbose: bool = False
) -> dict:
    """Genomic prediction performance metrics.

    Mirrors the R implementation (R/metrics.R:fn_prediction_performance_metrics).

    Returns
    -------
    dict with keys:
        mbe, mae, rmse, r2, corr, corr_rank,
        power_t10, power_b10,
        var_additive, var_residual, h2
    """
    yt = np.asarray(y_true, dtype=float)
    yp = np.asarray(y_pred, dtype=float)

    mask = np.isfinite(yt) & np.isfinite(yp)
    yt = yt[mask]
    yp = yp[mask]
    n = yt.size

    if n == 0:
        raise ValueError("No valid y_true/y_pred pairs to score")

    diff = yt - yp

    # ── Basic error metrics ──────────────────────────────────────────
    mbe = float(np.mean(diff))
    mae = float(np.mean(np.abs(diff)))
    rmse = float(np.sqrt(np.mean(diff ** 2)))

    # R²
    ss_res = float(np.sum(diff ** 2))
    ss_tot = float(np.sum((yt - np.mean(yt)) ** 2))
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")

    # ── Correlations ─────────────────────────────────────────────────
    if n > 1:
        corr = float(np.corrcoef(yt, yp)[0, 1])
    else:
        corr = float("nan")

    # Spearman rank correlation (corr_rank in R)
    if n > 1:
        spearman_result = sp_stats.spearmanr(yt, yp)
        corr_rank = float(spearman_result.correlation)
    else:
        corr_rank = float("nan")

    # ── Power to select true top/bottom 10% ─────────────────────────
    n_top_bottom = max(1, int(round(0.1 * n)))

    top10_true_idx = np.argsort(-yt)[:n_top_bottom]
    bottom10_true_idx = np.argsort(yt)[:n_top_bottom]
    top10_pred_idx = np.argsort(-yp)[:n_top_bottom]
    bottom10_pred_idx = np.argsort(yp)[:n_top_bottom]

    power_t10 = float(
        np.isin(top10_true_idx, top10_pred_idx).sum() / n_top_bottom
    )
    power_b10 = float(
        np.isin(bottom10_true_idx, bottom10_pred_idx).sum() / n_top_bottom
    )

    # ── Variance components & heritability ───────────────────────────
    var_additive = float(np.var(yp, ddof=1))
    var_residual = float(np.var(diff, ddof=1))

    # Narrow-sense heritability scaled by prediction accuracy
    denom = var_additive + var_residual
    if denom > 0:
        h2 = round(var_additive / denom, 10)
        if h2 < 0.0 or h2 > 1.0:
            h2 = float("nan")
    else:
        h2 = float("nan")

    if verbose:
        print(f"Mean bias error (mbe)              = {mbe:.6f}")
        print(f"Mean absolute error (mae)            = {mae:.6f}")
        print(f"Root mean square error (rmse)        = {rmse:.6f}")
        print(f"Coefficient of determination (r2)    = {r2:.6f}")
        print(f"Pearson's correlation (corr)         = {corr:.6f}")
        print(f"Spearman's correlation (corr_rank)   = {corr_rank:.6f}")
        print(f"Power to identify top 10 (power_t10) = {power_t10:.6f}")
        print(f"Power to identify bottom 10 (power_b10)= {power_b10:.6f}")
        print(f"Variance of predicted (var_additive)  = {var_additive:.6f}")
        print(f"Variance of residuals (var_residual)  = {var_residual:.6f}")
        print(f"Narrow-sense heritability (h2)       = {h2}")

    return {
        "mbe": mbe,
        "mae": mae,
        "rmse": rmse,
        "r2": r2,
        "corr": corr,
        "corr_rank": corr_rank,
        "power_t10": power_t10,
        "power_b10": power_b10,
        "var_additive": var_additive,
        "var_residual": var_residual,
        "h2": h2,
    }

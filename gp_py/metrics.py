from __future__ import annotations

import numpy as np


def fn_prediction_performance_metrics(y_true, y_pred, *, verbose: bool = False) -> dict:
    yt = np.asarray(y_true, dtype=float)
    yp = np.asarray(y_pred, dtype=float)

    mask = np.isfinite(yt) & np.isfinite(yp)
    yt = yt[mask]
    yp = yp[mask]
    if yt.size == 0:
        raise ValueError("No valid y_true/y_pred pairs to score")

    diff = yt - yp
    mae = float(np.mean(np.abs(diff)))
    rmse = float(np.sqrt(np.mean(diff**2)))
    mbe = float(np.mean(diff))
    corr = float(np.corrcoef(yt, yp)[0, 1]) if yt.size > 1 else float("nan")

    ss_res = float(np.sum((yt - yp) ** 2))
    ss_tot = float(np.sum((yt - np.mean(yt)) ** 2))
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")

    return {
        "mbe": mbe,
        "mae": mae,
        "rmse": rmse,
        "r2": r2,
        "corr": corr,
    }

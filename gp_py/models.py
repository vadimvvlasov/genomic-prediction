from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.linear_model import BayesianRidge, ElasticNetCV, LassoCV, RidgeCV

from gp_py.metrics import fn_prediction_performance_metrics
from gp_py.schema import MergedData


def _split_xy(list_merged: MergedData, vec_idx_training, vec_idx_validation):
    X_train = list_merged.G.iloc[vec_idx_training, :].to_numpy(dtype=float)
    y_train = list_merged.y.iloc[vec_idx_training].to_numpy(dtype=float)
    X_valid = list_merged.G.iloc[vec_idx_validation, :].to_numpy(dtype=float)
    y_valid = list_merged.y.iloc[vec_idx_validation].to_numpy(dtype=float)

    train_mask = np.isfinite(y_train)
    return X_train[train_mask], y_train[train_mask], X_valid, y_valid


def _wrap_output(
    list_merged: MergedData, vec_idx_validation, y_validation, y_pred, coefs, intercept, model_name
):
    df_y_validation = pd.DataFrame(
        {
            "id": list_merged.G.index.to_series().iloc[vec_idx_validation].astype(str).to_list(),
            "pop": list_merged.pop.iloc[vec_idx_validation].astype(str).to_list(),
            "y_true": y_validation,
            "y_pred": y_pred,
        }
    )
    perf = fn_prediction_performance_metrics(df_y_validation["y_true"], df_y_validation["y_pred"])
    coef_names = ["intercept", *list_merged.G.columns.astype(str).to_list()]
    vec_effects = pd.Series([intercept, *coefs], index=coef_names, dtype=float)
    n_non_zero = int(np.sum(np.abs(vec_effects.to_numpy()) >= np.finfo(float).eps))
    return {
        "list_perf": perf,
        "df_y_validation": df_y_validation,
        "vec_effects": vec_effects,
        "n_non_zero": n_non_zero,
        "model": model_name,
    }


def fn_ridge(
    list_merged: MergedData,
    vec_idx_training,
    vec_idx_validation,
    other_params=None,
    *,
    verbose: bool = False,
) -> dict:
    X_train, y_train, X_valid, y_valid = _split_xy(
        list_merged, vec_idx_training, vec_idx_validation
    )
    model = RidgeCV(alphas=np.logspace(-6, 6, 25), fit_intercept=True)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_valid)
    return _wrap_output(
        list_merged, vec_idx_validation, y_valid, y_pred, model.coef_, model.intercept_, "ridge"
    )


def fn_lasso(
    list_merged: MergedData,
    vec_idx_training,
    vec_idx_validation,
    other_params=None,
    *,
    verbose: bool = False,
) -> dict:
    X_train, y_train, X_valid, y_valid = _split_xy(
        list_merged, vec_idx_training, vec_idx_validation
    )
    n_folds = 10 if other_params is None else int(other_params.get("n_folds", 10))
    model = LassoCV(cv=n_folds, random_state=123)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_valid)
    return _wrap_output(
        list_merged, vec_idx_validation, y_valid, y_pred, model.coef_, model.intercept_, "lasso"
    )


def fn_elastic_net(
    list_merged: MergedData,
    vec_idx_training,
    vec_idx_validation,
    other_params=None,
    *,
    verbose: bool = False,
) -> dict:
    X_train, y_train, X_valid, y_valid = _split_xy(
        list_merged, vec_idx_training, vec_idx_validation
    )
    model = ElasticNetCV(l1_ratio=[0.1, 0.5, 0.9], cv=10, random_state=123)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_valid)
    return _wrap_output(
        list_merged,
        vec_idx_validation,
        y_valid,
        y_pred,
        model.coef_,
        model.intercept_,
        "elastic_net",
    )


def _run_bglr_via_rpy2(model_name: str, X_train, y_train, X_valid):
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import numpy2ri
    except ImportError as exc:
        raise ImportError(
            "Install optional dependency group 'rbridge' for Bayes A/B/C bridge"
        ) from exc

    numpy2ri.activate()

    ro.r("library(BGLR)")
    ro.globalenv["X_train"] = X_train
    ro.globalenv["y_train"] = y_train
    ro.globalenv["X_valid"] = X_valid
    ro.globalenv["model_name"] = model_name
    ro.r(
        """
        yNA <- c(y_train, rep(NA_real_, nrow(X_valid)))
        Xall <- rbind(X_train, X_valid)
        ETA <- list(list(X=Xall, model=model_name))
        fit <- BGLR::BGLR(y=yNA, ETA=ETA, nIter=12000, burnIn=2000, verbose=FALSE)
        y_pred_bridge <- fit$yHat[(length(y_train)+1):length(yNA)]
        b_bridge <- fit$ETA[[1]]$b
        """
    )
    y_pred = np.array(ro.globalenv["y_pred_bridge"], dtype=float)
    coefs = np.array(ro.globalenv["b_bridge"], dtype=float)
    return y_pred, coefs, 0.0


def _run_bayes_native(X_train, y_train, X_valid):
    model = BayesianRidge(compute_score=False)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_valid)
    return y_pred, model.coef_, float(model.intercept_)


def _fn_bayes_bridge(
    list_merged: MergedData,
    vec_idx_training,
    vec_idx_validation,
    *,
    model_name: str,
    backend: str = "auto",
):
    X_train, y_train, X_valid, y_valid = _split_xy(
        list_merged, vec_idx_training, vec_idx_validation
    )
    bglr_map = {"Bayes_A": "BayesA", "Bayes_B": "BayesB", "Bayes_C": "BayesC"}
    if backend not in {"auto", "native", "rbridge"}:
        raise ValueError("bayes_backend must be one of: auto, native, rbridge")

    if backend == "native":
        y_pred, coefs, intercept = _run_bayes_native(X_train, y_train, X_valid)
    elif backend == "rbridge":
        y_pred, coefs, intercept = _run_bglr_via_rpy2(
            bglr_map[model_name], X_train, y_train, X_valid
        )
    else:
        try:
            y_pred, coefs, intercept = _run_bglr_via_rpy2(
                bglr_map[model_name], X_train, y_train, X_valid
            )
        except Exception:
            y_pred, coefs, intercept = _run_bayes_native(X_train, y_train, X_valid)
    return _wrap_output(
        list_merged, vec_idx_validation, y_valid, y_pred, coefs, intercept, model_name
    )


def fn_Bayes_A(
    list_merged: MergedData,
    vec_idx_training,
    vec_idx_validation,
    other_params=None,
    *,
    verbose: bool = False,
) -> dict:
    backend = "auto" if other_params is None else str(other_params.get("bayes_backend", "auto"))
    return _fn_bayes_bridge(
        list_merged,
        vec_idx_training,
        vec_idx_validation,
        model_name="Bayes_A",
        backend=backend,
    )


def fn_Bayes_B(
    list_merged: MergedData,
    vec_idx_training,
    vec_idx_validation,
    other_params=None,
    *,
    verbose: bool = False,
) -> dict:
    backend = "auto" if other_params is None else str(other_params.get("bayes_backend", "auto"))
    return _fn_bayes_bridge(
        list_merged,
        vec_idx_training,
        vec_idx_validation,
        model_name="Bayes_B",
        backend=backend,
    )


def fn_Bayes_C(
    list_merged: MergedData,
    vec_idx_training,
    vec_idx_validation,
    other_params=None,
    *,
    verbose: bool = False,
) -> dict:
    backend = "auto" if other_params is None else str(other_params.get("bayes_backend", "auto"))
    return _fn_bayes_bridge(
        list_merged,
        vec_idx_training,
        vec_idx_validation,
        model_name="Bayes_C",
        backend=backend,
    )


def fn_gBLUP(
    list_merged: MergedData,
    vec_idx_training,
    vec_idx_validation,
    other_params=None,
    *,
    verbose: bool = False,
) -> dict:
    raise NotImplementedError("gBLUP scaffold only. Implement GRM + mixed model next.")

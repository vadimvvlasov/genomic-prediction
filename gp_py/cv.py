from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold

from gp_py.models import (
    fn_Bayes_A,
    fn_Bayes_B,
    fn_Bayes_C,
    fn_elastic_net,
    fn_gBLUP,
    fn_lasso,
    fn_LightGBM,
    fn_ridge,
    fn_RandomForest,
    fn_SVR,
    fn_XGBoost,
)
from gp_py.schema import MergedData


MODEL_REGISTRY = {
    "ridge": fn_ridge,
    "lasso": fn_lasso,
    "elastic_net": fn_elastic_net,
    "Bayes_A": fn_Bayes_A,
    "Bayes_B": fn_Bayes_B,
    "Bayes_C": fn_Bayes_C,
    "gBLUP": fn_gBLUP,
    # Machine learning models
    "SVR": fn_SVR,
    "RandomForest": fn_RandomForest,
    "XGBoost": fn_XGBoost,
    "LightGBM": fn_LightGBM,
}


@dataclass
class CVRun:
    model: str
    rep: int
    fold: int
    train_idx: np.ndarray
    valid_idx: np.ndarray


def fn_cross_validation_preparation(
    list_merged: MergedData,
    *,
    cv_type: int = 1,
    n_folds: int = 10,
    n_reps: int = 10,
    vec_models_to_test: tuple[str, ...] = ("ridge", "lasso", "elastic_net"),
    max_mem_gb: float = 15.0,
    verbose: bool = False,
) -> list[CVRun]:
    if cv_type != 1:
        raise NotImplementedError(
            "Scaffold currently implements only cv_type=1 (within-pop k-fold)"
        )

    n = list_merged.G.shape[0]
    if n < 4:
        raise ValueError("Need at least 4 samples for k-fold CV scaffold")

    runs: list[CVRun] = []
    for rep in range(1, n_reps + 1):
        kf = KFold(n_splits=n_folds, shuffle=True, random_state=1000 + rep)
        for fold_idx, (train_idx, valid_idx) in enumerate(kf.split(np.arange(n)), start=1):
            for model in vec_models_to_test:
                runs.append(
                    CVRun(
                        model=model,
                        rep=rep,
                        fold=fold_idx,
                        train_idx=train_idx,
                        valid_idx=valid_idx,
                    )
                )
    return runs


def _run_single_cv(
    run: CVRun,
    list_merged: MergedData,
    bayes_backend: str,
    verbose: bool,
) -> tuple[dict, pd.DataFrame]:
    """Worker function for a single CV fold/rep/model combination."""
    fn_model = MODEL_REGISTRY.get(run.model)
    if fn_model is None:
        return {}, pd.DataFrame()
    out = fn_model(
        list_merged,
        run.train_idx,
        run.valid_idx,
        other_params={"n_folds": 10, "bayes_backend": bayes_backend},
        verbose=verbose,
    )
    perf = out["list_perf"]
    metrics_row = {
        "rep": run.rep,
        "fold": run.fold,
        "model": run.model,
        **perf,
    }
    dfv = out["df_y_validation"].copy()
    dfv["rep"] = run.rep
    dfv["fold"] = run.fold
    dfv["model"] = run.model
    return metrics_row, dfv


def fn_cross_validation_within_population(
    list_merged: MergedData,
    *,
    n_folds: int = 10,
    n_reps: int = 10,
    vec_models_to_test: tuple[str, ...] = ("ridge", "lasso", "elastic_net"),
    bool_parallel: bool = True,
    bayes_backend: str = "auto",
    max_mem_gb: float = 15.0,
    n_threads: int = 2,
    dir_output: str | None = None,
    verbose: bool = False,
) -> dict:
    runs = fn_cross_validation_preparation(
        list_merged,
        cv_type=1,
        n_folds=n_folds,
        n_reps=n_reps,
        vec_models_to_test=vec_models_to_test,
        max_mem_gb=max_mem_gb,
        verbose=verbose,
    )

    if bool_parallel and n_threads > 1:
        from joblib import Parallel, delayed

        results = Parallel(n_jobs=n_threads, verbose=max(0, verbose * 10))(
            delayed(_run_single_cv)(run, list_merged, bayes_backend, verbose)
            for run in runs
        )
    else:
        try:
            from tqdm import tqdm
        except ImportError:
            tqdm = iter  # fallback: no progress bar

        results = [
            _run_single_cv(run, list_merged, bayes_backend, verbose)
            for run in tqdm(runs, desc="CV folds", disable=not verbose)
        ]

    metrics_rows = []
    ypred_rows = []
    for metrics_row, dfv in results:
        if metrics_row:
            metrics_rows.append(metrics_row)
            ypred_rows.append(dfv)

    metrics = pd.DataFrame(metrics_rows)
    ypred = pd.concat(ypred_rows, ignore_index=True) if ypred_rows else pd.DataFrame()
    return {
        "METRICS_WITHIN_POP": metrics,
        "YPRED_WITHIN_POP": ypred,
    }

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNetCV, LassoCV, RidgeCV

from gp_py.bayes import bayes_a_gibbs, bayes_b_gibbs, bayes_c_gibbs
from gp_py.distances import (
    add_diagonal_load,
    genomic_relationship_matrix,
)
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

    yt = df_y_validation["y_true"]
    yp = df_y_validation["y_pred"]
    # Filter out NaN pairs (can happen if MCMC diverges)
    valid_mask = np.isfinite(yt) & np.isfinite(yp)
    has_valid_true = np.any(valid_mask)

    if has_valid_true:
        perf = fn_prediction_performance_metrics(
            yt[valid_mask], yp[valid_mask]
        )
    else:
        perf = {
            "mbe": float("nan"),
            "mae": float("nan"),
            "rmse": float("nan"),
            "r2": float("nan"),
            "corr": float("nan"),
        }

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
    model = ElasticNetCV(
        l1_ratio=[0.1, 0.5, 0.9], cv=10, random_state=123,
        max_iter=100_000, tol=1e-4,
    )
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


def _run_bayes_a_native(X_train, y_train, X_valid, *, seed: int = 42):
    """BayesA via MCMC Gibbs sampler (no R dependency)."""
    result = bayes_a_gibbs(X_train, y_train, n_iter=12000, burn_in=2000, thin=10, seed=seed)
    beta = result["beta_mean"]
    intercept = result["intercept"]
    y_pred = X_valid @ beta + intercept
    return y_pred, beta, intercept


def _run_bayes_b_native(X_train, y_train, X_valid, *, seed: int = 42):
    """BayesB via MCMC Gibbs sampler with spike-and-slab prior."""
    result = bayes_b_gibbs(X_train, y_train, n_iter=12000, burn_in=2000, thin=10, seed=seed)
    beta = result["beta_mean"]
    intercept = result["intercept"]
    y_pred = X_valid @ beta + intercept
    return y_pred, beta, intercept


def _run_bayes_c_native(X_train, y_train, X_valid, *, seed: int = 42):
    """BayesC via MCMC Gibbs sampler with common marker variance."""
    result = bayes_c_gibbs(X_train, y_train, n_iter=12000, burn_in=2000, thin=10, seed=seed)
    beta = result["beta_mean"]
    intercept = result["intercept"]
    y_pred = X_valid @ beta + intercept
    return y_pred, beta, intercept


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

    # Dispatch to correct native MCMC sampler
    native_dispatch = {
        "Bayes_A": _run_bayes_a_native,
        "Bayes_B": _run_bayes_b_native,
        "Bayes_C": _run_bayes_c_native,
    }

    if backend == "native":
        y_pred, coefs, intercept = native_dispatch[model_name](X_train, y_train, X_valid)
    elif backend == "rbridge":
        y_pred, coefs, intercept = _run_bglr_via_rpy2(
            bglr_map[model_name], X_train, y_train, X_valid
        )
    else:
        # auto: try R bridge first, fall back to native MCMC
        try:
            y_pred, coefs, intercept = _run_bglr_via_rpy2(
                bglr_map[model_name], X_train, y_train, X_valid
            )
        except Exception:
            y_pred, coefs, intercept = native_dispatch[model_name](X_train, y_train, X_valid)
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
    """Genomic BLUP via mixed linear model.

    Mirrors the R implementation (R/models.R:fn_gBLUP):
    1. Build GRM (genomic relationship matrix) from training genotypes.
    2. Fit mixed model: y = Xβ + Zu + ε, where
       - X: fixed effects (intercept + optional covariates)
       - Z: identity for random genetic effects (one per sample)
       - u ~ N(0, σ²_g * GRM): random additive genetic effects
       - ε ~ N(0, σ²_e * I): residual error
    3. Predict BLUPs for validation samples via relationship-based kriging.
    """
    try:
        import statsmodels.api as sm
    except ImportError as exc:
        raise ImportError(
            "Install 'statsmodels' for gBLUP: uv add statsmodels"
        ) from exc

    # ── Prepare data ──────────────────────────────────────────────
    y_full = list_merged.y.copy()
    y_training = y_full.iloc[vec_idx_training].to_numpy(dtype=float)
    y_validation = y_full.iloc[vec_idx_validation].to_numpy(dtype=float)

    # Remove missing phenotypes from training
    train_mask = np.isfinite(y_training)
    idx_clean = np.array(vec_idx_training)[train_mask]
    y_train = y_training[train_mask]
    G_train = list_merged.G.iloc[idx_clean, :].copy()

    n_train_clean = len(idx_clean)
    ids_train = G_train.index.tolist()

    # ── Build GRM ─────────────────────────────────────────────────
    K = genomic_relationship_matrix(G_train)

    # Ensure positive-definiteness via diagonal load if needed
    try:
        np.linalg.cholesky(K.values)
    except np.linalg.LinAlgError:
        K = add_diagonal_load(K, eps=0.001)
        if verbose:
            print("[gBLUP] GRM was not positive-definite; added diagonal load.")

    # ── Build mixed model: y = Xβ + Zu + ε ────────────────────────
    # Fixed effects: intercept + covariates (if present)
    X_fixed = np.ones((n_train_clean, 1), dtype=float)
    covar_train = None
    if list_merged.covar is not None:
        covar_raw = list_merged.covar.iloc[idx_clean, :].to_numpy(dtype=float)
        # If >10 covariates, use first 10 PCs (mirrors R logic)
        if covar_raw.shape[1] > 10:
            from sklearn.decomposition import PCA

            covar_scaled = (covar_raw - covar_raw.mean(axis=0)) / (
                covar_raw.std(axis=0) + 1e-10
            )
            pca = PCA(n_components=10)
            covar_train = pca.fit_transform(covar_scaled)
        else:
            covar_train = covar_raw
        X_fixed = np.hstack([X_fixed, covar_train])

    b_names = ["intercept"]
    if covar_train is not None:
        b_names += [f"covar_{i}" for i in range(covar_train.shape[1])]

    try:
        # exog_vc: pass Cholesky factor of GRM so that Var(Zu) = L @ L' = GRM
        L = np.linalg.cholesky(K.values)
        # statsmodels 0.14+ expects {group_name: {vc_name: array}}
        group_name = "genetic"
        exog_vc = {group_name: {"grm": L}}

        md = sm.MixedLM(
            endog=y_train,
            exog=X_fixed,
            exog_vc=exog_vc,
            groups=[group_name] * n_train_clean,
        )
        mdf = md.fit(method="bfgs", maxiter=200, disp=verbose)

        # Extract fixed effects (numpy array in statsmodels 0.14+)
        b_hat = np.asarray(mdf.fe_params, dtype=float)

        # Extract random effects (BLUPs)
        re = mdf.random_effects
        if isinstance(re, dict):
            u_hat = np.asarray(re[group_name], dtype=float).flatten()
        else:
            u_hat = np.asarray(re, dtype=float).flatten()

        # ── Predict for validation samples ────────────────────────
        # Build GRM including validation samples for relationship-based prediction
        G_val = list_merged.G.iloc[vec_idx_validation, :].copy()
        G_all = pd.concat([G_train, G_val])
        K_all = genomic_relationship_matrix(G_all)

        try:
            np.linalg.cholesky(K_all.values)
        except np.linalg.LinAlgError:
            K_all = add_diagonal_load(K_all, eps=0.001)

        n_val = len(vec_idx_validation)

        # Fixed effect prediction
        X_val = np.ones((n_val, 1))
        if covar_train is not None and list_merged.covar is not None:
            covar_val_raw = list_merged.covar.iloc[
                vec_idx_validation, :
            ].to_numpy(dtype=float)
            if covar_train.shape[1] <= 10:
                covar_val = covar_val_raw
            else:
                from sklearn.decomposition import PCA

                # Re-fit PCA on all data for consistent projection
                covar_all_raw = np.vstack([covar_raw, covar_val_raw])
                covar_all_scaled = (covar_all_raw - covar_all_raw.mean(axis=0)) / (
                    covar_all_raw.std(axis=0) + 1e-10
                )
                pca_all = PCA(n_components=10)
                covar_all_pca = pca_all.fit_transform(covar_all_scaled)
                covar_val = covar_all_pca[n_train_clean:]
            X_val = np.hstack([X_val, covar_val])

        y_pred_fixed = X_val @ b_hat

        # Random effect prediction via kriging: u_v = K_vt @ K_tt^{-1} @ u_train
        K_vt = K_all.iloc[n_train_clean:, :n_train_clean].values
        try:
            L_tt = np.linalg.cholesky(K_all.iloc[:n_train_clean, :n_train_clean].values)
            alpha = np.linalg.solve(L_tt, u_hat)
            u_val = K_vt @ np.linalg.solve(L_tt.T, alpha)
        except np.linalg.LinAlgError:
            if verbose:
                print("[gBLUP] Cholesky failed for validation; using mean BLUP.")
            u_val = np.mean(u_hat) * np.ones(n_val)

        y_pred = y_pred_fixed + u_val
        coefs_combined = np.concatenate([b_hat, u_hat])
        coef_names = b_names + ids_train

    except Exception as e:
        if verbose:
            print(f"[gBLUP] Model fitting failed: {e}. Falling back to mean prediction.")
        y_pred = np.full(len(vec_idx_validation), np.mean(y_train))
        coefs_combined = np.array([np.mean(y_train)])
        coef_names = ["intercept"]

    # ── Wrap output ───────────────────────────────────────────────
    vec_effects = pd.Series(coefs_combined, index=coef_names, dtype=float)
    n_non_zero = int(np.sum(np.abs(vec_effects.to_numpy()) >= np.finfo(float).eps))

    df_y_validation = pd.DataFrame(
        {
            "id": list_merged.G.index.to_series()
            .iloc[vec_idx_validation]
            .astype(str)
            .to_list(),
            "pop": list_merged.pop.iloc[vec_idx_validation].astype(str).to_list(),
            "y_true": y_validation,
            "y_pred": y_pred,
        }
    )

    yt = df_y_validation["y_true"]
    has_valid_true = np.any(np.isfinite(yt))

    if has_valid_true:
        perf = fn_prediction_performance_metrics(
            df_y_validation["y_true"], df_y_validation["y_pred"]
        )
    else:
        perf = {
            "mbe": float("nan"),
            "mae": float("nan"),
            "rmse": float("nan"),
            "r2": float("nan"),
            "corr": float("nan"),
        }

    return {
        "list_perf": perf,
        "df_y_validation": df_y_validation,
        "vec_effects": vec_effects,
        "n_non_zero": n_non_zero,
        "model": "gBLUP",
    }

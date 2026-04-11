"""
Optimized native MCMC Gibbs samplers for Bayes A, Bayes B, and Bayes C.

Key optimizations:
1. Precompute Gram matrix G = X'X (once, O(np²))
2. Standardize X internally for numerical stability
3. Vectorized posterior calculations
4. Numerical safety: clamped variance, NaN/Inf detection, divergence guards

All three models share the common structure:
    y = Xβ + ε,   ε ~ N(0, σ²_e I)

but differ in their priors on β:
    BayesA: β_j ~ N(0, σ²_j), σ²_j ~ S⁻²χ²(dfa, S²a)
    BayesB: β_j = 0 w.p. (1-π); β_j ~ N(0, σ²_j) w.p. π, σ²_j ~ S⁻²χ²(dfa, S²a)
    BayesC: β_j = 0 w.p. (1-π); β_j ~ N(0, σ²) w.p. π (common σ² for all non-zero)
"""

from __future__ import annotations

import numpy as np
from numpy.random import Generator

# ── Numerical safety constants ──────────────────────────────────────
_SIGMA2_E_MIN_FRAC = 1e-4    # σ²_e minimum = this fraction of var(y)
_SIGMA2_E_MAX_FRAC = 10.0    # σ²_e maximum = this multiple of var(y)
_VAR_POST_MAX = 1e10         # cap on posterior variance
_MU_POST_MAX = 1e6           # cap on posterior mean magnitude
_BETA_MAX = 1e6              # cap on individual β
_RSS_MAX = 1e20              # cap on residual sum of squares
_SIGMA2_MAX_FRAC = 100.0     # common marker variance cap (multiple of var(y))


def _safe(val: float, default: float = 0.0) -> float:
    """Replace NaN/Inf with default."""
    if not np.isfinite(val):
        return default
    return val


def _clamp(val: float, lo: float, hi: float) -> float:
    """Clamp scalar to [lo, hi], returning lo if NaN/Inf."""
    if not np.isfinite(val):
        return lo
    return max(lo, min(hi, val))


def _sample_scaled_inv_chisq(rng: Generator, df: float, scale: float) -> float:
    """Sample from scaled inverse chi-square: X ~ S⁻²χ²(df, scale)."""
    if df <= 0:
        return max(scale, 1e-300)
    y = rng.chisquare(df)
    result = df * max(scale, 1e-300) / max(y, 1e-300)
    # Guard against extreme overflow from chisquare ≈ 0
    return min(result, 1e300)


def _standardize(X: np.ndarray, y: np.ndarray):
    """Standardize X columns and center y. Returns (X_std, y_c, means, sds, y_mean)."""
    X_mean = X.mean(axis=0)
    X_sd = X.std(axis=0)
    # Prevent division by near-zero
    X_sd = np.where(X_sd < 1e-10, 1.0, X_sd)
    X_std = (X - X_mean) / X_sd
    y_mean = float(np.mean(y))
    y_c = y - y_mean
    return X_std, y_c, X_mean, X_sd, y_mean


def _back_transform_beta(beta_std: np.ndarray, X_sd: np.ndarray, y_sd: float,
                         y_mean: float = 0.0, X_mean: np.ndarray | None = None) -> tuple[np.ndarray, float]:
    """Transform standardized-space β back to original scale."""
    beta_orig = beta_std * y_sd / X_sd
    intercept = 0.0
    if X_mean is not None:
        intercept = y_mean - X_mean @ beta_orig
        intercept = _safe(intercept, y_mean)
    return beta_orig, float(intercept)


def _init_priors(var_y: float, p: int, r2: float = 0.5, df_a: int = 5) -> dict:
    """Compute prior hyperparameters."""
    if var_y < 1e-10:
        var_y = 1.0
    df_e = 5.0
    s2_e = var_y * (1.0 - r2)
    s2_a = r2 * var_y * max(df_a - 2, 1) / (p * df_a) if (p > 0 and df_a > 2) else 1.0
    return {
        "var_y": var_y, "df_e": df_e, "s2_e": max(s2_e, 1e-10),
        "df_a": df_a, "s2_a": max(s2_a, 1e-10),
    }


def _safe_rss(y_c: np.ndarray, X: np.ndarray, beta: np.ndarray,
              rss_max: float = _RSS_MAX) -> float:
    """Compute RSS = ||y - Xβ||² with overflow protection."""
    resid = y_c - X @ beta
    if not np.all(np.isfinite(resid)):
        return rss_max
    rss = float(resid @ resid)
    return min(rss, rss_max)


# ── Bayes A ─────────────────────────────────────────────────────────

def bayes_a_gibbs(
    X: np.ndarray, y: np.ndarray,
    n_iter: int = 12000, burn_in: int = 2000, thin: int = 10,
    seed: int = 42, r2: float = 0.5, df_a: int = 5,
) -> dict:
    """BayesA MCMC Gibbs sampler (vectorized, standardized, numerically stable)."""
    rng = np.random.default_rng(seed)
    n, p = X.shape

    X_std, y_c, X_mean, X_sd, y_mean = _standardize(X, y)
    y_var = float(np.var(y_c, ddof=1))
    priors = _init_priors(y_var, p, r2=r2, df_a=df_a)

    # Minimum σ²_e: don't let it collapse below fraction of phenotypic variance
    sigma2_e_min = y_var * _SIGMA2_E_MIN_FRAC
    sigma2_e_max = y_var * _SIGMA2_E_MAX_FRAC
    sigma2_j_max = y_var * _SIGMA2_MAX_FRAC

    G = X_std.T @ X_std
    d = np.diag(G).copy()
    Xty = X_std.T @ y_c

    sigma2_e = priors["var_y"]
    sigma2_j = np.full(p, priors["s2_a"])
    beta = np.zeros(p)

    n_saved = max(1, (n_iter - burn_in) // thin)
    beta_samples = np.zeros((n_saved, p))
    sigma2_e_samples = np.zeros(n_saved)
    idx = 0
    diverged = False

    for it in range(n_iter):
        # Vectorized sweep: Xtr = X'y - X'Xβ + d*β_j
        G_beta = G @ beta
        if not np.all(np.isfinite(G_beta)):
            # Reset beta if it exploded
            beta = np.zeros(p)
            sigma2_e = priors["s2_e"]
            sigma2_j = np.full(p, priors["s2_a"])
            continue

        Xtr = Xty - G_beta + d * beta

        prec = d / sigma2_e + 1.0 / sigma2_j
        var_post = 1.0 / np.maximum(prec, 1e-300)
        var_post = np.minimum(var_post, _VAR_POST_MAX)

        mu_post = var_post * Xtr / sigma2_e
        mu_post = np.clip(mu_post, -_MU_POST_MAX, _MU_POST_MAX)

        beta = mu_post + np.sqrt(np.maximum(var_post, 1e-300)) * rng.normal(size=p)
        beta = np.clip(beta, -_BETA_MAX, _BETA_MAX)

        # Update marker variances
        s2_post = (priors["df_a"] * priors["s2_a"] + beta ** 2) / (priors["df_a"] + 1)
        for j in range(p):
            sigma2_j[j] = _sample_scaled_inv_chisq(rng, priors["df_a"] + 1, max(s2_post[j], 1e-300))
            sigma2_j[j] = max(sigma2_j[j], 1e-300)
            sigma2_j[j] = min(sigma2_j[j], sigma2_j_max)

        # Update residual variance
        rss = _safe_rss(y_c, X_std, beta)
        sigma2_e = _sample_scaled_inv_chisq(
            rng, priors["df_e"] + n,
            (priors["df_e"] * priors["s2_e"] + rss) / (priors["df_e"] + n),
        )
        # Clamp σ²_e to prevent collapse and explosion
        sigma2_e = max(sigma2_e, sigma2_e_min)
        sigma2_e = min(sigma2_e, sigma2_e_max)

        if it >= burn_in and (it - burn_in) % thin == 0:
            if idx < n_saved:
                beta_samples[idx] = beta.copy()
                sigma2_e_samples[idx] = sigma2_e
                idx += 1

        # Divergence check: if σ²_e hit the floor, the chain is unstable
        if sigma2_e <= sigma2_e_min * 1.01:
            diverged = True

    # Back-transform to original scale
    beta_mean = beta_samples.mean(axis=0)
    if not np.all(np.isfinite(beta_mean)):
        beta_mean = np.where(np.isfinite(beta_mean), beta_mean, 0.0)
    beta_orig, intercept = _back_transform_beta(beta_mean, X_sd, 1.0, y_mean, X_mean)

    return {
        "beta_mean": beta_orig,
        "intercept": intercept,
        "sigma2_e_mean": float(sigma2_e_samples.mean()),
        "sigma2_j_mean": sigma2_j,
        "samples_beta": beta_samples,
        "n_saved": idx,
        "diverged": diverged,
    }


# ── Bayes B ─────────────────────────────────────────────────────────

def bayes_b_gibbs(
    X: np.ndarray, y: np.ndarray,
    n_iter: int = 12000, burn_in: int = 2000, thin: int = 10,
    seed: int = 42, r2: float = 0.5, df_a: int = 5,
    pi_prior_a: float = 1.0, pi_prior_b: float = 1.0,
) -> dict:
    """BayesB MCMC Gibbs sampler with spike-and-slab prior (numerically stable)."""
    rng = np.random.default_rng(seed)
    n, p = X.shape

    X_std, y_c, X_mean, X_sd, y_mean = _standardize(X, y)
    y_var = float(np.var(y_c, ddof=1))
    priors = _init_priors(y_var, p, r2=r2, df_a=df_a)

    sigma2_e_min = y_var * _SIGMA2_E_MIN_FRAC
    sigma2_e_max = y_var * _SIGMA2_E_MAX_FRAC
    sigma2_j_max = y_var * _SIGMA2_MAX_FRAC

    G = X_std.T @ X_std
    d = np.diag(G).copy()
    Xty = X_std.T @ y_c

    sigma2_e = priors["var_y"]
    sigma2_j = np.full(p, priors["s2_a"])
    beta = np.zeros(p)
    delta = np.zeros(p, dtype=int)
    pi = 0.5

    n_saved = max(1, (n_iter - burn_in) // thin)
    beta_samples = np.zeros((n_saved, p))
    sigma2_e_samples = np.zeros(n_saved)
    pi_samples = np.zeros(n_saved)
    delta_samples = np.zeros((n_saved, p))
    idx = 0
    diverged = False

    for it in range(n_iter):
        G_beta = G @ beta
        if not np.all(np.isfinite(G_beta)):
            beta = np.zeros(p)
            sigma2_e = priors["s2_e"]
            sigma2_j = np.full(p, priors["s2_a"])
            delta = np.zeros(p, dtype=int)
            pi = 0.5
            continue

        Xtr = Xty - G_beta + d * beta

        n_included = 0
        for j in range(p):
            if d[j] < 1e-300:
                continue

            v_j = max(sigma2_j[j], 1e-300)
            prec_j = d[j] / sigma2_e + 1.0 / v_j
            var_post_j = 1.0 / max(prec_j, 1e-300)
            var_post_j = min(var_post_j, _VAR_POST_MAX)

            mu_post_j = _safe(var_post_j * Xtr[j] / sigma2_e, 0.0)
            mu_post_j = _clamp(mu_post_j, -_MU_POST_MAX, _MU_POST_MAX)

            mu_sq_2v = mu_post_j ** 2 / (2.0 * max(var_post_j, 1e-300))

            w_j = v_j * d[j] / sigma2_e
            w_j = min(w_j, 1e300)
            log_bf = -0.5 * np.log1p(w_j) + mu_sq_2v
            log_bf = _clamp(log_bf, -500, 500)

            log_odds = (np.log(max(pi, 1e-300))
                        - np.log(max(1.0 - pi, 1e-300))
                        + log_bf)
            prob_include = 1.0 / (1.0 + np.exp(-np.clip(log_odds, -500, 500)))

            if rng.random() < prob_include:
                delta[j] = 1
                beta[j] = rng.normal(mu_post_j, np.sqrt(max(var_post_j, 1e-300)))
                beta[j] = _clamp(beta[j], -_BETA_MAX, _BETA_MAX)
                n_included += 1

                s2_post = (priors["df_a"] * priors["s2_a"] + beta[j] ** 2) / (priors["df_a"] + 1)
                sigma2_j[j] = _sample_scaled_inv_chisq(rng, priors["df_a"] + 1, max(s2_post, 1e-300))
                sigma2_j[j] = max(sigma2_j[j], 1e-300)
                sigma2_j[j] = min(sigma2_j[j], sigma2_j_max)
            else:
                delta[j] = 0
                beta[j] = 0.0
                sigma2_j[j] = _sample_scaled_inv_chisq(rng, priors["df_a"], priors["s2_a"])

        pi = float(rng.beta(pi_prior_a + n_included, pi_prior_b + p - n_included))

        rss = _safe_rss(y_c, X_std, beta)
        sigma2_e = _sample_scaled_inv_chisq(
            rng, priors["df_e"] + n,
            (priors["df_e"] * priors["s2_e"] + rss) / (priors["df_e"] + n),
        )
        sigma2_e = max(sigma2_e, sigma2_e_min)
        sigma2_e = min(sigma2_e, sigma2_e_max)

        if it >= burn_in and (it - burn_in) % thin == 0:
            if idx < n_saved:
                beta_samples[idx] = beta.copy()
                sigma2_e_samples[idx] = sigma2_e
                pi_samples[idx] = pi
                delta_samples[idx] = delta.copy()
                idx += 1

        if sigma2_e <= sigma2_e_min * 1.01:
            diverged = True

    beta_mean = beta_samples.mean(axis=0)
    if not np.all(np.isfinite(beta_mean)):
        beta_mean = np.where(np.isfinite(beta_mean), beta_mean, 0.0)
    beta_orig, intercept = _back_transform_beta(beta_mean, X_sd, 1.0, y_mean, X_mean)

    return {
        "beta_mean": beta_orig,
        "intercept": intercept,
        "sigma2_e_mean": float(sigma2_e_samples.mean()),
        "pi_mean": float(pi_samples.mean()),
        "delta_mean": delta_samples.mean(axis=0),
        "samples_beta": beta_samples,
        "n_saved": idx,
        "diverged": diverged,
    }


# ── Bayes C ─────────────────────────────────────────────────────────

def bayes_c_gibbs(
    X: np.ndarray, y: np.ndarray,
    n_iter: int = 12000, burn_in: int = 2000, thin: int = 10,
    seed: int = 42, r2: float = 0.5, df_a: int = 5,
    pi_prior_a: float = 1.0, pi_prior_b: float = 1.0,
) -> dict:
    """BayesC MCMC Gibbs sampler with spike-and-slab + common variance (numerically stable)."""
    rng = np.random.default_rng(seed)
    n, p = X.shape

    X_std, y_c, X_mean, X_sd, y_mean = _standardize(X, y)
    y_var = float(np.var(y_c, ddof=1))
    priors = _init_priors(y_var, p, r2=r2, df_a=df_a)

    sigma2_e_min = y_var * _SIGMA2_E_MIN_FRAC
    sigma2_e_max = y_var * _SIGMA2_E_MAX_FRAC
    sigma2_max = y_var * _SIGMA2_MAX_FRAC

    G = X_std.T @ X_std
    d = np.diag(G).copy()
    Xty = X_std.T @ y_c

    sigma2_e = priors["var_y"]
    sigma2 = priors["s2_a"]
    beta = np.zeros(p)
    delta = np.zeros(p, dtype=int)
    pi = 0.5

    n_saved = max(1, (n_iter - burn_in) // thin)
    beta_samples = np.zeros((n_saved, p))
    sigma2_e_samples = np.zeros(n_saved)
    sigma2_samples = np.zeros(n_saved)
    pi_samples = np.zeros(n_saved)
    delta_samples = np.zeros((n_saved, p))
    idx = 0
    diverged = False

    for it in range(n_iter):
        G_beta = G @ beta
        if not np.all(np.isfinite(G_beta)):
            beta = np.zeros(p)
            sigma2_e = priors["s2_e"]
            sigma2 = priors["s2_a"] * p
            delta = np.zeros(p, dtype=int)
            pi = 0.5
            continue

        Xtr = Xty - G_beta + d * beta

        n_included = 0
        sum_beta_sq = 0.0
        for j in range(p):
            if d[j] < 1e-300:
                continue

            sv = max(sigma2, 1e-300)
            prec_j = d[j] / sigma2_e + 1.0 / sv
            var_post_j = 1.0 / max(prec_j, 1e-300)
            var_post_j = min(var_post_j, _VAR_POST_MAX)

            mu_post_j = _safe(var_post_j * Xtr[j] / sigma2_e, 0.0)
            mu_post_j = _clamp(mu_post_j, -_MU_POST_MAX, _MU_POST_MAX)

            mu_sq_2v = mu_post_j ** 2 / (2.0 * max(var_post_j, 1e-300))

            w_j = sv * d[j] / sigma2_e
            w_j = min(w_j, 1e300)
            log_bf = -0.5 * np.log1p(w_j) + mu_sq_2v
            log_bf = _clamp(log_bf, -500, 500)

            log_odds = (np.log(max(pi, 1e-300))
                        - np.log(max(1.0 - pi, 1e-300))
                        + log_bf)
            prob_include = 1.0 / (1.0 + np.exp(-np.clip(log_odds, -500, 500)))

            if rng.random() < prob_include:
                delta[j] = 1
                beta[j] = rng.normal(mu_post_j, np.sqrt(max(var_post_j, 1e-300)))
                beta[j] = _clamp(beta[j], -_BETA_MAX, _BETA_MAX)
                sum_beta_sq += beta[j] ** 2
                n_included += 1
            else:
                delta[j] = 0
                beta[j] = 0.0

        if n_included > 0:
            df_post = df_a + n_included
            s2_post = (df_a * priors["s2_a"] + sum_beta_sq) / df_post
            sigma2 = _sample_scaled_inv_chisq(rng, df_post, max(s2_post, 1e-300))
        else:
            sigma2 = _sample_scaled_inv_chisq(rng, df_a, priors["s2_a"])

        # Clamp common variance to prevent divergence
        sigma2 = max(sigma2, 1e-300)
        sigma2 = min(sigma2, sigma2_max)

        pi = float(rng.beta(pi_prior_a + n_included, pi_prior_b + p - n_included))

        rss = _safe_rss(y_c, X_std, beta)
        sigma2_e = _sample_scaled_inv_chisq(
            rng, priors["df_e"] + n,
            (priors["df_e"] * priors["s2_e"] + rss) / (priors["df_e"] + n),
        )
        sigma2_e = max(sigma2_e, sigma2_e_min)
        sigma2_e = min(sigma2_e, sigma2_e_max)

        if it >= burn_in and (it - burn_in) % thin == 0:
            if idx < n_saved:
                beta_samples[idx] = beta.copy()
                sigma2_e_samples[idx] = sigma2_e
                sigma2_samples[idx] = sigma2
                pi_samples[idx] = pi
                delta_samples[idx] = delta.copy()
                idx += 1

        if sigma2_e <= sigma2_e_min * 1.01:
            diverged = True

    beta_mean = beta_samples.mean(axis=0)
    if not np.all(np.isfinite(beta_mean)):
        beta_mean = np.where(np.isfinite(beta_mean), beta_mean, 0.0)
    beta_orig, intercept = _back_transform_beta(beta_mean, X_sd, 1.0, y_mean, X_mean)

    return {
        "beta_mean": beta_orig,
        "intercept": intercept,
        "sigma2_e_mean": float(sigma2_e_samples.mean()),
        "sigma2_mean": float(sigma2_samples.mean()),
        "pi_mean": float(pi_samples.mean()),
        "delta_mean": delta_samples.mean(axis=0),
        "samples_beta": beta_samples,
        "n_saved": idx,
        "diverged": diverged,
    }

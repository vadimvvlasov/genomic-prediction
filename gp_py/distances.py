from __future__ import annotations

import numpy as np
import pandas as pd


def genomic_relationship_matrix(G: pd.DataFrame) -> pd.DataFrame:
    """Simple GRM: (X @ X.T) / p, where X is centered genotype matrix."""
    X = G.to_numpy(dtype=float)
    X = X - np.nanmean(X, axis=0, keepdims=True)
    denom = X.shape[1] if X.shape[1] > 0 else 1
    K = (X @ X.T) / denom
    return pd.DataFrame(K, index=G.index, columns=G.index)


def genomic_relationship_matrix_van_raden(
    G: pd.DataFrame, ploidy: int = 2
) -> pd.DataFrame:
    """Ploidy-aware GRM following VanRaden (2008) / Bell et al. (2017).

    Steps:
        1. G_star = ploidy * (G - 0.5)
        2. q = column means of G (allele frequencies)
        3. Z = G_star - ploidy * (q - 0.5)   # centered
        4. GRM = Z @ Z.T / (ploidy * sum(q * (1 - q)))
    """
    G_np = G.to_numpy(dtype=float)
    q = np.nanmean(G_np, axis=0)  # allele frequencies
    G_star = ploidy * (G_np - 0.5)
    Z = G_star - ploidy * (q - 0.5)
    denom = ploidy * np.sum(q * (1 - q))
    if denom <= 0:
        raise ValueError(
            "VanRaden GRM denominator is zero; "
            "check for monomorphic markers or invalid ploidy."
        )
    K = (Z @ Z.T) / denom
    return pd.DataFrame(K, index=G.index, columns=G.index)


def add_diagonal_load(K: pd.DataFrame, eps: float = 0.001) -> pd.DataFrame:
    """Add a small diagonal load to a symmetric matrix to ensure positive-definiteness.

    This mimics the R pattern: if solve() fails, add eps * diag(n) and retry.
    """
    arr = K.to_numpy(copy=True)  # ensure writable copy
    np.fill_diagonal(arr, np.diag(arr) + eps)
    return pd.DataFrame(arr, index=K.index, columns=K.columns)

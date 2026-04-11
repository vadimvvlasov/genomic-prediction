from __future__ import annotations

import numpy as np
import pandas as pd


def genomic_relationship_matrix(G: pd.DataFrame) -> pd.DataFrame:
    X = G.to_numpy(dtype=float)
    X = X - np.nanmean(X, axis=0, keepdims=True)
    denom = X.shape[1] if X.shape[1] > 0 else 1
    K = (X @ X.T) / denom
    return pd.DataFrame(K, index=G.index, columns=G.index)

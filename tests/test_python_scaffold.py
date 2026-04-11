import numpy as np
import pandas as pd
import pytest

from gp_py.models import fn_Bayes_A
from gp_py.schema import GPArgs, MergedData


def test_gp_args_defaults() -> None:
    args = GPArgs(fname_geno="geno.rds", fname_pheno="pheno.tsv", population="pop_1")
    assert args.n_folds == 10
    assert args.n_reps == 10
    assert args.bayes_backend == "auto"
    assert "ridge" in args.vec_models_to_test


def test_bayes_falls_back_without_r() -> None:
    rng = np.random.default_rng(123)
    G = pd.DataFrame(rng.normal(size=(20, 5)), index=[f"id_{i:03d}" for i in range(20)])
    beta = np.array([0.4, -0.2, 0.1, 0.0, 0.3])
    y = pd.Series(G.to_numpy() @ beta + rng.normal(scale=0.05, size=20), index=G.index)
    pop = pd.Series(["pop_1"] * 20, index=G.index)
    merged = MergedData(G=G, y=y, pop=pop, trait_name="trait")

    out = fn_Bayes_A(
        merged,
        vec_idx_training=list(range(15)),
        vec_idx_validation=list(range(15, 20)),
        other_params={"bayes_backend": "native"},
    )

    assert out["df_y_validation"].shape[0] == 5
    assert "corr" in out["list_perf"]


def test_bayes_rbridge_raises_without_r() -> None:
    rng = np.random.default_rng(123)
    G = pd.DataFrame(rng.normal(size=(12, 4)), index=[f"id_{i:03d}" for i in range(12)])
    y = pd.Series(rng.normal(size=12), index=G.index)
    pop = pd.Series(["pop_1"] * 12, index=G.index)
    merged = MergedData(G=G, y=y, pop=pop, trait_name="trait")

    with pytest.raises(Exception):
        fn_Bayes_A(
            merged,
            vec_idx_training=list(range(8)),
            vec_idx_validation=list(range(8, 12)),
            other_params={"bayes_backend": "rbridge"},
        )

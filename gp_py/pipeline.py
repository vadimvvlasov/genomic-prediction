from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from gp_py.cv import MODEL_REGISTRY, fn_cross_validation_within_population
from gp_py.io import (
    fn_filter_genotype,
    fn_filter_phenotype,
    fn_load_genotype,
    fn_load_phenotype,
    fn_merge_genotype_and_phenotype,
)
from gp_py.schema import GPArgs


def gp(args: GPArgs) -> str:
    G = fn_load_genotype(args.fname_geno, verbose=args.verbose)
    list_pheno = fn_load_phenotype(
        args.fname_pheno,
        sep=args.pheno_sep,
        header=args.pheno_header,
        idx_col_id=args.pheno_idx_col_id,
        idx_col_pop=args.pheno_idx_col_pop,
        idx_col_y=args.pheno_idx_col_y,
        na_strings=args.pheno_vec_na_strings,
        verbose=args.verbose,
    )

    Gf = fn_filter_genotype(G, verbose=args.verbose)
    phf = fn_filter_phenotype(
        list_pheno,
        remove_outliers=args.pheno_bool_remove_outliers,
        remove_na=False,
        verbose=args.verbose,
    )
    merged = fn_merge_genotype_and_phenotype(Gf, phf, COVAR=None, verbose=args.verbose)

    known_mask = merged.y.notna()
    if known_mask.sum() < 2:
        raise ValueError("Phenotype data is missing or too sparse for model training")

    merged_known = type(merged)(
        G=merged.G.loc[known_mask].copy(),
        y=merged.y.loc[known_mask].copy(),
        pop=merged.pop.loc[known_mask].copy(),
        trait_name=merged.trait_name,
        covar=merged.covar,
    )

    cv_out = fn_cross_validation_within_population(
        merged_known,
        n_folds=args.n_folds,
        n_reps=args.n_reps,
        vec_models_to_test=args.vec_models_to_test,
        bool_parallel=args.bool_parallel,
        bayes_backend=args.bayes_backend,
        max_mem_gb=args.max_mem_gb,
        n_threads=args.n_threads,
        dir_output=args.dir_output,
        verbose=args.verbose,
    )

    metrics = cv_out["METRICS_WITHIN_POP"]
    if metrics.empty:
        best_model = None
    else:
        best_model = (
            metrics.groupby("model", as_index=False)["corr"]
            .mean()
            .sort_values("corr", ascending=False)
            .iloc[0]["model"]
        )

    genomic_predictions = pd.DataFrame()
    if best_model is not None:
        known_idx = np.where(known_mask.to_numpy())[0].tolist()
        missing_idx = np.where(~known_mask.to_numpy())[0].tolist()
        if len(missing_idx) > 0:
            fn_model = MODEL_REGISTRY.get(str(best_model))
            if fn_model is not None:
                pred_out = fn_model(
                    merged,
                    known_idx,
                    missing_idx,
                    other_params={"n_folds": 10, "bayes_backend": args.bayes_backend},
                    verbose=args.verbose,
                )
                genomic_predictions = pred_out["df_y_validation"].copy()
                genomic_predictions["model"] = best_model

    output = {
        "TRAIT_NAME": merged.trait_name,
        "POPULATION": args.population,
        "METRICS_WITHIN_POP": cv_out["METRICS_WITHIN_POP"],
        "YPRED_WITHIN_POP": cv_out["YPRED_WITHIN_POP"],
        "METRICS_ACROSS_POP_BULK": pd.DataFrame(),
        "YPRED_ACROSS_POP_BULK": pd.DataFrame(),
        "METRICS_ACROSS_POP_PAIRWISE": pd.DataFrame(),
        "YPRED_ACROSS_POP_PAIRWISE": pd.DataFrame(),
        "METRICS_ACROSS_POP_LOPO": pd.DataFrame(),
        "YPRED_ACROSS_POP_LOPO": pd.DataFrame(),
        "GENOMIC_PREDICTIONS": genomic_predictions,
        "ADDITIVE_GENETIC_EFFECTS": {},
    }

    out_dir = Path(args.dir_output or ".")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"GENOMIC_PREDICTIONS_OUTPUT-{merged.trait_name}-{args.population}.parquet"
    pd.DataFrame({"key": list(output.keys())}).to_parquet(out_path)
    return str(out_path)

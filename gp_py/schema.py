from dataclasses import dataclass
from typing import Optional

import pandas as pd


@dataclass
class GPArgs:
    fname_geno: str
    fname_pheno: str
    population: str
    dir_output: Optional[str] = None
    pheno_sep: str = "\t"
    pheno_header: bool = True
    pheno_idx_col_id: int = 1
    pheno_idx_col_pop: int = 2
    pheno_idx_col_y: int = 3
    pheno_vec_na_strings: tuple[str, ...] = ("", "-", "NA", "na", "NaN", "missing", "MISSING")
    pheno_bool_remove_outliers: bool = False
    pheno_bool_remove_na: bool = False
    bool_within: bool = True
    bool_across: bool = False
    n_folds: int = 10
    n_reps: int = 10
    vec_models_to_test: tuple[str, ...] = (
        "ridge",
        "lasso",
        "elastic_net",
        "Bayes_A",
        "Bayes_B",
        "Bayes_C",
        "gBLUP",
    )
    bool_parallel: bool = True
    max_mem_gb: float = 15.0
    n_threads: int = 2
    verbose: bool = True


@dataclass
class MergedData:
    G: pd.DataFrame
    y: pd.Series
    pop: pd.Series
    trait_name: str
    covar: Optional[pd.DataFrame] = None

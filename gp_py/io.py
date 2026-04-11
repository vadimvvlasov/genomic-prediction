from __future__ import annotations

from pathlib import Path

import pandas as pd

from gp_py.schema import MergedData


def fn_load_genotype(fname_geno: str, *, verbose: bool = False) -> pd.DataFrame:
    path = Path(fname_geno)
    if not path.exists():
        raise FileNotFoundError(f"Genotype file does not exist: {fname_geno}")

    suffix = path.suffix.lower()
    if suffix == ".rds":
        try:
            import pyreadr  # type: ignore
        except ImportError as exc:
            raise ImportError("Install optional dependency group 'io' for RDS support") from exc

        obj = pyreadr.read_r(str(path))
        if len(obj.keys()) == 0:
            raise ValueError(f"No objects found in RDS file: {fname_geno}")
        G = next(iter(obj.values()))
        if not isinstance(G, pd.DataFrame):
            raise TypeError("Expected genotype RDS to contain a data.frame/matrix-like object")
        return G

    if suffix in {".tsv", ".txt", ".csv"}:
        sep = "\t" if suffix in {".tsv", ".txt"} else ","
        return pd.read_csv(path, sep=sep)

    raise ValueError(f"Unsupported genotype format: {suffix}")


def fn_load_phenotype(
    fname_pheno: str,
    *,
    sep: str = "\t",
    header: bool = True,
    idx_col_id: int = 1,
    idx_col_pop: int = 2,
    idx_col_y: int = 3,
    na_strings: tuple[str, ...] = ("", "-", "NA", "na", "NaN", "missing", "MISSING"),
    verbose: bool = False,
) -> dict:
    path = Path(fname_pheno)
    if not path.exists():
        raise FileNotFoundError(f"Phenotype file does not exist: {fname_pheno}")

    df = pd.read_csv(path, sep=sep, header=0 if header else None, na_values=list(na_strings))
    id_col = df.columns[idx_col_id - 1]
    pop_col = df.columns[idx_col_pop - 1]
    y_col = df.columns[idx_col_y - 1]

    return {
        "df": df,
        "id_col": id_col,
        "pop_col": pop_col,
        "y_col": y_col,
        "trait_name": str(y_col),
    }


def fn_filter_genotype(
    G: pd.DataFrame, *, maf: float = 0.01, sdev_min: float = 1e-4, verbose: bool = False
) -> pd.DataFrame:
    numeric = G.select_dtypes(include=["number"])
    if numeric.shape[1] == 0:
        raise ValueError("Genotype matrix has no numeric marker columns")

    maf_mask = numeric.apply(lambda c: min(c.mean(skipna=True), 1 - c.mean(skipna=True)) >= maf)
    sdev_mask = numeric.std(axis=0, skipna=True) >= sdev_min
    keep = maf_mask & sdev_mask
    Gf = numeric.loc[:, keep]
    if Gf.shape[1] == 0:
        raise ValueError("All genotype loci were filtered out")
    return Gf


def fn_filter_phenotype(
    list_pheno: dict,
    *,
    remove_outliers: bool = False,
    remove_na: bool = False,
    verbose: bool = False,
) -> dict:
    df = list_pheno["df"].copy()
    y_col = list_pheno["y_col"]
    if remove_na:
        df = df[df[y_col].notna()].copy()
    if remove_outliers:
        q1 = df[y_col].quantile(0.25)
        q3 = df[y_col].quantile(0.75)
        iqr = q3 - q1
        low = q1 - 1.5 * iqr
        high = q3 + 1.5 * iqr
        df = df[df[y_col].between(low, high) | df[y_col].isna()].copy()
    list_pheno["df"] = df
    return list_pheno


def fn_merge_genotype_and_phenotype(
    G: pd.DataFrame, list_pheno: dict, COVAR=None, *, verbose: bool = False
) -> MergedData:
    df = list_pheno["df"]
    id_col = list_pheno["id_col"]
    pop_col = list_pheno["pop_col"]
    y_col = list_pheno["y_col"]

    if id_col not in df.columns:
        raise ValueError(f"Phenotype ID column is missing: {id_col}")

    if "id" in G.columns:
        G_idx = G.set_index("id")
    else:
        G_idx = G.copy()
        if G_idx.index.name is None:
            G_idx.index = G_idx.index.astype(str)

    ph = df[[id_col, pop_col, y_col]].copy()
    ph[id_col] = ph[id_col].astype(str)
    merged = ph.merge(G_idx, left_on=id_col, right_index=True, how="inner")
    if merged.empty:
        raise ValueError("No overlapping sample IDs between phenotype and genotype")

    y = merged[y_col]
    pop = merged[pop_col]
    Gm = merged.drop(columns=[id_col, pop_col, y_col])

    return MergedData(
        G=Gm,
        y=y,
        pop=pop,
        trait_name=list_pheno["trait_name"],
        covar=COVAR,
    )


def fn_estimate_memory_footprint(
    X: MergedData,
    *,
    n_models: int,
    n_folds: int,
    n_reps: int,
    memory_requested_gb: float,
    verbose: bool = False,
) -> dict:
    bytes_per_value = 8
    x_size = X.G.shape[0] * X.G.shape[1] * bytes_per_value
    multiplier = max(1, n_models * n_folds * n_reps)
    total = x_size * multiplier
    requested = int(memory_requested_gb * (1024**3))
    n_threads = max(1, requested // max(x_size, 1))
    return {
        "size_x_bytes": x_size,
        "size_total_bytes": total,
        "n_threads": n_threads,
    }

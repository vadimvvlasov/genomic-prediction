"""Visualization utilities for genomic prediction results."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

try:
    import matplotlib
    matplotlib.use("Agg")  # non-interactive backend
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


def _compute_regression(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    """Compute OLS regression line and correlation."""
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    yt = y_true[mask]
    yp = y_pred[mask]
    if yt.size < 2:
        return {"slope": np.nan, "intercept": np.nan, "r2": np.nan, "r": np.nan, "n": 0}

    slope, intercept = np.polyfit(yt, yp, 1)
    r = float(np.corrcoef(yt, yp)[0, 1])
    ss_res = float(np.sum((yp - (slope * yt + intercept)) ** 2))
    ss_tot = float(np.sum((yp - np.mean(yp)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    return {"slope": slope, "intercept": intercept, "r2": r2, "r": r, "n": yt.size}


def plot_predicted_vs_actual(
    df: pd.DataFrame,
    *,
    y_true_col: str = "y_true",
    y_pred_col: str = "y_pred",
    color_by: str | None = None,
    title: str = "Predicted vs Actual Breeding Values",
    add_regression: bool = True,
    add_identity: bool = True,
    figsize: tuple[float, float] = (7, 6),
    dpi: int = 150,
    point_alpha: float = 0.6,
    point_size: float = 12,
) -> "matplotlib.figure.Figure":
    """Scatter plot of predicted vs actual breeding values with regression line.

    Mirrors the Shiny app observed-vs-predicted scatter
    (inst/plot_gs_gp/app.R: fn_within_scatterplot_server).

    Parameters
    ----------
    df :
        DataFrame with at least y_true and y_pred columns.
        Optionally: ``model``, ``pop``, ``rep``, ``fold`` for annotations.
    y_true_col, y_pred_col :
        Column names for observed and predicted values.
    color_by :
        Column name to color points by (e.g. ``"model"``, ``"pop"``).
    title :
        Plot title.
    add_regression :
        Draw OLS regression line with annotation.
    add_identity :
        Draw y=x reference line (dashed).
    figsize, dpi, point_alpha, point_size :
        Plot aesthetics.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if not HAS_MPL:
        raise ImportError(
            "Install matplotlib: uv add matplotlib"
        )

    yt = df[y_true_col].to_numpy(dtype=float)
    yp = df[y_pred_col].to_numpy(dtype=float)
    valid = np.isfinite(yt) & np.isfinite(yp)
    yt, yp = yt[valid], yp[valid]

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # ── Scatter ─────────────────────────────────────────────────────
    if color_by and color_by in df.columns:
        labels = df.loc[valid, color_by].astype(str)
        unique_labels = sorted(labels.unique())
        cmap = plt.get_cmap("tab10" if len(unique_labels) <= 10 else "tab20")
        for i, label in enumerate(unique_labels):
            mask = labels == label
            ax.scatter(
                yt[mask], yp[mask],
                c=[cmap(i / max(len(unique_labels) - 1, 1))],
                s=point_size, alpha=point_alpha,
                label=str(label),
                edgecolors="none",
                zorder=2,
            )
        if len(unique_labels) <= 15:
            ax.legend(title=color_by, fontsize=6, title_fontsize=7,
                      loc="upper left", framealpha=0.8)
    else:
        ax.scatter(yt, yp, c="steelblue", s=point_size, alpha=point_alpha,
                   edgecolors="none", zorder=2)

    # ── Identity line ───────────────────────────────────────────────
    if add_identity:
        lims = [min(yt.min(), yp.min()), max(yt.max(), yp.max())]
        ax.plot(lims, lims, "k--", linewidth=1, alpha=0.4, zorder=1,
                label="y = x")

    # ── Regression line + annotation ────────────────────────────────
    if add_regression:
        stats = _compute_regression(yt, yp)
        x_range = np.linspace(yt.min(), yt.max(), 200)
        y_line = stats["slope"] * x_range + stats["intercept"]
        ax.plot(x_range, y_line, "r-", linewidth=1.5, alpha=0.8, zorder=3,
                label=f"OLS fit")

        # Annotation box
        n_unique_models = df["model"].nunique() if "model" in df.columns else 1
        n_unique_pops = df["pop"].nunique() if "pop" in df.columns else 1
        n_reps = df["rep"].nunique() if "rep" in df.columns else 1
        n_folds = df["fold"].nunique() if "fold" in df.columns else 1

        annot_lines = [
            f"n = {stats['n']}",
            f"r = {stats['r']:.4f}",
            f"R\u00B2 = {stats['r2']:.4f}",
            f"slope = {stats['slope']:.4f}",
            f"intercept = {stats['intercept']:.4f}",
        ]
        if n_unique_models > 1:
            annot_lines.append(f"models = {n_unique_models}")
        if n_unique_pops > 1:
            annot_lines.append(f"pops = {n_unique_pops}")
        if n_reps > 1:
            annot_lines.append(f"reps = {n_reps}, folds = {n_folds}")

        textstr = "\n".join(annot_lines)
        props = dict(boxstyle="round,pad=0.4", facecolor="white",
                     edgecolor="gray", alpha=0.9)
        ax.text(0.02, 0.98, textstr, transform=ax.transAxes,
                fontsize=8, verticalalignment="top", fontfamily="monospace",
                bbox=props, zorder=4)

    # ── Labels & title ──────────────────────────────────────────────
    ax.set_xlabel("Observed (Actual Breeding Value)", fontsize=10)
    ax.set_ylabel("Predicted (GEBV)", fontsize=10)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(True, alpha=0.3, zorder=0)
    ax.set_axisbelow(True)

    fig.tight_layout()
    return fig


def plot_predicted_vs_actual_by_model(
    df: pd.DataFrame,
    *,
    y_true_col: str = "y_true",
    y_pred_col: str = "y_pred",
    model_col: str = "model",
    ncols: int = 3,
    figsize_per_panel: tuple[float, float] = (4.5, 4),
    **kwargs,
) -> "matplotlib.figure.Figure":
    """Faceted scatter plot: one panel per model.

    Parameters
    ----------
    df :
        Must contain ``model_col``, ``y_true_col``, ``y_pred_col``.
    ncols :
        Number of columns in the facet grid.
    figsize_per_panel :
        (width, height) per subplot in inches.
    **kwargs :
        Passed to :func:`plot_predicted_vs_actual`.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if not HAS_MPL:
        raise ImportError("Install matplotlib: uv add matplotlib")

    models = sorted(df[model_col].dropna().unique())
    n_models = len(models)
    if n_models == 0:
        raise ValueError(f"No models found in column '{model_col}'")

    nrows = int(np.ceil(n_models / ncols))
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows),
        squeeze=False,
    )

    for idx, model in enumerate(models):
        row, col = divmod(idx, ncols)
        ax = axes[row, col]

        sub = df[df[model_col] == model]
        yt = sub[y_true_col].to_numpy(dtype=float)
        yp = sub[y_pred_col].to_numpy(dtype=float)
        valid = np.isfinite(yt) & np.isfinite(yp)
        yt, yp = yt[valid], yp[valid]

        ax.scatter(yt, yp, c="steelblue", s=kwargs.get("point_size", 10),
                   alpha=kwargs.get("point_alpha", 0.5),
                   edgecolors="none", zorder=2)

        if kwargs.get("add_identity", True):
            lims = [min(yt.min(), yp.min()), max(yt.max(), yp.max())]
            ax.plot(lims, lims, "k--", linewidth=0.8, alpha=0.3, zorder=1)

        if kwargs.get("add_regression", True):
            stats = _compute_regression(yt, yp)
            x_range = np.linspace(yt.min(), yt.max(), 200)
            ax.plot(x_range, stats["slope"] * x_range + stats["intercept"],
                    "r-", linewidth=1.2, alpha=0.7, zorder=3)

            annot_text = f"r = {stats['r']:.3f}\nR\u00B2 = {stats['r2']:.3f}\nn = {stats['n']}"
            props = dict(boxstyle="round,pad=0.3", facecolor="white",
                         edgecolor="gray", alpha=0.85)
            ax.text(0.97, 0.03, annot_text, transform=ax.transAxes,
                    fontsize=7, verticalalignment="bottom",
                    horizontalalignment="right", fontfamily="monospace",
                    bbox=props, zorder=4)

        ax.set_xlabel("Observed", fontsize=8)
        ax.set_ylabel("Predicted", fontsize=8)
        ax.set_title(model, fontsize=9, fontweight="bold")
        ax.grid(True, alpha=0.25, zorder=0)
        ax.set_axisbelow(True)

    # Hide unused subplots
    for idx in range(n_models, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row, col].set_visible(False)

    fig.suptitle("Predicted vs Actual Breeding Values by Model",
                 fontsize=12, fontweight="bold", y=1.01)
    fig.tight_layout()
    return fig


def save_plot(
    fig: "matplotlib.figure.Figure",
    path: str | Path,
    *,
    dpi: int = 150,
    bbox_inches: str = "tight",
) -> str:
    """Save a matplotlib figure to file.

    Supports PNG, SVG, PDF based on file extension.
    """
    if not HAS_MPL:
        raise ImportError("Install matplotlib: uv add matplotlib")

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(path), dpi=dpi, bbox_inches=bbox_inches)
    return str(path)

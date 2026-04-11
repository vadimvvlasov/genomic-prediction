from __future__ import annotations

import typer

from gp_py.pipeline import gp
from gp_py.schema import GPArgs

app = typer.Typer(help="gp Python pipeline scaffold")

AVAILABLE_MODELS = (
    "ridge",
    "lasso",
    "elastic_net",
    "Bayes_A",
    "Bayes_B",
    "Bayes_C",
    "gBLUP",
)


@app.callback(invoke_without_command=True)
def main(ctx: typer.Context) -> None:
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
        raise typer.Exit(code=0)


@app.command()
def run(
    fname_geno: str = typer.Option(..., help="Genotype file path"),
    fname_pheno: str = typer.Option(..., help="Phenotype file path"),
    population: str = typer.Option(..., help="Population for within-pop workflow"),
    dir_output: str = typer.Option("./output", help="Output directory"),
    pheno_idx_col_y: int = typer.Option(3, help="1-based trait column index in phenotype file"),
    models: list[str] = typer.Option(
        list(AVAILABLE_MODELS),
        "--models",
        "-m",
        help=f"Models to test. Repeatable. [default: all {len(AVAILABLE_MODELS)} models]",
    ),
    n_folds: int = typer.Option(10, help="K-fold count"),
    n_reps: int = typer.Option(10, help="Replication count"),
    n_threads: int = typer.Option(2, help="Number of parallel threads for CV"),
    parallel: bool = typer.Option(True, "--parallel/--no-parallel", help="Enable parallel CV"),
    bayes_backend: str = typer.Option("auto", help="Bayes backend: auto|native|rbridge"),
    pheno_sep: str = typer.Option("\t", help="Phenotype file separator"),
    verbose: bool = typer.Option(True, help="Verbose logs"),
) -> None:
    args = GPArgs(
        fname_geno=fname_geno,
        fname_pheno=fname_pheno,
        population=population,
        dir_output=dir_output,
        pheno_idx_col_y=pheno_idx_col_y,
        pheno_sep=pheno_sep,
        n_folds=n_folds,
        n_reps=n_reps,
        n_threads=n_threads,
        bool_parallel=parallel,
        vec_models_to_test=tuple(models),
        bayes_backend=bayes_backend,
        verbose=verbose,
    )
    out = gp(args)
    typer.echo(out)


if __name__ == "__main__":
    app()

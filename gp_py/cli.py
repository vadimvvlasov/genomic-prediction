from __future__ import annotations

import typer

from gp_py.pipeline import gp
from gp_py.schema import GPArgs

app = typer.Typer(help="gp Python pipeline scaffold")


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
    n_folds: int = typer.Option(2, help="K-fold count"),
    n_reps: int = typer.Option(2, help="Replication count"),
    verbose: bool = typer.Option(True, help="Verbose logs"),
) -> None:
    args = GPArgs(
        fname_geno=fname_geno,
        fname_pheno=fname_pheno,
        population=population,
        dir_output=dir_output,
        pheno_idx_col_y=pheno_idx_col_y,
        n_folds=n_folds,
        n_reps=n_reps,
        verbose=verbose,
    )
    out = gp(args)
    typer.echo(out)


if __name__ == "__main__":
    app()

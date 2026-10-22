import typer


from resdel.cli.prepare import app as prepare_app
from resdel.cli.transform import app as transform_app
from resdel.cli.analyze import app as analyze_app

app = typer.Typer()

app.add_typer(
    prepare_app,
    name="prepare",
    help="Build wild type and mutant GROMACS topologies."
)

app.add_typer(
    transform_app,
    name="transform",
    help="Build the alchemical transformation topology."
)

app.add_typer(
    analyze_app,
    name="analyze",
    help="Estimate the free energy difference with BAR."
)

#app.add_typer(
#    extract_top_app,
#    name="extract_topology"
#)

if __name__ == "__main__":
    app()

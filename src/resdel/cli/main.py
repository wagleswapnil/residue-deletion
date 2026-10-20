import typer

from resdel.config.loader import load_config

from resdel.cli.prepare import app as prepare_app
from resdel.cli.transform import app as transform_app
from resdel.cli.extract_top import app as extract_top_app

app = typer.Typer()

app.add_typer(
    prepare_app,
    name="prepare"
)

app.add_typer(
    transform_app,
    name="transform"
)

app.add_typer(
    extract_top_app,
    name="extract_topology"
)

if __name__ == "__main__":
    app()

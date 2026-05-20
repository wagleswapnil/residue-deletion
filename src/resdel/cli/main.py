import typer

from resdel.config.loader import load_config

from resdel.cli.prepare import app as prepare_app
from resdel.cli.transform import app as transform_app

app = typer.Typer()

app.add_typer(
    prepare_app,
    name="prepare"
)

app.add_typer(
    transform_app,
    name="transform"
)

if __name__ == "__main__":
    app()

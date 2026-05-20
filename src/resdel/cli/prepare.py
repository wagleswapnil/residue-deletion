import typer
from resdel.preparation.workflow import run_prepare_workflow
from resdel.config.loader import load_config


app = typer.Typer()

@app.callback(invoke_without_command=True)
def prepare(config_file: str):

    config = load_config(config_file)

    typer.echo(
        f"Preparing system {config.system.name}"
    )

    run_prepare_workflow(config)
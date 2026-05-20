import typer

from resdel.config.loader import load_config

app = typer.Typer()

@app.callback(invoke_without_command=True)
def tranform(config_file: str):

    config = load_config(config_file)

    typer.echo(
        f"Transforming topology {config.system.name}"
    )

    

import typer
from resdel.preparation.workflow import run_prepare_workflow
from resdel.config.loader import load_config
from resdel.workflow.paths import OutputPaths

app = typer.Typer()

@app.callback(invoke_without_command=True)
def prepare(config_file: str):

    config = load_config(config_file)

    typer.echo(
        f"Preparing system {config.system.name}"
    )

    output_dir = (config.io.output_dir
                  or (config.system.name + "_output" if config.system.name else None)
                  or "./residue_deletion_output")

    paths = OutputPaths(output_dir)

    run_prepare_workflow(config, paths)
import os
import typer
from resdel.config.loader import load_config
from resdel.analysis.workflow import run_analyze_workflow 

app = typer.Typer()

@app.callback(invoke_without_command=True)
def prepare(config_file: str):

    config = load_config(config_file)

    typer.echo(
        f"Analyzing system {config.system.name}"
    )

    output_dir = (config.io.output_dir
                  or (config.system.name + "_output" if config.system.name else None)
                  or "./residue_deletion_output")

    run_analyze_workflow(os.path.join(output_dir, "transformation"))

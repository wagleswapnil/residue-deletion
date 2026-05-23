import typer
from resdel.transformations.workflow import run_transform_workflow
from resdel.config.loader import load_config
from resdel.workflow.paths import OutputPaths


app = typer.Typer()

@app.callback(invoke_without_command=True)
def transform(config_file: str):

    config = load_config(config_file)

    typer.echo(
        f"Transforming topology {config.system.name}"
    )

    output_dir = (config.io.output_dir
                  or (config.system.name + "_output" if config.system.name else None)
                  or "./residue_deletion_output")
    
    paths = OutputPaths(output_dir)

    run_transform_workflow(config, paths)
    

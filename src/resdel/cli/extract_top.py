import typer
from resdel.config.loader import load_config
from resdel.workflow.paths import OutputPaths
from resdel.extract_topology.workflow import run_extract_topology_workflow

app = typer.Typer()

@app.callback(invoke_without_command=True)
def extract_topology(config_file: str):
    """
    Extracts the topology of a single residue from a molecular system based on the provided configuration file.
    
    Args:
        config_file (str): Path to the configuration file containing system and residue information.
    """
    # Load the configuration
    config = load_config(config_file)

    typer.echo(
        f"Extracting topology for residue {config.system.residue_to_delete} from system {config.system.name}"
    )   

    output_dir = (config.io.output_dir
                  or (config.system.name + "_output" if config.system.name else None)
                  or "./residue_deletion_output")
    
    paths = OutputPaths(output_dir)

    run_extract_topology_workflow(config, paths)
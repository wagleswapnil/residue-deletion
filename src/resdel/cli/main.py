import typer

from resdel.topology.parser import Parser
from resdel.transformations.topology_transformer import TopologyTransformer

app = typer.Typer()

@app.command()
def transform(
    topA: str,
    topB: str,
    residue: int,
    output: str = "output.top",
    edge1_steps: int = 15
):

    parserA = Parser(file_path=topA)
    parserA.parse_topology()

    parserB = Parser(file_path=topB)
    parserB.parse_topology()

    transformer = TopologyTransformer(
        topA=parserA.top,
        topB=parserB.top,
        residue_to_delete=str(residue),
        edge1_steps=edge1_steps
    )

    transformer.generate_resdel_topology()

    typer.echo("Topology transformation complete.")


if __name__ == "__main__":
    app()

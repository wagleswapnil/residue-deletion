from resdel.transformations.topology_transformer import TopologyTransformer
from resdel.topology.parser import Parser
from resdel.topology.writer import Writer

def run_transform_workflow(config, paths):
    topA = Parser(file_path=paths.topology_file("wt"))
    topA.parse_topology()

    topB = Parser(file_path=paths.topology_file("mutant"))
    topB.parse_topology()

    topology_transformation = TopologyTransformer(topA.top, topB.top, config, paths)
    topology_transformation.generate_resdel_topology()
    topology_transformation.write_topology_output(str(paths.topology_file("transformation")))
    return
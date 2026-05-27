from resdel.transformations.topology_transformer import TopologyTransformer
from resdel.topology.parser import Parser
from resdel.topology.writer import Writer
from resdel.utils.utils import *

def run_transform_workflow(config, paths):
    topA = Parser(file_path=paths.topology_file("wt/solvated"))
    topA.parse_topology()

    topB = Parser(file_path=paths.topology_file("mutant/solvated"))
    topB.parse_topology()

    topology_transformation = TopologyTransformer(topA.top, topB.top, config, paths)
    topology_transformation.generate_resdel_topology()
    topology_transformation.write_topology_output(str(paths.topology_file("transformation/solvated")))

    copy_structure(str(paths.structure_PDBfile("wt/solvated")), str(paths.structure_PDBfile("transformation/solvated")))

    generate_vaccuum_topology_from_solvent_topology(str(paths.topology_file("transformation/solvated")), str(paths.topology_file("transformation/vaccuum")))
    generate_vaccuum_structure_from_solvent_structure(str(paths.structure_PDBfile("transformation/solvated")), str(paths.structure_PDBfile("transformation/vaccuum")))
    return
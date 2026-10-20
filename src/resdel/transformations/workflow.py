from resdel.transformations.topology_transformer import TopologyTransformer
from resdel.topology.parser import Parser
from resdel.utils.utils import copy_structure, generate_vacuum_structure_from_solvent_structure, generate_vacuum_topology_from_solvent_topology

# The function takes in "config" and "paths" objects, infers the paths of the wt and mutant topologies from the paths object, and runs a 
# topology transformation on these.  
def run_transform_workflow(config, paths):
    top_wt = Parser(file_path=paths.topology_file("wt/solvated"))
    top_wt.parse_topology()

    top_mutant = Parser(file_path=paths.topology_file("mutant/solvated"))
    top_mutant.parse_topology()

    topology_transformation = TopologyTransformer(top_wt.top, top_mutant.top, config, paths)
    topology_transformation.generate_resdel_topology()
    #topology_transformation.add_posre_section_to_topology()
    topology_transformation.write_topology_output(str(paths.topology_file("transformation/solvated")))

    copy_structure(str(paths.structure_PDBfile("wt/solvated")), str(paths.structure_PDBfile("transformation/solvated")))

    generate_vacuum_topology_from_solvent_topology(str(paths.topology_file("transformation/solvated")), str(paths.topology_file("transformation/vacuum")))
    generate_vacuum_structure_from_solvent_structure(str(paths.structure_PDBfile("transformation/solvated")), str(paths.structure_PDBfile("transformation/vacuum")))
    return
from resdel.extract_topology.extract_single_res_topology import TopologyExtractor
from resdel.utils.utils import generate_vacuum_topology_from_solvent_topology

def run_extract_topology_workflow(config, paths):
    """
    Extracts the topology of a single residue from a molecular system based on the provided configuration and paths.
    
    Args:
        config: Configuration object containing system and residue information.
        paths: OutputPaths object specifying where to save the extracted topology.
    """
    residue_index = config.system.residue_to_delete
    molecule_name = config.transform.states["A"].molecule_name
    output_path = str(paths.topology_file("extracted_residue/solvated"))

    extractor = TopologyExtractor(str(paths.topology_file("wt/solvated")))
    extractor.extract_molecule_topology(molecule_name, residue_index)
    extractor.build_residue_topology()
    extractor.write_residue_topology(output_path)

    generate_vacuum_topology_from_solvent_topology(str(paths.topology_file("extracted_residue/solvated")), str(paths.topology_file("extracted_residue/vacuum")))
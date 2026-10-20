from resdel.extract_topology.extract_top_section import ExtractTopSection
from resdel.utils.utils import generate_vacuum_topology_from_solvent_topology
import sys

in_path = sys.argv[1]
begin_idx = sys.argv[2]
end_idx = sys.argv[3]
first_residue_idx = sys.argv[4]
output_path = sys.argv[5]
output_path_2 = sys.argv[6]

"""
Extracts a topology section, reindexing atoms and residues.

Args:
    config: Configuration object containing system and residue information.
    paths: OutputPaths object specifying where to save the extracted topology.
"""

molecule_name = None  # Default to None, will be set to "system1" if not provided

extractor = ExtractTopSection(
    topology_path=in_path,
    molecule_name=molecule_name,
    begin_idx=begin_idx,
    end_idx=end_idx,
    first_residue_idx=first_residue_idx
)
extractor.extract_topology_sections()
#extractor.write_residue_topology(output_path)

extractor.make_last_terminal_residue_dual_topology(last_residue=3)  # Create dual topology for the terminal residue
extractor.write_residue_topology(output_path_2)

#generate_vacuum_topology_from_solvent_topology(str(paths.topology_file("extracted_topology/solvated")), str(paths.topology_file("extracted_topology/vacuum")))
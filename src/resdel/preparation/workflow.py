from .peptide_builder import PeptideSystemBuilder
from resdel.utils.utils import generate_vacuum_structure_from_solvent_structure, generate_vacuum_topology_from_solvent_topology

# This function is executed when only a sequence of peptide/protein, along with the residue number to be deleted, 
# is given as an input. It generates the 3D structures and topoogies of the (wt) protein and the mutant, in solvent and in vacumm.  
def run_prepare_workflow(config, paths):
    sequence = config.system.sequence
    residue_to_delete = config.system.residue_to_delete

    system_builder = PeptideSystemBuilder(sequence, residue_to_delete)
    system_builder.generate_wt_peptide_structure(paths.structure_PDBfile("wt"))
    system_builder.generate_mutant_structure(paths.structure_PDBfile("mutant"))

    system_builder.generate_topology_from_structure(str(paths.structure_PDBfile("wt")), str(paths.structure_PDBfile("wt/solvated")), str(paths.topology_file("wt/solvated")))
    #system_builder.add_posre_section_to_topology(str(paths.topology_file("wt/solvated")))
    system_builder.generate_topology_from_structure(str(paths.structure_PDBfile("mutant")), str(paths.structure_PDBfile("mutant/solvated")), str(paths.topology_file("mutant/solvated")))
    #system_builder.add_posre_section_to_topology(str(paths.topology_file("mutant/solvated")))

    generate_vacuum_structure_from_solvent_structure(str(paths.structure_PDBfile("wt/solvated")), str(paths.structure_PDBfile("wt/vacuum")))
    generate_vacuum_structure_from_solvent_structure(str(paths.structure_PDBfile("mutant/solvated")), str(paths.structure_PDBfile("mutant/vacuum")))

    generate_vacuum_topology_from_solvent_topology(str(paths.topology_file("wt/solvated")), str(paths.topology_file("wt/vacuum")))
    generate_vacuum_topology_from_solvent_topology(str(paths.topology_file("mutant/solvated")), str(paths.topology_file("mutant/vacuum")))
    return
    



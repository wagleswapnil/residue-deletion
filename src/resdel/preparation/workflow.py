from .peptide_builder import PeptideSystemBuilder
from resdel.utils.utils import * 

def run_prepare_workflow(config, paths):
    sequence = config.system.sequence
    residue_to_delete = config.system.residue_to_delete

    system_builder = PeptideSystemBuilder(sequence, residue_to_delete)
    system_builder.generate_wt_peptide_structure(paths.structure_PDBfile("wt"))
    system_builder.generate_mutant_structure(paths.structure_PDBfile("mutant"))

    system_builder.generate_topology_from_structure(str(paths.structure_PDBfile("wt")), str(paths.structure_PDBfile("wt/solvated")), str(paths.topology_file("wt/solvated")))
    system_builder.generate_topology_from_structure(str(paths.structure_PDBfile("mutant")), str(paths.structure_PDBfile("mutant/solvated")), str(paths.topology_file("mutant/solvated")))

    generate_vaccuum_structure_from_solvent_structure(str(paths.structure_PDBfile("wt/solvated")), str(paths.structure_PDBfile("wt/vaccuum")))
    generate_vaccuum_structure_from_solvent_structure(str(paths.structure_PDBfile("mutant/solvated")), str(paths.structure_PDBfile("mutant/vaccuum")))


    generate_vaccuum_topology_from_solvent_topology(str(paths.topology_file("wt/solvated")), str(paths.topology_file("wt/vaccuum")))
    generate_vaccuum_topology_from_solvent_topology(str(paths.topology_file("mutant/solvated")), str(paths.topology_file("mutant/vaccuum")))
    return
    



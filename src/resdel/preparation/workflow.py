from .peptide_builder import PeptideSystemBuilder

def run_prepare_workflow(config, paths):
    sequence = config.system.sequence
    residue_to_delete = config.system.residue_to_delete

    system_builder = PeptideSystemBuilder(sequence, residue_to_delete)
    system_builder.generate_wt_peptide_structure(paths.structure_PDBfile("wt"))
    system_builder.generate_mutant_structure(paths.structure_PDBfile("mutant"))

    system_builder.generate_topology_from_structure(str(paths.structure_PDBfile("wt")), str(paths.topology_file("wt")))
    system_builder.generate_topology_from_structure(str(paths.structure_PDBfile("mutant")), str(paths.topology_file("mutant")))
    return
    



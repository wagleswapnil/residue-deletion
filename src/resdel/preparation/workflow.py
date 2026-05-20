from .peptide_builder import PeptideSystemBuilder

def run_prepare_workflow(config):
    sequence = config.system.sequence
    residue_to_delete = config.system.residue_to_delete

    system_builder = PeptideSystemBuilder(sequence, residue_to_delete)
    system_builder.generate_wt_peptide_structure("full_peptide.pdb")
    system_builder.generate_mutant_structure("mutant.pdb")

    system_builder.generate_topology_from_structure("full_peptide.pdb")



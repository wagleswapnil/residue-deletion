# We are going to build model peptide systems using PeptideBuilder.
# Originally, PeptideBuilder was written by Wilke Lab at UT Austin.  
# But we have found a fork (Bio2byte :: PeptideBuilder) of the original GitHub repo,
# which we will be using here. It supports some additional functionalities from the
# original PeptideBuilder, such as, terminal residues and three letter amino acid codes.

import re
from Bio.PDB import PDBIO
import bio2byte.PeptideBuilder as PeptideBuilder
from bio2byte.PeptideBuilder import Geometry
from resdel.preparation.create_openMM_topology import create_receptor_system

class PeptideSystemBuilder:
    def __init__(self, sequence, residue_to_delete):
        self.sequence = sequence
        self.residue_to_delete = int(residue_to_delete)
        self.full_sequence, self.mutated_sequence = self._get_indexed_sequences()


    def generate_wt_peptide_structure(self, output_path):
        self.generate_structure_from_sequence(self.full_sequence, output_path)
        return


    def generate_mutant_structure(self, output_path):
        self.generate_structure_from_sequence(self.mutated_sequence, output_path)
        return

    def generate_topology_from_structure(self, structure_path):
        pmd_receptor_struct = create_receptor_system(structure_path)
        pmd_receptor_struct.save("output1.pdb", overwrite=True)
        pmd_receptor_struct.save("output1.top", overwrite=True)
        return


    def generate_structure_from_sequence(self, sequence, output_path):
        extended_sheet_PhiPsi = (-135., 135.)
        structure = None
        if sequence[0] == "ACE":
            structure = PeptideBuilder.initialize_ACE()
        else:
            geo = Geometry.geometry(sequence[0])
            geo.phi, geo.psi_im1 = extended_sheet_PhiPsi
            structure = PeptideBuilder.initialize_res(geo)
        
        for res in sequence[1:-1]:
            geo = Geometry.geometry(res)
            geo.phi, geo.psi_im1 = extended_sheet_PhiPsi
            PeptideBuilder.add_residue(structure, geo)

        if sequence[-1] == "NME":
            PeptideBuilder.add_terminal_NME(structure)
        else:
            geo = Geometry.geometry(sequence[-1])
            geo.phi, geo.psi_im1 = extended_sheet_PhiPsi
            structure = PeptideBuilder.initialize_res(structure, geo)

        pdbwriter = PDBIO()
        pdbwriter.set_structure(structure)
        pdbwriter.save(output_path)
        return



    def _get_indexed_sequences(self):
        full_sequence = None
        mutated_sequence = None
        if self.sequence.isalpha():
            full_sequence = list(self.sequence)
            mutated_sequence = list(self.sequence[:self.residue_to_delete] + self.sequence[self.residue_to_delete :])
        else:
            sequence = re.split(r'[^`\=-~!@#$%^&*()_+\[\]{};\'\\:"|<,./<>?]', self.sequence)
            full_sequence = [s for s in sequence]
            mutated_sequence = full_sequence[:self.residue_to_delete - 1] + full_sequence[self.residue_to_delete :]
        print(self.sequence)
        return full_sequence, mutated_sequence




    

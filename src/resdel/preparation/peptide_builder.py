# We are going to build model peptide systems using PeptideBuilder.
# Originally, PeptideBuilder was written by Wilke Lab at UT Austin.  
# But we have found a fork (Bio2byte :: PeptideBuilder) of the original GitHub repo,
# which we will be using here. It supports some additional functionalities from the
# original PeptideBuilder, such as, terminal residues and three letter amino acid codes.

import re
import subprocess
from resdel.topology.parser import Parser
from resdel.topology.section import Section
from resdel.topology.writer import Writer
from resdel.topology.formatter import GromacsFormatter
from Bio.PDB import PDBIO
import bio2byte.PeptideBuilder as PeptideBuilder
from bio2byte.PeptideBuilder import Geometry
from resdel.preparation.create_openMM_topology import create_receptor_system

class PeptideSystemBuilder:
    def __init__(self, sequence, residue_to_delete):
        self.sequence = sequence
        self.residue_to_delete = int(residue_to_delete)
        self.wt_sequence, self.mutant_sequence = self._get_indexed_sequences()


    def generate_wt_peptide_structure(self, output_path):
        self.generate_structure_from_sequence(self.wt_sequence, output_path)
        return


    def generate_mutant_structure(self, output_path):
        self.generate_structure_from_sequence(self.mutant_sequence, output_path)
        return

    def generate_topology_from_structure(self, input_structure, structure_path, topology_path):
        pmd_receptor_struct = create_receptor_system(input_structure)
        pmd_receptor_struct.save(structure_path, overwrite=True)
        pmd_receptor_struct.save(topology_path, overwrite=True)
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
        pdbwriter.save(str(output_path))
        return


    def _get_indexed_sequences(self):
        wt_sequence = None
        mutant_sequence = None
        if self.sequence.isalpha():
            wt_sequence = list(self.sequence)
            mutant_sequence = list(self.sequence[:self.residue_to_delete] + self.sequence[self.residue_to_delete :])
        else:
            sequence = re.split(r'[^`\=-~!@#$%^&*()_+\[\]{};\'\\:"|<,./<>?]', self.sequence)
            wt_sequence = [s for s in sequence]
            mutant_sequence = wt_sequence[:self.residue_to_delete - 1] + wt_sequence[self.residue_to_delete :]
        print(self.sequence)
        return wt_sequence, mutant_sequence


    def generate_vaccuum_structure_from_solvent_structure(self, in_PDB_path, out_PDB_path):
        resnames = ["HOH", "SOL", "WAT", "NA", "CL", "Na", "Cl", "K"]
        cmd = f"pdb_delresname -{','.join(resnames)} {in_PDB_path} > {out_PDB_path}"
        subprocess.run(cmd, shell=True, check=True)
        return

    def generate_vaccuum_topology_from_solvent_topology(self, in_topology_path, out_topology_path):
        resnames = ["HOH", "SOL", "WAT", "NA", "CL", "Na", "Cl", "K"]
        parser = Parser(in_topology_path)
        topology = parser.parse_topology()

        topology.replace_tail_section_by_name("molecules", self.updated_molecules_section(topology.get_tail_section_by_name("molecules"), resnames)) 

        writer = Writer(topology, out_topology_path, GromacsFormatter())
        writer.write_topology()
        return

    def updated_molecules_section(self, molecules_section, resnames):
        new_molecules_section = Section("molecules")
        for line in molecules_section.lines:
            if line.tokens:
                molecule = line.tokens[0]
                if molecule in resnames:
                    pass
                else:
                    new_molecules_section.add_line(line)
            else:
                new_molecules_section.add_line(line)
        return new_molecules_section





    

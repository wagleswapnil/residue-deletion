# We are going to build model peptide systems using PeptideBuilder.
# Originally, PeptideBuilder was written by Wilke Lab at UT Austin.  
# But we have found a fork (Bio2byte :: PeptideBuilder) of the original GitHub repo,
# which we will be using here. It supports some additional functionalities from the
# original PeptideBuilder, such as, terminal residues and three letter amino acid codes.

import re
from typing import Optional
from Bio.PDB import PDBIO
import bio2byte.PeptideBuilder as PeptideBuilder
from bio2byte.PeptideBuilder import Geometry
from resdel.topology.parser import Parser, Line
from resdel.topology.writer import Writer
from resdel.preparation.create_openMM_topology import create_receptor_system
from resdel.topology.formatter import GromacsFormatter

class PeptideSystemBuilder:
    def __init__(self, sequence, residue_to_delete):
        self.sequence = sequence
        self.residue_to_delete = int(residue_to_delete)
        self.wt_sequence, self.mutant_sequence = self._get_indexed_sequences()

    def generate_wt_peptide_structure(self, output_path):
        self._generate_structure_from_sequence(self.wt_sequence, output_path)
        return

    def generate_mutant_structure(self, output_path):
        self._generate_structure_from_sequence(self.mutant_sequence, output_path)
        return

    def generate_topology_from_structure(self, input_structure, structure_path, topology_path):
        pmd_receptor_struct = create_receptor_system(input_structure)
        pmd_receptor_struct.save(structure_path, overwrite=True)
        pmd_receptor_struct.save(topology_path, overwrite=True)
        return

    def write_topology(self, topology, out_path):
        topology_writer = Writer(topology, out_path, formatter=GromacsFormatter())
        topology_writer.write_topology()
        return

    def add_posre_section_to_topology(self, in_path, out_path: Optional[str] = None, mol_name: Optional[str] = "system1"):
        if not out_path:
            out_path = in_path
        parser = Parser(in_path)
        topology = parser.parse_topology()
        
        for mol in topology.molecules:
            if mol.name == mol_name:
                line = Line(f"")
                mol.sections[-1].add_line(line)
                line = Line(f"#ifdef POSRES")
                mol.sections[-1].add_line(line)
                line = Line(f'#include "posre.itp"')
                mol.sections[-1].add_line(line)
                line = Line(f"#endif")
                mol.sections[-1].add_line(line)
                line = Line(f"")
                mol.sections[-1].add_line(line)
        self.write_topology(topology, out_path=out_path)
        return

    def _generate_structure_from_sequence(self, sequence, output_path):
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
            structure = PeptideBuilder.add_residue(structure, geo)
            PeptideBuilder.add_terminal_OXT(structure)

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

    





    

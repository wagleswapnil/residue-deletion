from resdel.topology import *
from typing import List, Optional

from resdel.topology.atom import Atom

class Transformation_Bonded_Params:
    def __init__(self, topA : Topology, topB : Topology, res_to_delete : int, molA_name : Optional[str] = "system1", molB_name : Optional[str] = "system1"):
        self.topA = topA
        self.topB = topB
        self.molA = None
        self.molB = None
        self.molA_name = molA_name
        self.molB_name = molB_name
        self.res_to_delete = res_to_delete
        for mol in self.topA.molecules:
            if mol.name == self.molA_name:
                self.molA = mol
                break

        for mol in self.topB.molecules:
            if mol.name == self.molB_name:
                self.molB = mol
                break
        if self.molA is None:
            raise ValueError(f"Could not find molecule with name {self.molA_name} in topology A")
        if self.molB is None:
            raise ValueError(f"Could not find molecule with name {self.molB_name} in topology B")
    
    
    def atom_key(self, tokens, top_name):
        resnr = tokens[2]
        residue = tokens[3]
        atom = tokens[4]
        if top_name == "B" and int(resnr) >= self.res_to_delete:
            resnr = str(int(resnr) + 1)
        return (resnr, residue, atom)


    def build_atom_mapping(self):
        atomsA_dict = {}
        atomsB_dict = {}
        self.mapping = {}

        for i, line in enumerate(self.molA.get_section("atoms").lines):
            if line.tokens:
                atomsA_dict[self.atom_key(line.tokens, "A")] = line.tokens[0]

        for i, line in enumerate(self.molB.get_section("atoms").lines):
            if line.tokens:
                atomsB_dict[self.atom_key(line.tokens, "B")] = line.tokens[0]
        
        for key in atomsB_dict:
            if key in atomsA_dict:
                self.mapping[atomsB_dict[key]] = atomsA_dict[key]
            else:
                raise ValueError(f"Warning: Atom {key} in topology B not found in topology A. This is not supposed to happen.")
        return
    
    def canonicalize_bond(self, atom1, atom2):
        if int(atom1) < int(atom2):
            return (atom1, atom2)
        else:
            return (atom2, atom1)

    def tranform_bonds(self):
        for i, line in enumerate(self.molA.get_section("bonds").lines):
            if line.tokens:
                atom1 = self.mapping[line.tokens[0]]
                atom2 = self.mapping[line.tokens[1]]
                canonical_bond = self.canonicalize_bond(atom1, atom2)
                line.tokens[0] = canonical_bond[0]
                line.tokens[1] = canonical_bond[1]
            breakpoint()

        for i, line in enumerate(self.molB.get_section("bonds").lines):
            if line.tokens:
                if line.tokens[0] in self.mapping and line.tokens[1] in self.mapping:
                    atom1 = self.mapping[line.tokens[0]]
                    atom2 = self.mapping[line.tokens[1]]
                    canonical_bond = self.canonicalize_bond(atom1, atom2)
                    line.tokens[0] = canonical_bond[0]
                    line.tokens[1] = canonical_bond[1]
                else:
                    raise ValueError(f"Warning: Bond between atoms {line.tokens[0]} and {line.tokens[1]} in topology B cannot be mapped to topology A. This is not supposed to happen.")
            
        return
    

    


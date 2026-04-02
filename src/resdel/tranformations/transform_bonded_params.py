from resdel.topology import *
from typing import List, Optional

from resdel.topology.atom import Atom

class Transformation_Bonded_Params:
    def __init__(self, topA : Topology, topB : Topology, molA_name : Optional[str] = "system1", molB_name : Optional[str] = "system1"):
        self.topA = topA
        self.topB = topB
        self.molA = None
        self.molB = None
        self.molA_name = molA_name
        self.molB_name = molB_name
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

    def reassign_atom_indices(self):
        atoms : List[Atom] = []
        for i, line in enumerate(self.molA.get_section("atoms").lines):
            if line.tokens:
                atom = Atom()
                atom.get_atom_info(line.tokens)
                atoms.append(atom)
                breakpoint()
        return atoms

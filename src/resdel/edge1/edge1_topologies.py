from resdel.topology import *
from typing import Optional, List
from resdel.tranformations.atom_mapping import build_atom_mapping

class Edge1_Topologies:
    def __init__(self, topA : Topology, topB : Topology, residue_to_delete : str, molA_name : Optional[str] = "system1", molB_name : Optional[str] = "system1"):
        self.topA = topA
        self.topB = topB
        self.molA = None
        self.molB = None
        self.molA_name = molA_name
        self.molB_name = molB_name
        self.residue_to_delete = residue_to_delete
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
        self.mapping = build_atom_mapping(self.molA.get_section("atoms").lines, self.molB.get_section("atoms").lines, int(self.residue_to_delete))

    def get_exclusions(self):
        exclusionsA = self.molA.get_section("exclusions")
        exclusionsB = self.molB.get_section("exclusions")
        breakpoint()
        return exclusionsA, exclusionsB


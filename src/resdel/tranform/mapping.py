from resdel.topology.topology import Topology
from typing import Optional

class AtomMapping:
    def __init__(self, topA: Topology, topB: Topology, mol_nameA : Optional[str] = "system1", mol_nameB : Optional[str] = "system1"):
        self.topA = topA
        self.topB = topB
        self.mol_nameA = mol_nameA
        self.mol_nameB = mol_nameB
        self.molA = None
        self.molB = None
        
        for mol in self.topA.molecules:
            if mol.name == self.mol_nameA:
                self.molA = mol
                break
        for mol in self.topB.molecules:
            if mol.name == self.mol_nameB:
                self.molB = mol
                break
        if self.molA is None:
            raise ValueError(f"Could not find molecule with name {self.mol_nameA} in topology A")
        if self.molB is None:   
            raise ValueError(f"Could not find molecule with name {self.mol_nameB} in topology B")

    def build_atom_mapping(self):
        """
        Build a mapping of atoms between two topologies based on their names, types and residue numbers.
        
        Args:
            topA (Topology): The first topology (i-1_i_i+1).
            topB (Topology): The second topology (i-1_i+1).
        """
        

        atomsA = self.topA.molecules[0].get_section("atoms").lines
        atomsB = self.topB.molecules[0].get_section("atoms").lines
        breakpoint()
        return
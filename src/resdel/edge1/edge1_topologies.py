from resdel.topology import *
from typing import Optional, List
from resdel.tranformations.atom_mapping import build_atom_mapping
import subprocess
import re

class Edge1_Topologies:
    def __init__(self, topA : Topology, topB : Topology, residue_to_delete : str, edge_steps : int, molA_name : Optional[str] = "system1", molB_name : Optional[str] = "system1"):
        self.topA = topA
        self.topB = topB
        self.edge_steps = edge_steps
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
        self.get_tpr_dump(mdp="./tests/MDP/em.mdp", structure="./tests/data/minimized.gro", topology="./tests/data/system_stage1.top", output_prefix="./tests/data/system_stage1")
        output_prefix="./tests/data/system_stage1"
        self.exclusionsA = self.extract_exclusions_from_tpr_dump(f"{output_prefix}.txt", f"{output_prefix}_exclusions.txt")
        self.pairsA = self.extract_pairs_from_topology()
        assert(self.pairsA.issubset(self.exclusionsA)), "Error: Not all pairs in topology A are present in the exclusions extracted from the tpr dump. This is not supposed to happen."


    def get_tpr_dump(self, mdp, structure, topology, output_prefix):
        cmd =["gmx", "grompp", "-f", mdp, "-c", structure, "-p", topology, "-o", f"{output_prefix}.tpr"]
        #result = subprocess.run(cmd, check=True, cwd="./")
        cmd = f"gmx dump -s {output_prefix}.tpr > {output_prefix}.txt"
        #result = subprocess.run(cmd, shell=True, check=True)
        return
    
    def extract_exclusions_from_tpr_dump(self, tpr_dump_file, output_file, molecule_name : Optional[str] = "system1"):
        exclusions = {}
        in_target_moltype = False
        read_exclusions = False
        buffer = None
        collecting = False
        f = open(tpr_dump_file, "r")
        for line in f:
            if line.strip().startswith("moltype"):
                in_target_moltype = False

            if f'name="{molecule_name}"' in line:
                in_target_moltype = True
                continue

            if in_target_moltype:
                if "Bond:" in line:
                    read_exclusions = False
                elif "excls:" in line:
                    read_exclusions = True
            
            if in_target_moltype and read_exclusions:
                if "numLists" in line or "numElements" in line:
                    continue
                if line.strip().startswith("excls["):
                    buffer = line.strip()
                    collecting = True
                    if "}" in line:
                        atom, nums = self.parse_excls_buffer(buffer)
                        exclusions[atom] = nums
                        buffer = None
                        collecting = False
                elif collecting:
                    buffer += " " + line.strip()
                    if "}" in line:
                        atom, nums = self.parse_excls_buffer(buffer)
                        exclusions[atom] = nums
                        buffer = None
                        collecting = False
        f.close()
        return self.make_exclusions_set(exclusions)
    
    def parse_excls_buffer(self, buffer):
        pattern = r"excls\[(\d+)\]\[num=\d+\]=\{([^}]*)\}"
        match = re.search(pattern, buffer)
        atom_index = self.make_exclusions_1_indexed(match.group(1))
        nums_str = self.make_exclusions_1_indexed(match.group(2).strip().split(','))
        return  atom_index, nums_str 
    
    
    def make_exclusions_1_indexed(self, num):
        if isinstance(num, list):
            return [str(int(n) + 1) for n in num]
        return str(int(num) + 1)

    def make_exclusions_set(self, exclusions_dict):
        exclusions_set = set()
        for key, value in exclusions_dict.items():
            print(f"Processing atom {key} with exclusions {value}")
            for num in value:
                pair = tuple(sorted([int(key), int(num)]))
                exclusions_set.add(pair)
        print(f"Extracted exclusions: {exclusions_set}")
        return exclusions_set

    def extract_pairs_from_topology(self):
        for mol in self.topA.molecules:
            if mol.name == self.molA_name:
                pairs_section = mol.get_section("pairs")
                if pairs_section is not None:
                    print(f"Extracted pairs: {self.generate_pairs_set(pairs_section.lines)}")
                    return self.generate_pairs_set(pairs_section.lines)
                else:
                    print(f"No pairs section found in molecule {self.molA_name} of topology A.")
        return

    def generate_pairs_set(self, pairs_lines):
        pairs_set = set()
        for line in pairs_lines:
            if line.tokens:
                pairs_set.add(tuple(sorted([int(line.tokens[0]), int(line.tokens[1])])))
        return pairs_set
        

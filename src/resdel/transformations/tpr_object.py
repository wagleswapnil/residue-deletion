from typing import Optional
import subprocess
import re

class TPR_Object:
    def __init__(self, topology, structure, mdp : Optional[str] = "./tests/MDP/em.mdp"):
        self.topology = topology
        self.structure = structure
        self.mdp = mdp

    def make_tpr_dump(self, output_prefix):
        self.output_prefix = output_prefix
        cmd =["gmx", "grompp", "-f", self.mdp, "-c", self.structure, "-p", self.topology, "-o", f"{output_prefix}.tpr"]
        subprocess.run(cmd, check=True, cwd="./")
        cmd = f"gmx dump -s {output_prefix}.tpr > {output_prefix}.txt"
        subprocess.run(cmd, shell=True, check=True)
        return
    
    def extract_exclusions_from_tpr_dump(self, molecule_name : Optional[str] = "system1"):
        exclusions = {}
        in_target_moltype = False
        read_exclusions = False
        buffer = None
        collecting = False
        tpr_dump_file = f"{self.output_prefix}.txt"
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
                        atom, nums = self._parse_excls_buffer(buffer)
                        exclusions[atom] = nums
                        buffer = None
                        collecting = False
                elif collecting:
                    buffer += " " + line.strip()
                    if "}" in line:
                        atom, nums = self._parse_excls_buffer(buffer)
                        exclusions[atom] = nums
                        buffer = None
                        collecting = False
        f.close()
        return self._make_exclusions_set(exclusions)

    def _parse_excls_buffer(self, buffer):
        pattern = r"excls\[(\d+)\]\[num=\d+\]=\{([^}]*)\}"
        match = re.search(pattern, buffer)
        atom_index = self._make_exclusions_1_indexed(match.group(1))
        nums_str = self._make_exclusions_1_indexed(match.group(2).strip().split(','))
        return  atom_index, nums_str 


    def _make_exclusions_1_indexed(self, num):
        if isinstance(num, list):
            return [str(int(n) + 1) for n in num]
        return str(int(num) + 1)

    def _make_exclusions_set(self, exclusions_dict):
        exclusions_set = set()
        for key, value in exclusions_dict.items():
            #print(f"Processing atom {key} with exclusions {value}")
            for num in value:
                pair = tuple(sorted([int(key), int(num)]))
                exclusions_set.add(pair)
        #print(f"Extracted exclusions: {exclusions_set}")
        return exclusions_set

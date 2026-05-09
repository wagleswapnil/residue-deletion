from typing import Optional
import re

def get_residue_atom_idxs(mol, residue_to_delete):
    idx_i_minus_1 = []
    idx_i = []
    idx_i_plus_1 = []
    idx_i_minus_1_C = None
    idx_i_plus_1_N = None
    idx_i_C = None
    idx_i_N = None
    idx_i_minus_1_N = None
    idx_i_plus_1_C = None
    for line in mol.get_section("atoms").lines:
        if line.tokens:
            atom_idx = line.tokens[0]
            resnr = line.tokens[2]
            atom_name = line.tokens[4]
            if resnr == residue_to_delete:
                idx_i.append(int(atom_idx))
                if atom_name == "C":
                    idx_i_C = int(atom_idx)
                elif atom_name == "N":
                    idx_i_N = int(atom_idx)
            elif int(resnr) == int(residue_to_delete) - 1:
                idx_i_minus_1.append(int(atom_idx))
                if atom_name == "C":
                    idx_i_minus_1_C = int(atom_idx)
                elif atom_name == "N":
                    idx_i_minus_1_N = int(atom_idx)
            elif int(resnr) == int(residue_to_delete) + 1:
                idx_i_plus_1.append(int(atom_idx))
                if atom_name == "N":
                    idx_i_plus_1_N = int(atom_idx)
                elif atom_name == "C":
                    idx_i_plus_1_C = int(atom_idx)
    if idx_i_minus_1_C is None or idx_i_plus_1_N is None or idx_i_C is None or idx_i_N is None or idx_i_minus_1_N is None or idx_i_plus_1_C is None:
        raise ValueError("Could not find all required atoms in the specified residues. This is not supposed to happen.")
    return idx_i_minus_1, idx_i, idx_i_plus_1, idx_i_minus_1_N, idx_i_minus_1_C, idx_i_N, idx_i_C, idx_i_plus_1_N,  idx_i_plus_1_C


def get_sigma_epsilon_charges(top_atomtypes, mol, idx_list):
    sigma_epsilon_charges = {}
    atomtypes = {}
    for line in top_atomtypes.lines:
        if line.tokens:
            atom_type = line.tokens[0]
            sigma = line.tokens[5]
            epsilon = line.tokens[6]
            atomtypes[atom_type] = (sigma, epsilon)
    
    for line in mol.get_section("atoms").lines:
        if line.tokens:
            atom_idx = line.tokens[0]
            atom_type = line.tokens[1]
            partial_charge = line.tokens[6]
            if int(atom_idx) in idx_list:
                sigma_epsilon_charges[atom_idx] = (atomtypes[atom_type][0] if atom_type in atomtypes else None, atomtypes[atom_type][1] if atom_type in atomtypes else None, partial_charge)
    return sigma_epsilon_charges


def extract_exclusions_from_tpr_dump(tpr_dump_file, output_file, molecule_name : Optional[str] = "system1"):
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
                    atom, nums = parse_excls_buffer(buffer)
                    exclusions[atom] = nums
                    buffer = None
                    collecting = False
            elif collecting:
                buffer += " " + line.strip()
                if "}" in line:
                    atom, nums = parse_excls_buffer(buffer)
                    exclusions[atom] = nums
                    buffer = None
                    collecting = False
    f.close()
    return make_exclusions_set(exclusions)

def parse_excls_buffer(buffer):
    pattern = r"excls\[(\d+)\]\[num=\d+\]=\{([^}]*)\}"
    match = re.search(pattern, buffer)
    atom_index = make_exclusions_1_indexed(match.group(1))
    nums_str = make_exclusions_1_indexed(match.group(2).strip().split(','))
    return  atom_index, nums_str 


def make_exclusions_1_indexed(num):
    if isinstance(num, list):
        return [str(int(n) + 1) for n in num]
    return str(int(num) + 1)

def make_exclusions_set(exclusions_dict):
    exclusions_set = set()
    for key, value in exclusions_dict.items():
        #print(f"Processing atom {key} with exclusions {value}")
        for num in value:
            pair = tuple(sorted([int(key), int(num)]))
            exclusions_set.add(pair)
    #print(f"Extracted exclusions: {exclusions_set}")
    return exclusions_set

def map_exclusions_pairs(pairsB, mapping):
    mapped_pairs = set()
    for pair in pairsB:
        #breakpoint()
        mapped_pair = tuple(sorted([mapping[pair[0]], mapping[pair[1]]]))
        mapped_pairs.add(mapped_pair)
    return mapped_pairs


def get_tpr_dump(mdp, structure, topology, output_prefix):
    cmd =["gmx", "grompp", "-f", mdp, "-c", structure, "-p", topology, "-o", f"{output_prefix}.tpr"]
    #result = subprocess.run(cmd, check=True, cwd="./")
    cmd = f"gmx dump -s {output_prefix}.tpr > {output_prefix}.txt"
    #result = subprocess.run(cmd, shell=True, check=True)
    return


def extract_pairs_from_topology(top, molecule_name):
    for mol in top.molecules:
        if mol.name == molecule_name:
            pairs_section = mol.get_section("pairs")
            if pairs_section is not None:
                #print(f"Extracted pairs: {self.generate_pairs_set(pairs_section.lines)}")
                return generate_pairs_set(pairs_section.lines)
            else:
                print(f"No pairs section found in molecule {molecule_name} of topology.")
    return

def generate_pairs_set(pairs_lines):
    pairs_set = set()
    for line in pairs_lines:
        if line.tokens:
            pairs_set.add(tuple(sorted([int(line.tokens[0]), int(line.tokens[1])])))
    return pairs_set
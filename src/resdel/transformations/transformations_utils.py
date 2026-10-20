from typing import Optional


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
            if int(resnr) == residue_to_delete:
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
    #breakpoint()
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



def map_exclusions_pairs(pairsB, mapping):
    mapped_pairs = set()
    for pair in pairsB:
        #breakpoint()
        mapped_pair = tuple(sorted([mapping[pair[0]], mapping[pair[1]]]))
        mapped_pairs.add(mapped_pair)
    return mapped_pairs

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
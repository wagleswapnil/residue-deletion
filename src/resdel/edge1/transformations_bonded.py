from resdel.topology import *
from typing import Optional, List

def add_peptide_bond(mol, idx_i_minus_1_C, idx_i_plus_1_N):
    line = Line(f"#ifdef NEW_PEPTIDE_BOND")
    mol.get_section("bonds").add_line(line)
    line = Line(f"\t{idx_i_minus_1_C} \t{idx_i_plus_1_N} \t 5")
    mol.get_section("bonds").add_line(line)
    line = Line(f"#endif")
    mol.get_section("bonds").add_line(line)
    line = Line(f"")
    mol.get_section("bonds").add_line(line)
    return

def extract_harmonic_bonds_from_topology(molecule):
    bonds_section = molecule.get_section("bonds")
    bonds_set = set()
    for line in bonds_section.lines:
        if line.tokens:
            if line.tokens[0].startswith("#"):
                continue
            else:
                idx1, idx2, bond_type = int(line.tokens[0]), int(line.tokens[1]), int(line.tokens[2])
                if bond_type == 1:  # Harmonic bond type
                    bonds_set.add(tuple(sorted((idx1, idx2))))
    return sorted(bonds_set)

def updated_bonds_section(bonds_lines, topB_bond_to_add, idx_i_minus_1_N, idx_i_minus_1_C, idx_i_N, idx_i_C, idx_i_plus_1_N,  idx_i_plus_1_C):
    new_bonds_section = Section("bonds")
    bond_i_minus_1_i = False
    bond_i_i_plus_1 = False
    bond_i_minus_1_i_plus_1 = False
    
    for i, line in enumerate(bonds_lines):
        if line.tokens:
            if line.tokens[0].startswith("#"):
                new_bonds_section.add_line(line)
            else:
                idx1, idx2 = sorted([int(line.tokens[0]), int(line.tokens[1])])
                if (idx1, idx2) == (idx_i_minus_1_C, idx_i_N):
                    new_bonds_section.add_line(Line(f"{line.raw} ; i-1 to i bond"))
                    bond_i_minus_1_i = True
                elif (idx1, idx2) == (idx_i_C, idx_i_plus_1_N):
                    b0, DA = line.tokens[3], line.tokens[4]
                    DB = 0
                    betaA, betaB = 1, 0
                    new_bonds_section.add_line(Line(f"#ifdef OLD_PEPTIDE_BOND"))
                    new_bonds_section.add_line(Line(f"{line.raw} ; i to i+1 bond"))
                    new_bonds_section.add_line(Line(f"#endif"))
                    new_bonds_section.add_line(Line(f"#ifdef NEW_MORSE_BOND"))
                    new_bonds_section.add_line(Line(f"{idx1}  {idx2}  3  {b0} {DA}  {betaA}  {b0} {DB} {betaB}; i to i+1 bond")) 
                    new_bonds_section.add_line(Line(f"#endif"))
                    bond_i_i_plus_1 = True
        else:
            new_bonds_section.add_line(line)
    return new_bonds_section

def map_bonds(bonds, mapping):
    mapped_bonds = set()
    for bond in bonds:
        idx1, idx2 = bond
        if idx1 in mapping and idx2 in mapping:
            mapped_bonds.add(tuple(sorted((mapping[idx1], mapping[idx2]))))
        else:
            print(f"Warning: Bond ({idx1}, {idx2}) cannot be mapped because one of the indices is missing in the mapping. This bond will be skipped.")
    return mapped_bonds


def get_topB_bond_parameters(molB, bondsB_minus_A, mapping):
    topB_bond_parameters = None
    for line in molB.get_section("bonds").lines:
        if line.tokens:
            if line.tokens[0].startswith("#"):
                continue
            else:
                idx1, idx2= int(mapping[int(line.tokens[0])]), int(mapping[int(line.tokens[1])])
                if tuple(sorted((idx1, idx2))) in bondsB_minus_A:
                    topB_bond_parameters = line
                    break
    if topB_bond_parameters is None:
        raise ValueError("Could not find the parameters for the new bond in topology B. This is not supposed to happen.")
    return topB_bond_parameters
    
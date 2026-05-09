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

def extract_bonds_from_topology(molecule):
    bonds_section = molecule.get_section("bonds")
    bonds_set = set()
    for line in bonds_section.lines:
        if line.tokens:
            if line.tokens[0].startswith("#"):
                break
            else:
                idx1, idx2 = int(line.tokens[0]), int(line.tokens[1])
                bonds_set.add(tuple(sorted((idx1, idx2))))
    return sorted(bonds_set)

def updated_bonds_section(bonds_lines, idx_i_minus_1_N, idx_i_minus_1_C, idx_i_N, idx_i_C, idx_i_plus_1_N,  idx_i_plus_1_C):
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
from resdel.topology import *

def add_peptide_bond(mol, idx_i_minus_1_C, idx_i_plus_1_N):
    line = Line(f"#ifdef ADD_PEPTIDE_BOND")
    mol.get_section("bonds").add_line(line)
    line = Line(f"\t{idx_i_minus_1_C} \t{idx_i_plus_1_N} \t 5")
    mol.get_section("bonds").add_line(line)
    line = Line(f"#endif\n")
    mol.get_section("bonds").add_line(line)
    return
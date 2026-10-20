from resdel.topology import *
from typing import Optional
from .calculate_bond_distances import get_min_max_distances

def add_fc0_bond(mol, idx_i_minus_1_C, idx_i_plus_1_N):
    line = Line(f"#ifdef NEW_FC0_BOND")
    mol.get_section("bonds").add_line(line)
    line = Line(f"\t{idx_i_minus_1_C} \t{idx_i_plus_1_N} \t 5")
    mol.get_section("bonds").add_line(line)
    line = Line(f"#endif")
    mol.get_section("bonds").add_line(line)
    return

def add_new_bond(mol, top_mutant_bond_to_add):
    restr_ftype = "10"
    restr_rmax = None
    if top_mutant_bond_to_add is not None:
        idx1, idx2, r0, fc = top_mutant_bond_to_add.tokens[0], top_mutant_bond_to_add.tokens[1], top_mutant_bond_to_add.tokens[3], top_mutant_bond_to_add.tokens[4]
        line = Line(f"#ifdef NEW_DIST_RESTR")
        mol.get_section("bonds").add_line(line)
        #min_d, max_d = get_min_max_distances(idx1, idx2)
        max_d = 0.3
        min_d = 0.3       
        restr_rmax = 1.65
        line = Line(f"\t{idx1}  {idx2} {restr_ftype}  {min_d}  {max_d}  {restr_rmax}  {float(fc) * 0.1}  {r0} {r0}  {restr_rmax}  {fc} ; Distance restraint for the new bond between atoms {idx1} and {idx2}")
        mol.get_section("bonds").add_line(line)
        line = Line(f"#endif")
        mol.get_section("bonds").add_line(line)
        line = Line(f"#ifdef NEW_PEPTIDE_BOND")
        mol.get_section("bonds").add_line(line)
        line = Line(f"\t{idx1}  {idx2} 1  {r0}  {fc} ; i-1 to i+1 peptide bond")
        mol.get_section("bonds").add_line(line)
        line = Line(f"#endif")
        mol.get_section("bonds").add_line(line)
        #line = Line(f"")
        #mol.get_section("bonds").add_line(line)
    else:
        print("Warning: No new bond parameters provided. Distance restraints for the new bond will not be added.")
    return

def _get_edge1_distance_restraint_line(idx1, idx2, fc, restr_ftype, restr_rmax, restr_rmin, lambda_value):
    fc = fc * lambda_value * 0.1
    line = Line(f"\t{idx1}  {idx2} {restr_ftype}  {restr_rmin}  {restr_rmax} 1.65 {fc}")
    return line

def add_edge1_distance_restraints(mol, top_mutant_bond_to_add, edge1_steps, edge1_lambda_vector):
    restr_ftype = "10"
    restr_rmax = 0.3
    restr_rmin = 0.3
    idx1, idx2, r0, fc = top_mutant_bond_to_add.tokens[0], top_mutant_bond_to_add.tokens[1], top_mutant_bond_to_add.tokens[3], top_mutant_bond_to_add.tokens[4]
    idx1, idx2, r0, fc = int(idx1), int(idx2), float(r0), float(fc)
    for step in range(0, len(edge1_lambda_vector)):
        line = Line(f"#ifdef EDGE1_{step}")
        mol.get_section("bonds").add_line(line)
        line = _get_edge1_distance_restraint_line(idx1, idx2, fc, restr_ftype, restr_rmax, restr_rmin, edge1_lambda_vector[step])
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

def updated_bonds_section(bonds_section, top_mutant_bond_to_add, idx_i, idx_i_minus_1_C, idx_i_N, idx_i_C, idx_i_plus_1_N):
    new_bonds_section = Section("bonds")
    bond_i_minus_1_i = False
    bond_i_i_plus_1 = False
    for i, line in enumerate(bonds_section.lines):
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
                elif idx1 in idx_i and idx2 in idx_i:
                    new_bonds_section.add_line(Line(f"{line.raw} ; internal residue i bonds"))
                else:
                    new_bonds_section.add_line(line)
        elif line.raw.strip() == "":
            pass
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


def get_top_mutant_bond_parameters(mol_mutant, bonds_mutant_minus_wt, mapping):
    top_mutant_bond_parameters = None
    for line in mol_mutant.get_section("bonds").lines:
        if line.tokens:
            if line.tokens[0].startswith("#"):
                continue
            else:
                idx1, idx2= int(mapping[int(line.tokens[0])]), int(mapping[int(line.tokens[1])])
                if tuple(sorted((idx1, idx2))) in bonds_mutant_minus_wt:
                    top_mutant_bond_parameters = Line(f"{idx1}  {idx2}  {' '.join(line.tokens[2:])}")
                    break
    if top_mutant_bond_parameters is None:
        raise ValueError("Could not find the parameters for the new bond in mutant topology. This is not supposed to happen.")
    return top_mutant_bond_parameters

def extract_angles_from_topology(molecule, mapping : Optional[dict] = None):
    angles_section = molecule.get_section("angles")
    angles_set = set()
    for line in angles_section.lines:
        if line.tokens:
            idx1, idx2, idx3 = int(line.tokens[0]), int(line.tokens[1]), int(line.tokens[2]) 
            if mapping:
                idx1, idx2, idx3 = mapping[idx1], mapping[idx2], mapping[idx3]
            if idx1 <= idx3:  # To avoid duplicates like (1,2,3) and (3,2,1)
                angles_set.add(tuple((idx1, idx2, idx3)))
            else:
                angles_set.add(tuple((idx3, idx2, idx1)))
    return angles_set

def updated_angles_section(angles_section, idx_i, top_mutant_angles_to_add):
    new_angles_section = Section("angles")
    residue_i_internal_angles = []
    angles_involving_residue_i = []
    for line in angles_section.lines:
        if line.tokens:
            idx1, idx2, idx3 = int(line.tokens[0]), int(line.tokens[1]), int(line.tokens[2])
            if idx1 in idx_i and idx2 in idx_i and idx3 in idx_i:
                residue_i_internal_angles.append(line)
            elif idx1 in idx_i or idx2 in idx_i or idx3 in idx_i:
                angles_involving_residue_i.append(line)
            else:
                new_angles_section.add_line(line)
        elif line.raw.strip() == "":
            pass
        else:
            new_angles_section.add_line(line)

    new_angles_section.add_line(Line(f"#ifdef INTERNAL_I_ANGLES"))
    for line in residue_i_internal_angles:
        new_angles_section.add_line(Line(f"{line.raw} ; internal angle for residue i"))
    new_angles_section.add_line(Line(f"#endif"))

    new_angles_section.add_line(Line(f"#ifdef DUAL_TOP_ANGLES"))
    for line in angles_involving_residue_i:
        idx1, idx2, idx3, ftype, theta0, fc = int(line.tokens[0]), int(line.tokens[1]), int(line.tokens[2]), line.tokens[3], line.tokens[4], line.tokens[5]
        new_angles_section.add_line(Line(f"\t{idx1}  {idx2}  {idx3}  {ftype}  {theta0}  {fc}  {theta0}  0.0 ; turning off angles involving residue i"))
    for line in top_mutant_angles_to_add:
        idx1, idx2, idx3, ftype, theta0, fc = int(line.tokens[0]), int(line.tokens[1]), int(line.tokens[2]), line.tokens[3], line.tokens[4], line.tokens[5]
        new_angles_section.add_line(Line(f"\t{idx1}  {idx2}  {idx3}  {ftype}  {theta0}  0.0  {theta0}   {fc} ; turning on Top B angles"))
    new_angles_section.add_line(Line(f"#endif"))

    new_angles_section.add_line(Line(f"#ifdef NEW_ANGLES"))
    for line in top_mutant_angles_to_add:
        new_angles_section.add_line(Line(f"{line.raw} ; topB angles"))
    new_angles_section.add_line(Line(f"#endif"))


    new_angles_section.add_line(Line(f"#ifdef ANGLES_OF_RES_I"))
    for line in angles_involving_residue_i:
        new_angles_section.add_line(Line(f"\t{line.raw} ; angle involving residue i"))
    new_angles_section.add_line(Line(f"#endif"))
    new_angles_section.add_line(Line(f""))
    return new_angles_section

def get_top_mutant_angles_parameters(mol_mutant, angles_mutant_minus_wt, mapping):
    mol_mutant_angles_to_add = []
    for angle in angles_mutant_minus_wt:
        for line in mol_mutant.get_section("angles").lines:
            if line.tokens:
                idx1, idx2, idx3, ftype, theta0, fc = int(mapping[int(line.tokens[0])]), int(mapping[int(line.tokens[1])]), int(mapping[int(line.tokens[2])]), line.tokens[3], line.tokens[4], line.tokens[5]
                if tuple((idx1, idx2, idx3)) == angle or tuple((idx3, idx2, idx1)) == angle:
                    mol_mutant_angles_to_add.append(Line(f"\t{idx1}  {idx2}  {idx3}  {ftype}  {theta0}  {fc}"))
                    break
    return mol_mutant_angles_to_add

def updated_dihedrals_section(dihedrals_section, idx_i, top_mutant_dihedrals_to_add):
    new_dihedrals_section = Section("dihedrals")
    residue_i_internal_dihedrals = []
    dihedrals_involving_residue_i = []
    for line in dihedrals_section.lines:
        if line.tokens:
            idx1, idx2, idx3, idx4 = int(line.tokens[0]), int(line.tokens[1]), int(line.tokens[2]), int(line.tokens[3])
            if idx1 in idx_i and idx2 in idx_i and idx3 in idx_i and idx4 in idx_i:
                residue_i_internal_dihedrals.append(line)
            elif idx1 in idx_i or idx2 in idx_i or idx3 in idx_i or idx4 in idx_i:
                dihedrals_involving_residue_i.append(line)
            else:
                new_dihedrals_section.add_line(line)
        elif line.raw.strip() == "":
            pass
        else:
            new_dihedrals_section.add_line(line)

    new_dihedrals_section.add_line(Line(f"#ifdef INTERNAL_I_DIHEDRALS"))
    for line in residue_i_internal_dihedrals:
        new_dihedrals_section.add_line(Line(f"{line.raw} ; internal dihedral for residue i"))
    new_dihedrals_section.add_line(Line(f"#endif"))

    new_dihedrals_section.add_line(Line(f"#ifdef DUAL_TOP_DIHEDRALS"))
    for line in dihedrals_involving_residue_i:
        idx1, idx2, idx3, idx4, ftype, phi0, fc, multiplicity = int(line.tokens[0]), int(line.tokens[1]), int(line.tokens[2]), int(line.tokens[3]), line.tokens[4], line.tokens[5], line.tokens[6], line.tokens[7]
        new_dihedrals_section.add_line(Line(f"\t{idx1}  {idx2}  {idx3}  {idx4}  {ftype}  {phi0}  {fc}  {multiplicity} {phi0}  0.0  {multiplicity} ; turning off dihedral involving residue i"))
    for line in top_mutant_dihedrals_to_add:
        idx1, idx2, idx3, idx4, ftype, phi0, fc, multiplicity = int(line.tokens[0]), int(line.tokens[1]), int(line.tokens[2]), int(line.tokens[3]), line.tokens[4], line.tokens[5], line.tokens[6], line.tokens[7]
        new_dihedrals_section.add_line(Line(f"\t{idx1}  {idx2}  {idx3}  {idx4}  {ftype}  {phi0}  0.0  {multiplicity} {phi0}  {fc}  {multiplicity}  ; turning on Top B dihedral"))
    new_dihedrals_section.add_line(Line(f"#endif"))

    new_dihedrals_section.add_line(Line(f"#ifdef NEW_DIHEDRALS"))
    for line in top_mutant_dihedrals_to_add:
        idx1, idx2, idx3, idx4, ftype, phi0, fc, multiplicity = int(line.tokens[0]), int(line.tokens[1]), int(line.tokens[2]), int(line.tokens[3]), line.tokens[4], line.tokens[5], line.tokens[6], line.tokens[7]
        new_dihedrals_section.add_line(Line(f"\t{idx1}  {idx2}  {idx3}  {idx4}  {ftype}  {phi0}  {fc}  {multiplicity}  ; Top B dihedral"))
    new_dihedrals_section.add_line(Line(f"#endif"))

    new_dihedrals_section.add_line(Line(f"#ifdef DIHEDRALS_OF_RES_I"))
    for line in dihedrals_involving_residue_i:
        new_dihedrals_section.add_line(Line(f"\t{line.raw} ; dihedral involving residue i"))
    new_dihedrals_section.add_line(Line(f"#endif"))
    new_dihedrals_section.add_line(Line(f""))
    return new_dihedrals_section


def get_top_mutant_dihedrals_parameters(mol_mutant, dihedrals_mutant_minus_wt, mapping):
    mol_mutant_dihedrals_to_add = []
    for dihedral in dihedrals_mutant_minus_wt:
        for line in mol_mutant.get_section("dihedrals").lines:
            if line.tokens:
                idx1, idx2, idx3, idx4, ftype, phi0, fc, multiplicity = int(mapping[int(line.tokens[0])]), int(mapping[int(line.tokens[1])]), int(mapping[int(line.tokens[2])]), int(mapping[int(line.tokens[3])]), line.tokens[4], line.tokens[5], line.tokens[6], line.tokens[7]
                if tuple((idx1, idx2, idx3, idx4, ftype, phi0, fc, multiplicity)) == dihedral or tuple((idx4, idx3, idx2, idx1, ftype, phi0, fc, multiplicity)) == dihedral:
                    mol_mutant_dihedrals_to_add.append(Line(f"\t{idx1}  {idx2}  {idx3}  {idx4}  {ftype}  {phi0}  {fc}  {multiplicity}"))
                    break
    return mol_mutant_dihedrals_to_add


def extract_dihedrals_from_topology(molecule, mapping : Optional[dict] = None):
    dihedrals_section = molecule.get_section("dihedrals")
    dihedrals_set = set()
    for line in dihedrals_section.lines:
        if line.tokens:
            idx1, idx2, idx3, idx4, ftype, phi0, fc, multiplicity = int(line.tokens[0]), int(line.tokens[1]), int(line.tokens[2]), int(line.tokens[3]), line.tokens[4], line.tokens[5], line.tokens[6], line.tokens[7]
            if mapping:
                idx1, idx2, idx3, idx4 = mapping[idx1], mapping[idx2], mapping[idx3], mapping[idx4]
            if idx1 <= idx4:  # To avoid duplicates like (1,2,3,4) and (4,3,2,1)
                dihedrals_set.add(tuple((idx1, idx2, idx3, idx4, ftype, phi0, fc, multiplicity)))
            else:
                dihedrals_set.add(tuple((idx4, idx3, idx2, idx1, ftype, phi0, fc, multiplicity)))
    return dihedrals_set

    
    
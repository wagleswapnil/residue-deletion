from resdel.topology import *
from resdel.transformations.pair_nb import Pair_nb_object

def add_exclusions_section_to_topology(exclusions_list):
    exclusions_dict = {}
    for x in exclusions_list:
        exclusions_dict[x[0]] = exclusions_dict[x[0]] + " " + str(x[1]) if x[0] in exclusions_dict else str(x[1])
        exclusions_dict[x[1]] = exclusions_dict[x[1]] + " " + str(x[0]) if x[1] in exclusions_dict else str(x[0])
    
    exclusions_section = Section("exclusions")
    exclusions_section.add_line(Line("#ifdef EXCLS_ON"))
    for key, value in sorted(exclusions_dict.items()):
        exclusions_section.add_line(Line(f"{key}  {value}"))
    exclusions_section.add_line(Line("#endif"))
    exclusions_section.add_line(Line(""))
    exclusions_section.add_line(Line(""))
    return exclusions_section           


def get_pairs_nb_objects(pairs_list, exclusions_list, edge1_steps, sigma_epsilon_charges, comb_rule, fudge_QQ):
    pairs_nb_list = []
    for x in sorted(pairs_list):
        pairs_nb_list.append(Pair_nb_object(x, sigma_epsilon_charges[str(x[0])], sigma_epsilon_charges[str(x[1])], stateA="full_interactions", stateB="scaled_interactions", comb_rule=comb_rule, fudge_QQ=fudge_QQ))
    for x in sorted(exclusions_list):
        pairs_nb_list.append(Pair_nb_object(x, sigma_epsilon_charges[str(x[0])], sigma_epsilon_charges[str(x[1])], stateA="full_interactions", stateB="no_interactions", comb_rule=comb_rule, fudge_QQ=fudge_QQ))
    return pairs_nb_list
       

def add_pairs_nb_section_to_topology(pairs_list, exclusions_list, edge1_steps, sigma_epsilon_charges, comb_rule, fudge_QQ):
    pairs_nb_list = get_pairs_nb_objects(pairs_list, exclusions_list, edge1_steps, sigma_epsilon_charges, comb_rule, fudge_QQ)
    pairs_nb_section = Section("pairs_nb")
    for step in range(0, edge1_steps + 1):
        pairs_nb_section.add_line(Line(f"#ifdef EDGE1_STEP{step}"))
        for pair_nb in pairs_nb_list:
            pairs_nb_section.add_line(pair_nb.get_pair_nb_line(step / edge1_steps))
        pairs_nb_section.add_line(Line("#endif"))
    pairs_nb_section.add_line(Line(f"#ifdef STAGE2"))
    for pair_nb in pairs_nb_list:
        pairs_nb_section.add_line(pair_nb.get_pair_nb_line(step / edge1_steps))
    pairs_nb_section.add_line(Line("#endif"))
    pairs_nb_section.add_line(Line(""))
    pairs_nb_section.add_line(Line(""))
    return pairs_nb_section
    

def get_topB_pairs_parameters(topB, pairsB_minus_A, mapping):
    topB_pairs_to_add = []
    for pair in pairsB_minus_A:
        idx1, idx2 = pair
        topB_pairs_to_add.append(Line(f"{idx1} {idx2}  1")) 
    return topB_pairs_to_add

def updated_pairs_section(pairs_section, idx_i, topB_pairs_to_add):
    new_pairs_section = Section("pairs")
    residue_i_internal_pairs = []
    pairs_involving_residue_i = []
    for line in pairs_section.lines:
        if line.tokens:
            idx1, idx2 = int(line.tokens[0]), int(line.tokens[1])
            if idx1 in idx_i and idx2 in idx_i:
                residue_i_internal_pairs.append(line)
            elif idx1 in idx_i or idx2 in idx_i:
                pairs_involving_residue_i.append(line)
            else:
                new_pairs_section.add_line(Line(f"\t{line.tokens[0]}  {line.tokens[1]} 1 ; {' '.join(line.tokens[2:])}"))
        elif line.raw.strip() == "":
            pass
        else:
            new_pairs_section.add_line(line)

    for line in residue_i_internal_pairs:
        new_pairs_section.add_line(Line(f"\t{line.tokens[0]}  {line.tokens[1]}  {line.tokens[2]}; {' '.join(line.tokens[3:])}  internal pair for residue i"))
    for line in pairs_involving_residue_i:
        new_pairs_section.add_line(Line(f"\t{line.tokens[0]}  {line.tokens[1]}  {line.tokens[2]}; {' '.join(line.tokens[3:])}  pair involving residue i"))
    
    new_pairs_section.add_line(Line(f"#ifdef EDGE_3"))
    for line in topB_pairs_to_add:
        idx1, idx2, ftype = line.tokens[0], line.tokens[1], line.tokens[2]
        new_pairs_section.add_line(Line(f"\t{idx1}  {idx2}  {ftype} ;  Top B pair"))
    new_pairs_section.add_line(Line("#endif"))
    new_pairs_section.add_line(Line(""))
    return new_pairs_section

def get_dual_state_atoms(molA, idx_i):
    dual_state_atoms = []
    atoms_section = molA.get_section("atoms")
    for line in atoms_section.lines:
        if line.tokens:
            nr, type, resnr, residue, atom, cgnr, charge, mass = line.tokens[0:8]
            if int(nr) in idx_i:
                dual_state_atoms.append(dual_state_atom(nr, type, resnr, residue, atom, cgnr, charge, mass))
    return dual_state_atoms

def dual_state_atom(nr, type, resnr, residue, atom, cgnr, charge, mass):
    if atom.startswith("H") or atom.startswith("O") or atom.startswith("N") or atom.startswith("S") or atom.startswith("C"):
        return Line(f"\t{nr}\t{type}\t {resnr}\t {residue}\t {atom} \t{cgnr} \t{charge} \t{mass}  \tdum_{atom[0]}  \t0.0 \t {mass}; dual state atom")
    else:
        return Line(f"\t{nr}\t{type}\t {resnr}\t {residue}\t {atom} \t{cgnr} \t{charge} \t{mass}  \tdum_X  \t0.0 \t {mass}; dual state atom")

def updated_atoms_section(atoms_section, idx_i, dual_state_atoms_to_add):
    new_atoms_section = Section("atoms")
    for line in atoms_section.lines:
        if line.tokens:
            nr, type, resnr, residue, atom, cgnr, charge, mass = line.tokens[0:8]
            if int(nr) == min(idx_i):
                new_atoms_section.add_line(Line(f"#ifdef EDGE_3"))
                for dual_state_atom in dual_state_atoms_to_add:
                    new_atoms_section.add_line(dual_state_atom)
                new_atoms_section.add_line(Line(f"#else"))
                new_atoms_section.add_line(line)
            elif int(nr) == max(idx_i):
                new_atoms_section.add_line(line)
                new_atoms_section.add_line(Line(f"#endif"))
            else:
                new_atoms_section.add_line(line)
        else:
            new_atoms_section.add_line(line)
    return new_atoms_section

def add_dummy_atomtypes_to_topology(atomtypes_section):
    new_atomtypes_section = Section("atomtypes")
    for line in atomtypes_section.lines:
        if line.raw.strip() == "":
            continue
        else:
            new_atomtypes_section.add_line(line)
    new_atomtypes_section.add_line(Line(f"dum_H     0.000000     0.000000   A     0.000000     0.000000"))
    new_atomtypes_section.add_line(Line(f"dum_O     0.000000     0.000000   A     0.000000     0.000000"))
    new_atomtypes_section.add_line(Line(f"dum_N     0.000000     0.000000   A     0.000000     0.000000"))
    new_atomtypes_section.add_line(Line(f"dum_C     0.000000     0.000000   A     0.000000     0.000000"))
    new_atomtypes_section.add_line(Line(f"dum_X     0.000000     0.000000   A     0.000000     0.000000"))
    new_atomtypes_section.add_line(Line(""))
    return new_atomtypes_section

def edit_defaults_section(defaults_section):
    new_defaults_section = Section("defaults")
    for line in defaults_section.lines:
        if line.tokens:
            nbfunc, comb_rule, gen_pairs, fudgeLJ, fudgeQQ = line.tokens[0:5]
            new_fudgeLJ = "0.5"
            new_gen_pairs = "yes"
            new_defaults_section.add_line(Line(f"{nbfunc} \t\t {comb_rule} \t\t {new_gen_pairs} \t\t {new_fudgeLJ}  \t\t{fudgeQQ}"))
        else:
            new_defaults_section.add_line(line)
    return new_defaults_section
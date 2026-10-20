def atom_key(tokens, top_name, residue_to_delete):
        resnr = tokens[2]
        residue = tokens[3]
        atom = tokens[4]
        if top_name == "B" and int(resnr) >= residue_to_delete:
            resnr = str(int(resnr) + 1)
        return (resnr, residue, atom)


def build_atom_mapping(linesA, linesB, residue_to_delete):
    atomsA_dict = {}
    atomsB_dict = {}
    mapping = {}

    for i, line in enumerate(linesA):
        if line.tokens:
            atomsA_dict[atom_key(line.tokens, "A", residue_to_delete)] = line.tokens[0]

    for i, line in enumerate(linesB):
        if line.tokens:
            atomsB_dict[atom_key(line.tokens, "B", residue_to_delete)] = line.tokens[0]
    
    for key in atomsB_dict:
        if key in atomsA_dict:
            #breakpoint()
            mapping[int(atomsB_dict[key])] = int(atomsA_dict[key])
        else:
            raise ValueError(f"Warning: Atom {key} in topology B not found in topology A. This is not supposed to happen.")
    return mapping
from .topology import Topology
from .section import Section
from .molecule import Molecule
from .line import Line

def parse_topology(file_path: str) -> Topology:
    top = Topology()
    current_section = None
    current_molecule = None

    state = "HEADER"  # Possible states: HEADER, MOLECULE, TAIL

    with open(file_path, 'r') as f:
        for raw_line in f:
            line = Line(raw_line)

            if line.is_section:
                name = line.section_name
                current_section = Section(name)
                current_section.add_line(line)

                if name == "moleculetype":
                    state = "MOLECULE"
                    current_molecule = None

                elif name in ["system", "molecules", "intermolecular_interactions"]:
                    state = "TAIL"

                if state == "HEADER":
                    top.header.append(current_section)
                elif state == "MOLECULE":

                    if name == "moleculetype":
                        pass

                    if current_molecule:
                        current_molecule.add_section(current_section)
                
                elif state == "TAIL":
                    top.tail.append(current_section)

                continue

            if current_section is None:
                top.preamble.append(line)
                continue

            current_section.add_line(line)

            if current_section.name == "moleculetype" and current_molecule is None:
                if line.tokens:
                    mol_name = line.tokens[0]
                    current_molecule = Molecule(mol_name)
                    top.molecules[mol_name] =current_molecule
                    
                    current_molecule.add_section(current_section)

    return top

def parse_line(raw_line):

    line = raw_line.rstrip("\n")
    stripped = line.strip()

    # empty line
    if not stripped:
        return Line(raw=line)

    # full comment
    if stripped.startswith(";"):
        return Line(raw=line)

    # section header
    if stripped.startswith("[") and stripped.endswith("]"):

        name = stripped.strip("[] ").lower()

        return Line(
            raw=line,
            is_section=True,
            section_name=name
        )

    # inline comment
    if ";" in line:
        data, comment = line.split(";", 1)
        comment = comment.strip()
    else:
        data = line
        comment = None

    tokens = data.split()

    return Line(
        raw=line,
        tokens=tokens if tokens else None,
        comment=comment
    )
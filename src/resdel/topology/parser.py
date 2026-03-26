from .topology import Topology
from .section import Section
from .molecule import Molecule
from .line import Line

def parse_topology(file_path: str) -> Topology:
    topology = Topology()
    current_section = None

    with open(file_path, 'r') as f:
        for line in f:
            stripped_line = line.strip()
            if not stripped_line:
                continue  # Skip empty lines

            if stripped_line.startswith('[') and stripped_line.endswith(']'):
                section_name = stripped_line[1:-1].strip().lower()
                current_section = Section(section_name)
                topology.header.append(current_section)

            elif current_section is not None:
                current_section.add_line(stripped_line)
            else:
                topology.preamble.append(stripped_line)

    return topology
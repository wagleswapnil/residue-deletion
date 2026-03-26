from .topology import Topology

def write_topology(topology: Topology, file_path: str):
    with open(file_path, 'w') as f:
        # Write preamble
        for line in topology.preamble:
            f.write(line + '\n')

        # Write header sections
        for section in topology.header:
            f.write(f'[{section.name}]\n')
            for line in section.lines:
                f.write(line.raw + '\n')

        # Write molecules
        for molecule in topology.molecules:
            f.write(f'[ moleculetype ]\n{molecule.name}\n')
            for section in molecule.sections:
                f.write(f'[{section.name}]\n')
                for line in section.lines:
                    f.write(line.raw + '\n')

        # Write tail sections
        for section in topology.tail:
            f.write(f'[{section.name}]\n')
            for line in section.lines:
                f.write(line.raw + '\n')
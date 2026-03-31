from .topology import Topology
from .line import Line
from .section import Section
from .molecule import Molecule
from .topology import Topology

class Writer:
    def __init__(self, topology: Topology, file_path: str):
        self.topology = topology
        self.file_path = file_path

        
    def write_topology(self):
        with open(self.file_path, 'w') as f:
            # Write preamble
            for line in self.topology.preamble:
                f.write(line.raw + '\n')

            # Write header sections
            for section in self.topology.header:
                f.write(f'[ {section.name} ]\n')
                for line in section.lines:
                    f.write(line.raw + '\n')

            # Write molecules
            for molecule in self.topology.molecules:
                #f.write(f'[ moleculetype ]\n{molecule.name}\n')
                for section in molecule.sections:
                    f.write(f'[ {section.name} ]\n')
                    for line in section.lines:
                        f.write(line.raw + '\n')

            # Write tail sections
            for section in self.topology.tail:
                f.write(f'[ {section.name} ]\n')
                for line in section.lines:
                    f.write(line.raw + '\n')
        return
from .topology import Topology
from .line import Line
from .section import Section
from .molecule import Molecule
from .topology import Topology
from .formatter import GromacsFormatter, Formatter

class Writer:
    def __init__(self, topology: Topology, file_path: str, formatter=None):
        self.topology = topology
        self.file_path = file_path
        self.formatter = formatter if formatter is not None else Formatter()

        
    def write_topology(self):
        with open(self.file_path, 'w') as f:
            # Write preamble
            for line in self.topology.preamble:
                formatted_line = self.formatter.format_line(line)
                f.write(formatted_line + '\n')

            # Write header sections
            for section in self.topology.header:
                f.write(f'[ {section.name} ]\n')
                for line in section.lines:
                    formatted_line = self.formatter.format_line(line, section.name)
                    f.write(formatted_line + '\n')

            # Write molecules
            for molecule in self.topology.molecules:
                #f.write(f'[ moleculetype ]\n{molecule.name}\n')
                for section in molecule.sections:
                    f.write(f'[ {section.name} ]\n')
                    for line in section.lines:
                        formatted_line = self.formatter.format_line(line, section.name)
                        f.write(formatted_line + '\n')

            # Write tail sections
            for section in self.topology.tail:
                f.write(f'[ {section.name} ]\n')
                for line in section.lines:
                    formatted_line = self.formatter.format_line(line, section.name)
                    f.write(formatted_line + '\n')
        return
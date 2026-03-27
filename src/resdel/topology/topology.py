from typing import List
from .section import Section
from .molecule import Molecule

class Topology():
    def __init__(self):
        self.preamble: List[str] = []
        self.header: List[Section] = []
        self.molecules: List[Molecule] = []
        self.tail: List[Section] = []

    def add_preamble_line(self, line: str):
        self.preamble.append(line)

    def add_header_section(self, section: Section):
        self.header.append(section)

    def add_molecule(self, molecule: Molecule):
        self.molecules.append(molecule)

    def add_tail_section(self, section: Section):
        self.tail.append(section)

    def get_section_by_name(self, name: str) -> Section:
        name = name.lower()

        for s in self.header_sections + self.sections:
            if s.name == name:
                return s
    
    def __repr__(self):
        return (f"Topology(preamble={self.preamble}, header={self.header}, "
                f"molecules={self.molecules}, tail={self.tail})")

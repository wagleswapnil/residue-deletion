from typing import List
from .line import Line
from .section import Section
from .molecule import Molecule

class Topology():
    def __init__(self):
        self.preamble: List[str] = []
        self.header: List[Section] = []
        self.molecules: List[Molecule] = []
        self.tail: List[Section] = []

    def add_preamble_line(self, line: Line):
        self.preamble.append(line)

    def add_header_section(self, section: Section):
        self.header.append(section)

    def add_molecule(self, molecule: Molecule):
        self.molecules.append(molecule)

    def add_tail_section(self, section: Section):
        self.tail.append(section)

    def get_header_section_by_name(self, name: str) -> Section:
        name = name.lower()
        for s in self.header:
            if s.name == name:
                return s
            
    def replace_header_section_by_name(self, name: str, new_section: Section):
        for s, section in enumerate(self.header):
            if section.name == name:
                self.header[s] = new_section
                return
        raise ValueError(f"Section with name {name} not found in topology header")

    def remove_molecule_by_name(self, name: str):
        for i, molecule in enumerate(self.molecules):
            if molecule.name == name:
                del self.molecules[i]
                return
        
    def get_tail_section_by_name(self, name: str):
        name = name.lower()
        for s in self.tail:
            if s.name == name:
                return s

    def replace_tail_section_by_name(self, name: str, new_section: Section):
        for i, section in enumerate(self.tail):
            if section.name == name:
                self.tail[i] = new_section
                return
        raise ValueError(f"Section with name {name} not found in topology tail")

    def get_comb_rule_fudgeQQ(self):
        comb_rule = None
        fudgeQQ = None
        for line in self.get_header_section_by_name("defaults").lines:
            if line.tokens:
                comb_rule = line.tokens[1]
                fudgeQQ = line.tokens[4]
        return comb_rule, fudgeQQ

    
    def __repr__(self):
        return (f"Topology(preamble={self.preamble}, header={self.header}, "
                f"molecules={self.molecules}, tail={self.tail})")

from typing import List
from .section import Section

class Molecule:
    def __init__(self, name: str):
        self.name = name
        self.sections: List[Section] = []

    def add_section(self, section: Section):
        self.sections.append(section)

    def __repr__(self):
        return f"Molecule(name={self.name}, sections={self.sections})"
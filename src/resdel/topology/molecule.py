from typing import List
from .section import Section

class Molecule:
    def __init__(self, name: str):
        self.name = name
        self.sections: List[Section] = []

    def add_section(self, section: Section):
        self.sections.append(section)

    def get_section(self, name):
        name = name.lower()
        for s in self.sections:
            if s.name == name:
                return s
        return None

    def replace_section(self, name, new_section):
        for s, section in enumerate(self.sections):
            if section.name == name:
                self.sections[s] = new_section
                return
        raise ValueError(f"Section with name {name} not found in molecule {self.name}")
        return



    def __repr__(self):
        return f"Molecule(name={self.name}, sections={self.sections})"
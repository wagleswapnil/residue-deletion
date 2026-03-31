from .line import Line
from .section import Section
from .molecule import Molecule
from .topology import Topology


class Parser:
    def __init__(self, file_path: str):
        self.file_path = file_path

    def parse_topology(self) -> Topology:
        self.top = Topology()
        self.current_section: Section = None
        self.current_molecule: Molecule = None
        self.current_state = None  # Possible states: HEADER, MOLECULE, TAIL
        count = 0
        with open(self.file_path, 'r') as f:
            for raw_line in f:
                line = Line(raw_line)
                if self.current_state is None:
                    if line.is_section:
                        self.current_section = Section(line.section_name)
                        self.current_state = "HEADER"
                    else:
                        self.top.add_preamble_line(line)
                        #continue
                elif self.current_state == "HEADER":
                    if line.is_section:
                        if line.section_name == "moleculetype":
                            self.store_section(self.current_section)
                            self.current_section = Section(line.section_name)
                            self.current_state = "MOLECULE"
                        else:
                            self.store_section(self.current_section)
                            self.current_section = Section(line.section_name)
                    else:
                        self.current_section.add_line(line)
                elif self.current_state == "MOLECULE":
                    if line.is_section:
                        if line.section_name in ["system", "molecules", "intermolecular_interactions"]:
                            self.store_section(self.current_section)
                            self.current_molecule.add_section(self.current_section)
                            self.top.add_molecule(self.current_molecule)
                            self.current_section = Section(line.section_name)
                            self.current_molecule = None
                            self.current_state = "TAIL"
                        elif line.section_name == "moleculetype":
                            self.store_section(self.current_section)
                            self.current_molecule.add_section(self.current_section)
                            self.top.add_molecule(self.current_molecule)
                            self.current_section = Section(line.section_name)
                            self.current_molecule = None
                        else:
                            self.store_section(self.current_section)
                            self.current_molecule.add_section(self.current_section)
                            self.current_section = Section(line.section_name)
                    else:
                        self.current_section.add_line(line)
                elif self.current_state == "TAIL":
                    if line.is_section:
                        self.store_section(self.current_section)
                        self.current_section = Section(line.section_name)
                    else:
                        self.current_section.add_line(line)
            self.store_section(self.current_section) # storing the final section in TAIL state
                #count += 1
                #if count > 340: 
                #    breakpoint()
        return self.top
                    

    def store_section(self, section: Section):
        if self.current_state == "HEADER":
            self.top.add_header_section(section)
        elif self.current_state == "MOLECULE":
            if section.name == "moleculetype":
                for line in section.lines:
                    if line.raw and line.raw[0] != ";":  # Skip comment and empty lines
                        mol_name = line.tokens[0]
                        break
                self.current_molecule = Molecule(mol_name)
        elif self.current_state == "TAIL":
            self.top.add_tail_section(section)
        return
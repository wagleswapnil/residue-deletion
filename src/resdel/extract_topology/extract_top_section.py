from resdel.topology import Parser, Line, Writer
from resdel.topology.formatter import GromacsFormatter
from resdel.topology.section import Section

class ExtractTopSection:
    def __init__(self, topology_path: str, molecule_name: str, begin_idx: str, end_idx: str, first_residue_idx: str):
        self.topology_path = topology_path
        self.molecule_name = molecule_name
        self.begin_idx = int(begin_idx)
        self.end_idx = int(end_idx)
        self.first_residue_idx = int(first_residue_idx)
        self.topology = Parser(file_path=self.topology_path).parse_topology()
        self.molecule_name = molecule_name if molecule_name else "system1"
        self.mol = None
        for mol in self.topology.molecules:
            if mol.name == self.molecule_name:
                self.mol = mol
        if self.mol is None:
            raise ValueError(f"Molecule '{self.molecule_name}' not found in the topology.")
        
    def _new_atom_idx(self, old_idx):
        return old_idx - self.begin_idx + 1

    def _new_residue_idx(self, old_resnr):
        return old_resnr - self.first_residue_idx + 1

    def _extract_atoms_section(self):
        new_atoms_section = Section("atoms")
        for line in self.mol.get_section("atoms").lines:
            if line.tokens:
                nr, type, resnr, residue, atom, cgnr, charge, mass = line.tokens[0:8]
                nr = int(nr)
                cgnr = int(cgnr)
                if self.begin_idx <= nr <= self.end_idx:
                    new_atoms_section.add_line(Line(f"\t{self._new_atom_idx(nr)}\t{type}\t {self._new_residue_idx(int(resnr))}\t {residue}\t {atom} \t{self._new_atom_idx(cgnr)} \t{charge} \t{mass}"))
            else:
                new_atoms_section.add_line(line)
        return new_atoms_section

    def _extract_bonds_section(self):
        new_bonds_section = Section("bonds")
        for line in self.mol.get_section("bonds").lines:
            if line.tokens:
                if line.tokens[0].startswith("#"):
                    new_bonds_section.add_line(line)
                else:
                    nr1, nr2, funct, d0, fc = line.tokens[0:5]
                    nr1, nr2 = int(nr1), int(nr2)
                    if self.begin_idx <= nr1 <= self.end_idx and self.begin_idx <= nr2 <= self.end_idx:
                        new_bonds_section.add_line(Line(f"\t{self._new_atom_idx(nr1)}\t{self._new_atom_idx(nr2)}\t{funct}\t{d0}\t{fc}"))
            else:
                new_bonds_section.add_line(line)
        return new_bonds_section
        
    def _extract_pairs_section(self):
        new_pairs_section = Section("pairs")
        for line in self.mol.get_section("pairs").lines:
            if line.tokens:
                if line.tokens[0].startswith("#"):
                    new_pairs_section.add_line(line)
                else:
                    nr1, nr2, funct, c0, c1 = line.tokens[0:5]
                    nr1, nr2 = int(nr1), int(nr2)
                    if self.begin_idx <= nr1 <= self.end_idx and self.begin_idx <= nr2 <= self.end_idx:
                        new_pairs_section.add_line(Line(f"\t{self._new_atom_idx(nr1)}\t{self._new_atom_idx(nr2)}\t{funct}\t{c0}\t{c1}"))
            else:
                new_pairs_section.add_line(line)
        return new_pairs_section    

    
    def _extract_angles_section(self):
        new_angles_section = Section("angles")
        for line in self.mol.get_section("angles").lines:
            if line.tokens:
                if line.tokens[0].startswith("#"):
                    new_angles_section.add_line(line)
                else:
                    nr1, nr2, nr3, funct, theta0, fc = line.tokens[0:6]
                    nr1, nr2, nr3 = int(nr1), int(nr2), int(nr3)
                    if self.begin_idx <= nr1 <= self.end_idx and self.begin_idx <= nr2 <= self.end_idx and self.begin_idx <= nr3 <= self.end_idx:
                        new_angles_section.add_line(Line(f"\t{self._new_atom_idx(nr1)}\t{self._new_atom_idx(nr2)}\t{self._new_atom_idx(nr3)}\t{funct}\t{theta0}\t{fc}"))
            else:
                new_angles_section.add_line(line)
        return new_angles_section

    def _extract_dihedrals_section(self):
        new_dihedrals_section = Section("dihedrals")
        for line in self.mol.get_section("dihedrals").lines:
            if line.tokens:
                if line.tokens[0].startswith("#"):
                    new_dihedrals_section.add_line(line)
                else:
                    nr1, nr2, nr3, nr4, funct, phi0, fc, multiplicity = line.tokens[0:8]
                    nr1, nr2, nr3, nr4 = int(nr1), int(nr2), int(nr3), int(nr4)
                    if self.begin_idx <= nr1 <= self.end_idx and self.begin_idx <= nr2 <= self.end_idx and self.begin_idx <= nr3 <= self.end_idx and self.begin_idx <= nr4 <= self.end_idx:
                        new_dihedrals_section.add_line(Line(f"\t{self._new_atom_idx(nr1)}\t{self._new_atom_idx(nr2)}\t{self._new_atom_idx(nr3)}\t{self._new_atom_idx(nr4)}\t{funct}\t{phi0}\t{fc}\t{multiplicity}"))
            else:
                new_dihedrals_section.add_line(line)
        return new_dihedrals_section

    def write_residue_topology(self, output_path):
        writer = Writer(file_path=output_path, topology=self.topology, formatter=GromacsFormatter())
        writer.write_topology()
        return

    def extract_topology_sections(self):
        self.mol.replace_section("atoms", self._extract_atoms_section())
        self.mol.replace_section("bonds", self._extract_bonds_section())
        self.mol.replace_section("pairs", self._extract_pairs_section())
        self.mol.replace_section("angles", self._extract_angles_section())
        self.mol.replace_section("dihedrals", self._extract_dihedrals_section())
        return

    def edit_defaults_section(self, defaults_section):
        new_defaults_section = Section("defaults")
        for line in defaults_section.lines:
            if line.tokens:
                nbfunc, comb_rule, gen_pairs, fudgeLJ, fudgeQQ = line.tokens[0:5]
                new_fudgeLJ = "0.5"
                new_gen_pairs = "yes"
                new_defaults_section.add_line(Line(f"{nbfunc} \t\t {comb_rule} \t\t {new_gen_pairs} \t\t {new_fudgeLJ}  \t\t{fudgeQQ}"))
            else:
                new_defaults_section.add_line(line)
        return new_defaults_section

    def add_dummy_atomtypes_to_topology(self, atomtypes_section):
        new_atomtypes_section = Section("atomtypes")
        for line in atomtypes_section.lines:
            if line.raw.strip() == "":
                continue
            else:
                new_atomtypes_section.add_line(line)
        new_atomtypes_section.add_line(Line(f"dum_H     0.000000     0.000000   A     0.000000     0.000000"))
        new_atomtypes_section.add_line(Line(f"dum_O     0.000000     0.000000   A     0.000000     0.000000"))
        new_atomtypes_section.add_line(Line(f"dum_N     0.000000     0.000000   A     0.000000     0.000000"))
        new_atomtypes_section.add_line(Line(f"dum_C     0.000000     0.000000   A     0.000000     0.000000"))
        new_atomtypes_section.add_line(Line(f"dum_X     0.000000     0.000000   A     0.000000     0.000000"))
        new_atomtypes_section.add_line(Line(""))
        return new_atomtypes_section


    def _get_last_residue_atom_idxs(self):
        for line in self.mol.get_section("atoms").lines:
            if line.tokens:
                nr, type, resnr, residue, atom, cgnr, charge, mass = line.tokens[0:8]
                if int(resnr) == self.last_residue:
                    return int(nr)
        return None  # Return None if the last residue is not found


    def _make_terminal_residue_atoms_section(self):
        new_atoms_section = Section("atoms")
        for line in self.mol.get_section("atoms").lines:
            if line.tokens:
                nr, type, resnr, residue, atom, cgnr, charge, mass = line.tokens[0:8]
                nr = int(nr)
                cgnr = int(cgnr)
                if int(resnr) == self.last_residue:
                    new_atoms_section.add_line(Line(f"\t{nr}\t{type}\t {resnr}\t {residue}\t {atom} \t{cgnr} \t{charge} \t{mass}\tdum_{atom[0]}  \t0.0 \t {mass}; dual state atom"))
                else:
                    new_atoms_section.add_line(line)
            else:
                new_atoms_section.add_line(line)
        return new_atoms_section
    

    def _make_terminal_residue_pairs_section(self):
        new_pairs_section = Section("pairs")
        for line in self.mol.get_section("pairs").lines:
            if line.tokens:
                nr1, nr2, funct, c0, c1 = line.tokens[0:5]
                new_pairs_section.add_line(Line(f"{nr1}\t{nr2}\t {funct} ; \t{c0}\t{c1}"))
            else:
                new_pairs_section.add_line(line)
        return new_pairs_section

    def _make_terminal_residue_angles_section(self):
        new_angles_section = Section("angles")
        for line in self.mol.get_section("angles").lines:
            if line.tokens:
                nr1, nr2, nr3, ftype, theta0, fc = line.tokens[:6]
                nr1, nr2, nr3 = int(nr1), int(nr2), int(nr3)
                if nr1 >= self.last_res_first_atom_idx and nr2 >= self.last_res_first_atom_idx and nr3 >= self.last_res_first_atom_idx:
                    new_angles_section.add_line(line)
                elif nr1 >= self.last_res_first_atom_idx or nr2 >= self.last_res_first_atom_idx or nr3 >= self.last_res_first_atom_idx:
                    new_angles_section.add_line(Line(f"{nr1}\t{nr2}\t{nr3}\t{ftype}\t{theta0}\t{fc}\t{theta0}\t 0.0"))
                else:
                    new_angles_section.add_line(line)
            else:
                new_angles_section.add_line(line)
        return new_angles_section

    def _make_terminal_residue_dihedrals_section(self):
        new_dihedrals_section = Section("dihedrals")
        for line in self.mol.get_section("dihedrals").lines:
            if line.tokens:
                nr1, nr2, nr3, nr4, ftype, theta0, fc, multiplicity = line.tokens[:9]
                nr1, nr2, nr3, nr4 = int(nr1), int(nr2), int(nr3), int(nr4)
                if nr1 >= self.last_res_first_atom_idx and nr2 >= self.last_res_first_atom_idx and nr3 >= self.last_res_first_atom_idx and nr4 >= self.last_res_first_atom_idx:
                    new_dihedrals_section.add_line(line)
                elif nr1 >= self.last_res_first_atom_idx or nr2 >= self.last_res_first_atom_idx or nr3 >= self.last_res_first_atom_idx or nr4 >= self.last_res_first_atom_idx:
                    new_dihedrals_section.add_line(Line(f"{nr1}\t{nr2}\t{nr3}\t{nr4}\t{ftype}\t{theta0}\t{fc}\t{multiplicity}\t{theta0}\t 0.0 \t {multiplicity}"))
                else:
                    new_dihedrals_section.add_line(line)
            else:
                new_dihedrals_section.add_line(line)
        return new_dihedrals_section
    

    def make_last_terminal_residue_dual_topology(self, last_residue): # Run this method after extract_topology_sections() to create a dual topology for the terminal residue on the updated topology.
        self.last_residue = last_residue
        self.topology.replace_header_section_by_name("defaults", self.edit_defaults_section(self.topology.get_header_section_by_name("defaults")))
        self.topology.replace_header_section_by_name("atomtypes", self.add_dummy_atomtypes_to_topology(self.topology.get_header_section_by_name("atomtypes")))
        self.last_res_first_atom_idx = self._get_last_residue_atom_idxs()
        self.mol.replace_section("atoms", self._make_terminal_residue_atoms_section())
        self.mol.replace_section("pairs", self._make_terminal_residue_pairs_section())
        self.mol.replace_section("angles", self._make_terminal_residue_angles_section())
        self.mol.replace_section("dihedrals", self._make_terminal_residue_dihedrals_section())
        return
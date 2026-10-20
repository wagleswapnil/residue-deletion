from resdel.topology import Parser, Line, Writer
from resdel.topology.formatter import GromacsFormatter

from resdel.topology.section import Section

class TopologyExtractor:
    def __init__(self, topology_path):
        self.topology_path = topology_path
        self.parser = Parser(topology_path)
        self.topology = self.parser.parse_topology()

    def extract_molecule_topology(self, molecule_name, residue_index):
        self.molecule_name = molecule_name if molecule_name else "system1"
        self.mol = None
        print(self.molecule_name)
        for mol in self.topology.molecules:
            if mol.name == self.molecule_name:
                self.mol = mol
        if self.mol is None:
            raise ValueError(f"Molecule '{self.molecule_name}' not found in the topology.")
        self.residue_atoms_indices = self._get_residue_atoms_indices(residue_index)
        self.idx_subtract_by = min(self.residue_atoms_indices) - 1
        self.resnr_subtract_by = residue_index - 1
        return

    def build_residue_topology(self):
        self.mol.replace_section("atoms", self._build_residue_atoms_section())
        self.mol.replace_section("bonds", self._build_residue_bonds_section())
        self.mol.replace_section("pairs", self._build_residue_pairs_section())
        self.mol.replace_section("angles", self._build_residue_angles_section())
        self.mol.replace_section("dihedrals", self._build_residue_dihedrals_section())
        return

    def write_residue_topology(self, output_path):
        writer = Writer(file_path=output_path, topology=self.topology, formatter=GromacsFormatter())
        writer.write_topology()
        return

    def _build_residue_atoms_section(self):
        new_atoms_section = Section("atoms")
        for line in self.mol.get_section("atoms").lines:
            if line.tokens:
                nr, type, resnr, residue, atom, cgnr, charge, mass = line.tokens[0:8]
                nr = int(nr)
                if nr in self.residue_atoms_indices:
                    new_atoms_section.add_line(Line(f"\t{self._new_atom_idx(nr)}\t{type}\t {(self._new_residue_idx(int(resnr)))}\t {residue}\t {atom} \t{self._new_atom_idx(int(cgnr))} \t{charge} \t{mass}"))
            else:
                new_atoms_section.add_line(line)
        return new_atoms_section

    def _build_residue_bonds_section(self):
        new_bonds_section = Section("bonds")
        for line in self.mol.get_section("bonds").lines:
            if line.tokens:
                if line.tokens[0].startswith("#"):
                    new_bonds_section.add_line(line)
                else:
                    nr1, nr2, funct, d0, fc = line.tokens[0:5]
                    nr1, nr2 = int(nr1), int(nr2)
                    if nr1 in self.residue_atoms_indices and nr2 in self.residue_atoms_indices:
                        new_bonds_section.add_line(Line(f"\t{self._new_atom_idx(nr1)}\t{self._new_atom_idx(nr2)}\t{funct}\t{d0}\t{fc}"))
            else:
                new_bonds_section.add_line(line)
        return new_bonds_section

    def _build_residue_pairs_section(self):
        new_pairs_section = Section("pairs")
        for line in self.mol.get_section("pairs").lines:
            if line.tokens:
                if line.tokens[0].startswith("#"):
                    new_pairs_section.add_line(line)
                else:
                    nr1, nr2, funct, d0, fc = line.tokens[0:5]
                    nr1, nr2 = int(nr1), int(nr2)
                    if nr1 in self.residue_atoms_indices and nr2 in self.residue_atoms_indices:
                        new_pairs_section.add_line(Line(f"\t{self._new_atom_idx(nr1)}\t{self._new_atom_idx(nr2)}\t{funct}\t{d0}\t{fc}"))
            else:
                new_pairs_section.add_line(line)
        return new_pairs_section

    def _build_residue_angles_section(self):
        new_angles_section = Section("angles")
        for line in self.mol.get_section("angles").lines:
            if line.tokens:
                if line.tokens[0].startswith("#"):
                    new_angles_section.add_line(line)
                else:
                    nr1, nr2, nr3, funct, theta0, fc = line.tokens[0:6]
                    nr1, nr2, nr3 = int(nr1), int(nr2), int(nr3)
                    if nr1 in self.residue_atoms_indices and nr2 in self.residue_atoms_indices and nr3 in self.residue_atoms_indices:
                        new_angles_section.add_line(Line(f"\t{self._new_atom_idx(nr1)}\t{self._new_atom_idx(nr2)}\t{self._new_atom_idx(nr3)}\t{funct}\t{theta0}\t{fc}"))
            else:
                new_angles_section.add_line(line)
        return new_angles_section

    def _build_residue_dihedrals_section(self):
        new_dihedrals_section = Section("dihedrals")
        for line in self.mol.get_section("dihedrals").lines:
            if line.tokens:
                if line.tokens[0].startswith("#"):
                    new_dihedrals_section.add_line(line)
                else:
                    nr1, nr2, nr3, nr4, funct, phi0, fc, multiplicity = line.tokens[0:8]
                    nr1, nr2, nr3, nr4 = int(nr1), int(nr2), int(nr3), int(nr4)
                    if nr1 in self.residue_atoms_indices and nr2 in self.residue_atoms_indices and nr3 in self.residue_atoms_indices and nr4 in self.residue_atoms_indices:
                        new_dihedrals_section.add_line(Line(f"\t{self._new_atom_idx(nr1)}\t{self._new_atom_idx(nr2)}\t{self._new_atom_idx(nr3)}\t{self._new_atom_idx(nr4)}\t{funct}\t{phi0}\t{fc}\t{multiplicity}"))
            else:
                new_dihedrals_section.add_line(line)
        return new_dihedrals_section


    def _get_residue_atoms_indices(self, residue_index):
        residue_atoms_indices = []
        for line in self.mol.get_section("atoms").lines:
            if line.tokens:
                nr, resnr = line.tokens[0], line.tokens[2]
                if int(resnr) == residue_index:
                    residue_atoms_indices.append(int(nr))
        return residue_atoms_indices

    def _new_atom_idx(self, idx):
        return idx - self.idx_subtract_by

    def _new_residue_idx(self, idx):
        return idx - self.resnr_subtract_by

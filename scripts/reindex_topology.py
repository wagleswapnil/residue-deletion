import sys
import resdel
from resdel.topology import Parser, Section, Line, Writer
from resdel.topology.formatter import GromacsFormatter

top_path = sys.argv[1]
mapping_path = sys.argv[2]
out_path = sys.argv[3]

parser = Parser(top_path)
topology = parser.parse_topology()

mapping = {}
with open (mapping_path, 'r') as f:
    for line in f.readlines():
        mapping[int(line.split()[0])] = int(line.split()[1].strip())
f.close()

def change_atoms(mol):
    new_atoms_section = Section("atoms")
    for line in mol.get_section("atoms").lines:
            if line.tokens:
                if line.raw.strip()[0] == "#":
                    new_atoms_section.add_line(line)
                else:    
                    nr, type, resnr, residue, atom, cgnr = line.tokens[0:6]
                    nr = int(nr)
                    cgnr = int(cgnr)
                    new_atoms_section.add_line(Line(f"\t{mapping[nr]}\t{type}\t {resnr}\t {residue}\t {atom} \t{mapping[cgnr]} \t{"\t".join(line.tokens[6:])}"))
            else:
                new_atoms_section.add_line(line)
    return new_atoms_section
    
def change_bonds(mol):
    new_bonds_section = Section("bonds")
    for line in mol.get_section("bonds").lines:
        if line.tokens:
            if line.tokens[0].startswith("#"):
                new_bonds_section.add_line(line)
            else:
                nr1, nr2 = line.tokens[0:2]
                nr1, nr2 = int(nr1), int(nr2)
                new_bonds_section.add_line(Line(f"\t{mapping[nr1]}\t{mapping[nr2]}\t{"\t".join(line.tokens[2:])}"))
        else:
            new_bonds_section.add_line(line)
    return new_bonds_section

def change_pairs(mol):
    new_bonds_section = Section("pairs")
    for line in mol.get_section("pairs").lines:
        if line.tokens:
            if line.tokens[0].startswith("#"):
                new_bonds_section.add_line(line)
            else:
                nr1, nr2 = line.tokens[0:2]
                nr1, nr2 = int(nr1), int(nr2)
                new_bonds_section.add_line(Line(f"\t{mapping[nr1]}\t{mapping[nr2]}\t{"\t".join(line.tokens[2:])}"))
        else:
            new_bonds_section.add_line(line)
    return new_bonds_section

def change_angles(mol):
    new_angles_section = Section("angles")
    for line in mol.get_section("angles").lines:
        if line.tokens:
            if line.tokens[0].startswith("#"):
                new_angles_section.add_line(line)
            else:
                nr1, nr2, nr3 = line.tokens[0:3]
                nr1, nr2, nr3 = int(nr1), int(nr2), int(nr3)
                new_angles_section.add_line(Line(f"\t{mapping[nr1]}\t{mapping[nr2]}\t{mapping[nr3]}\t{"\t".join(line.tokens[3:])}"))
        else:
            new_angles_section.add_line(line)
    return new_angles_section

def change_dihedrals(mol):
    new_dihedrals_section = Section("dihedrals")
    for line in mol.get_section("dihedrals").lines:
        if line.tokens:
            if line.tokens[0].startswith("#"):
                new_dihedrals_section.add_line(line)
            else:
                nr1, nr2, nr3, nr4 = line.tokens[0:4]
                nr1, nr2, nr3, nr4 = int(nr1), int(nr2), int(nr3), int(nr4)
                new_dihedrals_section.add_line(Line(f"\t{mapping[nr1]}\t{mapping[nr2]}\t{mapping[nr3]}\t{mapping[nr4]}\t{"\t".join(line.tokens[4:])}"))
        else:
            new_dihedrals_section.add_line(line)
    return new_dihedrals_section

def write_residue_topology(topology):
        writer = Writer(file_path=out_path, topology=topology, formatter=GromacsFormatter())
        writer.write_topology()
        return

def engine():
    mol = None
    for molecule in topology.molecules:
        if molecule.name == "system1":
            mol = molecule

    mol.replace_section("atoms", change_atoms(mol))
    mol.replace_section("bonds", change_bonds(mol))
    mol.replace_section("pairs", change_pairs(mol))
    mol.replace_section("angles", change_angles(mol))
    mol.replace_section("dihedrals", change_dihedrals(mol))
    write_residue_topology(topology)
    return

engine()
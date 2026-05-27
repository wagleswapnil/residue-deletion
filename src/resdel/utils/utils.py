import subprocess
from resdel.topology.parser import Parser
from resdel.topology.writer import Writer
from resdel.topology.section import Section
from resdel.topology.formatter import GromacsFormatter

def generate_vaccuum_structure_from_solvent_structure(in_PDB_path, out_PDB_path):
    resnames = ["HOH", "SOL", "WAT", "NA", "CL", "Na", "Cl", "K"]
    cmd = f"pdb_delresname -{','.join(resnames)} {in_PDB_path} > {out_PDB_path}"
    subprocess.run(cmd, shell=True, check=True)
    return

def generate_vaccuum_topology_from_solvent_topology(in_topology_path, out_topology_path):
    resnames = ["HOH", "SOL", "WAT", "NA", "CL", "Na", "Cl", "K"]
    parser = Parser(in_topology_path)
    topology = parser.parse_topology()

    topology.replace_tail_section_by_name("molecules", updated_molecules_section(topology.get_tail_section_by_name("molecules"), resnames)) 

    writer = Writer(topology, out_topology_path, GromacsFormatter())
    writer.write_topology()
    return

def updated_molecules_section(molecules_section, resnames):
        new_molecules_section = Section("molecules")
        for line in molecules_section.lines:
            if line.tokens:
                molecule = line.tokens[0]
                if molecule in resnames:
                    pass
                else:
                    new_molecules_section.add_line(line)
            else:
                new_molecules_section.add_line(line)
        return new_molecules_section

def copy_structure(in_path, out_path):
    cmd = f"cp {in_path} {out_path}"
    subprocess.run(cmd, shell=True, check=True)
    return
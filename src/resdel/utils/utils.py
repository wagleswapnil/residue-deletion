import subprocess
from resdel.topology.parser import Parser
from resdel.topology.writer import Writer
from resdel.topology.section import Section
from resdel.topology.formatter import GromacsFormatter

def generate_vacuum_structure_from_solvent_structure(in_PDB_path, out_PDB_path):
    moltypes = ["HOH", "SOL", "WAT", "NA", "CL", "Na", "Cl", "K"]
    cmd = f"pdb_delresname -{','.join(moltypes)} {in_PDB_path} > {out_PDB_path}"
    subprocess.run(cmd, shell=True, check=True)
    return

# Generates vacuum topology of a system from its solvent topology by removing moleculetypes in moltypes list.
def generate_vacuum_topology_from_solvent_topology(in_topology_path, out_topology_path):
    moltypes = ["HOH", "SOL", "WAT", "NA", "CL", "Na", "Cl", "K"]
    parser = Parser(in_topology_path)
    topology = parser.parse_topology()

    topology.replace_tail_section_by_name("molecules", updated_molecules_section(topology.get_tail_section_by_name("molecules"), moltypes)) 

    writer = Writer(topology, out_topology_path, GromacsFormatter())
    writer.write_topology()
    return

#This function removes all moleculetypes in moltypes from the [molecules] section of a topology.
def updated_molecules_section(molecules_section, moltypes):
        new_molecules_section = Section("molecules")
        for line in molecules_section.lines:
            if line.tokens:
                molecule = line.tokens[0]
                if molecule in moltypes:
                    pass
                else:
                    new_molecules_section.add_line(line)
            else:
                new_molecules_section.add_line(line)
        return new_molecules_section

# This function is used by transformations module. 
def copy_structure(in_path, out_path):
    cmd = f"cp {in_path} {out_path}"
    subprocess.run(cmd, shell=True, check=True)
    return
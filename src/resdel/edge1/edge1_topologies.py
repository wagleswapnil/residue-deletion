from resdel.topology import *
from typing import Optional
from resdel.tranformations.atom_mapping import build_atom_mapping
from resdel.edge1.edge1_utils import *
from resdel.edge1.transformations_nonbonded import *
from resdel.edge1.transformations_bonded import *
from resdel.topology.formatter import GromacsFormatter

class Edge1_Topologies:
    def __init__(self, topA : Topology, topB : Topology, residue_to_delete : str, molA_name : Optional[str] = "system1", molB_name : Optional[str] = "system1", edge1_steps : Optional[int] = 15):
        self.topA = topA
        self.topB = topB
        self.edge1_steps = edge1_steps
        self.residue_to_delete = residue_to_delete
        self.molA = None
        self.molB = None
        self.molA_name = molA_name
        self.molB_name = molB_name
        for mol in self.topA.molecules:
            if mol.name == self.molA_name:
                self.molA = mol
                break
        for mol in self.topB.molecules:
            if mol.name == self.molB_name:
                self.molB = mol
                break
        if self.molA is None:
            raise ValueError(f"Could not find molecule with name {self.molA_name} in topology A")
        if self.molB is None:
            raise ValueError(f"Could not find molecule with name {self.molB_name} in topology B")
        self.mapping = build_atom_mapping(self.molA.get_section("atoms").lines, self.molB.get_section("atoms").lines, int(self.residue_to_delete))
        self.comb_rule, self.fudge_QQ = self.topA.get_comb_rule_fudgeQQ()


    def get_residue_atom_info(self):
        self.idx_i_minus_1, self.idx_i, self.idx_i_plus_1, self.idx_i_minus_1_N, self.idx_i_minus_1_C, self.idx_i_N, self.idx_i_C, self.idx_i_plus_1_N,  self.idx_i_plus_1_C = get_residue_atom_idxs(self.molA, self.residue_to_delete)
        self.sigma_epsilon_charges = get_sigma_epsilon_charges(self.topA.get_header_section_by_name("atomtypes"), self.molA, self.idx_i_minus_1 + self.idx_i + self.idx_i_plus_1)
        return 


    def compute_pairs_exclusions(self):
        get_tpr_dump(mdp="./tests/MDP/em.mdp", structure="./tests/data/minimized_stage1.gro", topology="./tests/data/system_stage1.top", output_prefix="./tests/data/system_stage1")
        output_prefix="./tests/data/system_stage1"
        exclusionsA = extract_exclusions_from_tpr_dump(f"{output_prefix}.txt", f"{output_prefix}_exclusions.txt")
        pairsA = extract_pairs_from_topology(self.topA, self.molA_name)
        assert(pairsA.issubset(exclusionsA)), "Error: Not all pairs in topology A are present in the exclusions extracted from the tpr dump. This is not supposed to happen."
        
        add_peptide_bond(self.molA, self.idx_i_minus_1_C, self.idx_i_plus_1_N)
        topology_writer = Writer(self.topA, "./tests/data/test.top")
        topology_writer.write_topology()
        get_tpr_dump(mdp="./tests/MDP/em_test.mdp", structure="./tests/data/minimized_stage1.gro", topology="./tests/data/test.top", output_prefix="./tests/data/test")
        output_prefix="./tests/data/test"
        exclusions_temp = extract_exclusions_from_tpr_dump(f"{output_prefix}.txt", f"{output_prefix}_exclusions.txt")
        
        get_tpr_dump(mdp="./tests/MDP/em.mdp", structure="./tests/data/minimized_stage5.gro", topology="./tests/data/system_stage5.top", output_prefix="./tests/data/system_stage5")
        output_prefix="./tests/data/system_stage5"
        exclusionsB = map_exclusions_pairs(extract_exclusions_from_tpr_dump(f"{output_prefix}.txt", f"{output_prefix}_exclusions.txt"), self.mapping)
        pairsB = map_exclusions_pairs(extract_pairs_from_topology(self.topB, self.molB_name), self.mapping)
        assert(pairsB.issubset(exclusionsB)), "Error: Not all pairs in topology B are present in the exclusions extracted from the tpr dump. This is not supposed to happen."
        return exclusionsA, pairsA, exclusions_temp, exclusionsB, pairsB
    
    def compute_pairs_exclusions_to_add(self):
        exclusionsA, pairsA, exclusions_temp, exclusionsB, pairsB = self.compute_pairs_exclusions()
        pairsB_minus_A = pairsB.difference(pairsA)
        exclusionsB_minus_A = exclusionsB.difference(exclusionsA).difference(pairsB_minus_A)
        exclusions_temp_minus_A_B = exclusions_temp.difference(exclusionsA).difference(exclusionsB)
        return pairsB_minus_A, exclusionsB_minus_A, exclusions_temp_minus_A_B
    
    def add_pairs_nb_exclusions_to_topology(self):
        pairsB_minus_A, exclusionsB_minus_A, exclusions_temp_minus_A_B = self.compute_pairs_exclusions_to_add()
        self.molA.add_section(add_exclusions_section_to_topology(sorted(pairsB_minus_A.union(exclusionsB_minus_A).union(exclusions_temp_minus_A_B))))
        self.molA.add_section(add_pairs_nb_section_to_topology(pairsB_minus_A, exclusionsB_minus_A.union(exclusions_temp_minus_A_B), self.edge1_steps, self.sigma_epsilon_charges, self.comb_rule, self.fudge_QQ))
        return
    
    def edit_header_sections(self):
        self.topA.replace_header_section_by_name("defaults", edit_defaults_section(self.topA.get_header_section_by_name("defaults")))
        self.topA.replace_header_section_by_name("atomtypes", add_dummy_atomtypes_to_topology(self.topA.get_header_section_by_name("atomtypes")))
        return
    
    def compute_dual_atoms(self):
        dual_state_atoms_to_add = get_dual_state_atoms(self.molA, self.idx_i)
        return dual_state_atoms_to_add
    
    def compute_pairs_to_add(self):
        topB_pairs_to_add = get_topB_pairs_parameters(self.topB, self.pairsB_minus_A, self.mapping)
        return topB_pairs_to_add

    def compute_bonds_to_transform(self):
        harmonic_bondsB_minus_A = map_bonds(extract_harmonic_bonds_from_topology(self.molB), self.mapping).difference(extract_harmonic_bonds_from_topology(self.molA))
        assert(len(self.harmonic_bondsB_minus_A) == 1), "Error: There should be only one harmonic bond in topology B that is not present in topology A. "
        topB_bond_to_add = get_topB_bond_parameters(self.molB, harmonic_bondsB_minus_A, self.mapping)
        return topB_bond_to_add

    def compute_angles_to_add(self):
        anglesB_minus_A = extract_angles_from_topology(self.molB, self.mapping).difference(extract_angles_from_topology(self.molA))
        topB_angles_to_add = get_topB_angles_parameters(self.molB, anglesB_minus_A, self.mapping)
        return topB_angles_to_add

    def compute_dihedrals_to_add(self):
        dihedralsB_minus_A = extract_dihedrals_from_topology(self.molB, self.mapping).difference(extract_dihedrals_from_topology(self.molA))
        topB_dihedrals_to_add = get_topB_dihedrals_parameters(self.molB, dihedralsB_minus_A, self.mapping)
        return topB_dihedrals_to_add

    def write_topology(self):
        topology_writer = Writer(self.topA, "./tests/data/test.top", formatter=GromacsFormatter())
        topology_writer.write_topology()
        return
    
    def replace_topology_sections(self):
        self.molA.replace_section("atoms", updated_atoms_section(self.molA.get_section("atoms"), self.idx_i, self.compute_dual_atoms()))
        self.molA.replace_section("pairs", updated_pairs_section(self.molA.get_section("pairs"), self.idx_i, self.compute_pairs_to_add()))

        topB_bond_to_add = self.compute_bonds_to_transform()
        self.molA.replace_section("bonds", updated_bonds_section(self.molA.get_section("bonds"), topB_bond_to_add, self.idx_i, self.idx_i_minus_1_C, self.idx_i_N, self.idx_i_C, self.idx_i_plus_1_N))
        add_distance_restraint_for_new_bond(self.molA, topB_bond_to_add)

        self.molA.replace_section("angles", updated_angles_section(self.molA.get_section("angles"), self.idx_i, self.compute_angles_to_add()))
        self.molA.replace_section("dihedrals", updated_dihedrals_section(self.molA.get_section("dihedrals"), self.idx_i, self.compute_dihedrals_to_add()))
        return
    
    def write_topology_output(self):
        output_file = "tests/data/test.top"
        writer = Writer(topology=self.topA, file_path=output_file)
        writer.write_topology()
        return
    
    def generate_resdel_topology(self):
        self.get_residue_atom_info()
        self.add_pairs_nb_exclusions_to_topology()
        self.edit_header_sections()
        self.replace_topology_sections()
        self.write_topology_output()
        return


    def __repr__(self):
        return (f"Edge1_Topologies with topology A (topA), topology B (topB), residue to delete: {self.residue_to_delete}, molA name: {self.molA_name}, molB name: {self.molB_name}, edge1 steps: {self.edge1_steps}")
        

    


    
    


    
    

    
        
    
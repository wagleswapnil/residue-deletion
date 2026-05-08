from resdel.topology import *
from typing import Optional, List
from resdel.tranformations.atom_mapping import build_atom_mapping
from resdel.edge1.transformations_bonded import add_peptide_bond
from resdel.edge1.edge1_utils import *
from resdel.edge1.pair_nb import Pair_nb_object

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
        
        idx_i_minus_1, idx_i, idx_i_plus_1, idx_i_minus_1_C, idx_i_plus_1_N = get_residue_atom_idxs(self.molA, self.residue_to_delete)
        self.sigma_epsilon_charges = get_sigma_epsilon_charges(self.topA.get_header_section_by_name("atomtypes"), self.molA, idx_i_minus_1 + idx_i + idx_i_plus_1)
        for key, value in self.sigma_epsilon_charges.items():
            print(f"Atom idx: {key}, sigma: {value[0]}, epsilon: {value[1]}, charge: {value[2]}")
        self.mapping = build_atom_mapping(self.molA.get_section("atoms").lines, self.molB.get_section("atoms").lines, int(self.residue_to_delete))
        
        get_tpr_dump(mdp="./tests/MDP/em.mdp", structure="./tests/data/minimized_stage1.gro", topology="./tests/data/system_stage1.top", output_prefix="./tests/data/system_stage1")
        output_prefix="./tests/data/system_stage1"
        self.exclusionsA = extract_exclusions_from_tpr_dump(f"{output_prefix}.txt", f"{output_prefix}_exclusions.txt")
        #print(f"Extracted exclusions from topology A: {sorted(self.exclusionsA)}")
        self.pairsA = extract_pairs_from_topology(self.topA, self.molA_name)
        #print(f"Extracted pairs from topology A: {sorted(self.pairsA)}")
        assert(self.pairsA.issubset(self.exclusionsA)), "Error: Not all pairs in topology A are present in the exclusions extracted from the tpr dump. This is not supposed to happen."
        #self.exclusionsA.difference_update(self.pairsA)
        #print(f"Extracted exclusions from topology A: {sorted(self.exclusionsA)}")
        #
        add_peptide_bond(self.molA, idx_i_minus_1_C, idx_i_plus_1_N)
        topology_writer = Writer(self.topA, "./tests/data/test.top")
        topology_writer.write_topology()


        get_tpr_dump(mdp="./tests/MDP/em_test.mdp", structure="./tests/data/minimized_stage1.gro", topology="./tests/data/test.top", output_prefix="./tests/data/test")
        output_prefix="./tests/data/test"
        self.exclusions_temp = extract_exclusions_from_tpr_dump(f"{output_prefix}.txt", f"{output_prefix}_exclusions.txt")
        #print(f"Extracted exclusions from topology Temp: {sorted(self.exclusions_temp)}")
        

        #print(f"Exclusions in topology Temp but not in topology A: {sorted(self.exclusions_temp.difference(self.exclusionsA))}")
        #
        get_tpr_dump(mdp="./tests/MDP/em.mdp", structure="./tests/data/minimized_stage5.gro", topology="./tests/data/system_stage5.top", output_prefix="./tests/data/system_stage5")
        output_prefix="./tests/data/system_stage5"
        self.exclusionsB = map_exclusions_pairs(extract_exclusions_from_tpr_dump(f"{output_prefix}.txt", f"{output_prefix}_exclusions.txt"), self.mapping)
        #print(f"Extracted exclusions from topology B: {sorted(self.exclusionsB)}")
        self.pairsB = map_exclusions_pairs(extract_pairs_from_topology(self.topB, self.molB_name), self.mapping)
        #print(f"Extracted pairs from topology B: {sorted(self.pairsB)}")
        assert(self.pairsB.issubset(self.exclusionsB)), "Error: Not all pairs in topology B are present in the exclusions extracted from the tpr dump. This is not supposed to happen."
        #self.exclusionsB.difference_update(self.pairsB)
        #print(f"Extracted exclusions from topology B: {sorted(self.exclusionsB)}")

        self.pairsB_minus_A = self.pairsB.difference(self.pairsA)
        self.exclusionsB_minus_A = self.exclusionsB.difference(self.exclusionsA).difference(self.pairsB_minus_A)
        print(f"Pairs in topology B but not in  A: {sorted(self.pairsB_minus_A)}")
        #print(f"Pairs from topology B not in topology A: {sorted(self.pairsB.difference(self.pairsA))}")
        print(f"Exclusions in topology B not in A: {sorted(self.exclusionsB_minus_A)}")

        self.exclusions_temp_minus_A_B = self.exclusions_temp.difference(self.exclusionsA).difference(self.exclusionsB)
        print(f"Exclusions in topology Temp but not in A or B: {sorted(self.exclusions_temp_minus_A_B)}")
        self.comb_rule, self.fudge_QQ = self.topA.get_comb_rule_fudgeQQ()
    
        self.add_exclusions_section_to_topology()
        self.add_pairs_nb_section_to_topology()
        #topology_writer = Writer(self.topA, "./tests/data/test.top")
        #topology_writer.write_topology()
        
    def __repr__(self):
        return (f"Edge1_Topologies with topology A (topA), topology B (topB), residue to delete: {self.residue_to_delete}, molA name: {self.molA_name}, molB name: {self.molB_name}, edge1 steps: {self.edge1_steps}")
        

    def add_exclusions_section_to_topology(self):
        exclusions_dict = {}
        for x in sorted(self.pairsB_minus_A.union(self.exclusionsB_minus_A).union(self.exclusions_temp_minus_A_B)):
            exclusions_dict[x[0]] = exclusions_dict[x[0]] + " " + str(x[1]) if x[0] in exclusions_dict else str(x[1])
            exclusions_dict[x[1]] = exclusions_dict[x[1]] + " " + str(x[0]) if x[1] in exclusions_dict else str(x[0])
        
        exclusions_section = Section("exclusions")
        exclusions_section.add_line(Line("#ifdef EXCLS_ON"))
        for key, value in sorted(exclusions_dict.items()):
            exclusions_section.add_line(Line(f"{key}  {value}"))
        exclusions_section.add_line(Line("#endif"))
        exclusions_section.add_line(Line(""))
        exclusions_section.add_line(Line(""))
        self.molA.add_section(exclusions_section)
        return


    def add_pairs_nb_section_to_topology(self):
        pairs_nb_list = self.get_pairs_nb_objects()
        pairs_nb_section = Section("pairs_nb")
        for step in range(0, self.edge1_steps):
            pairs_nb_section.add_line(Line(f"#ifdef EDGE1_STEP{step}"))
            for pair_nb in pairs_nb_list:
                pairs_nb_section.add_line(pair_nb.get_pair_nb_line(step / self.edge1_steps))
            pairs_nb_section.add_line(Line("#endif"))
        pairs_nb_section.add_line(Line(f"#ifdef STAGE2"))
        for pair_nb in pairs_nb_list:
            pairs_nb_section.add_line(pair_nb.get_pair_nb_line(step / self.edge1_steps))
        pairs_nb_section.add_line(Line("#endif"))
        pairs_nb_section.add_line(Line(""))
        pairs_nb_section.add_line(Line(""))
        self.molA.add_section(pairs_nb_section)
        return
        
    def get_pairs_nb_objects(self):
        pairs_nb_list = []
        for x in sorted(self.pairsB_minus_A):
            pairs_nb_list.append(Pair_nb_object(x, self.sigma_epsilon_charges[str(x[0])], self.sigma_epsilon_charges[str(x[1])], stateA="full_interactions", stateB="scaled_interactions", comb_rule=self.comb_rule, fudge_QQ=self.fudge_QQ))
        for x in sorted(self.exclusionsB_minus_A):
            pairs_nb_list.append(Pair_nb_object(x, self.sigma_epsilon_charges[str(x[0])], self.sigma_epsilon_charges[str(x[1])], stateA="full_interactions", stateB="no_interactions", comb_rule=self.comb_rule, fudge_QQ=self.fudge_QQ))
        for x in sorted(self.exclusions_temp_minus_A_B):
            pairs_nb_list.append(Pair_nb_object(x, self.sigma_epsilon_charges[str(x[0])], self.sigma_epsilon_charges[str(x[1])], stateA="full_interactions", stateB="no_interactions", comb_rule=self.comb_rule, fudge_QQ=self.fudge_QQ))
        return pairs_nb_list
    


    
    

    
        
    
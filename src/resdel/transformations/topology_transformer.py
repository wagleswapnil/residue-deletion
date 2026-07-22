from resdel.topology import Topology
from typing import Optional
from resdel.transformations.atom_mapping import build_atom_mapping
from resdel.transformations.transformations_utils import *
from resdel.transformations.transformations_nonbonded import *
from resdel.transformations.transformations_bonded import *
from resdel.topology.formatter import GromacsFormatter
from resdel.transformations.tpr_object import TPR_Object

class TopologyTransformer:
    #def __init__(self, topA : Topology, topB : Topology, residue_to_delete : str, molA_name : Optional[str] = "system1", molB_name : Optional[str] = "system1", edge1_steps : Optional[int] = 15):
    def __init__(self, top_wt: Topology, top_mutant: Topology, config, paths):
        self.top_wt = top_wt
        self.top_mutant = top_mutant
        self.config = config
        self.paths = paths
        self.edge1_steps = self.config.transform.number_of_lambdas
        self.residue_to_delete = self.config.system.residue_to_delete
        self.mol_wt = None
        self.mol_mutant = None
        self.mol_wt_name = (config.transform.molecule_name or "system1")
        self.mol_mutant_name = (config.transform.molecule_name or "system1")
        self.config = config
        for mol in self.top_wt.molecules:
            if mol.name == self.mol_wt_name:
                self.mol_wt = mol
                break
        for mol in self.top_mutant.molecules:
            if mol.name == self.mol_mutant_name:
                self.mol_mutant = mol
                break
        if self.mol_wt is None:
            raise ValueError(f"Could not find molecule with name {self.mol_wt_name} in the wt topology")
        if self.mol_mutant is None:
            raise ValueError(f"Could not find molecule with name {self.mol_mutant_name} in the mutant topology")
        self.mapping = build_atom_mapping(self.mol_wt.get_section("atoms").lines, self.mol_mutant.get_section("atoms").lines, int(self.residue_to_delete))
        print(self.mapping)
        self.comb_rule, self.fudge_QQ = self.top_wt.get_comb_rule_fudgeQQ()

    def add_posre_section_to_topology(self):
        for mol in self.top_wt.molecules:
            if mol.name == self.mol_wt_name:
                line = Line(f"")
                mol.sections[-1].add_line(line)
                line = Line(f"#ifdef POSRES")
                mol.sections[-1].add_line(line)
                line = Line(f'#include "posre.itp"')
                mol.sections[-1].add_line(line)
                line = Line(f"#endif")
                mol.sections[-1].add_line(line)
                line = Line(f"")
                mol.sections[-1].add_line(line)
        return

    def get_residue_atom_info(self):
        self.idx_i_minus_1, self.idx_i, self.idx_i_plus_1, self.idx_i_minus_1_N, self.idx_i_minus_1_C, self.idx_i_N, self.idx_i_C, self.idx_i_plus_1_N,  self.idx_i_plus_1_C = get_residue_atom_idxs(self.mol_wt, self.residue_to_delete)
        self.sigma_epsilon_charges = get_sigma_epsilon_charges(self.top_wt.get_header_section_by_name("atomtypes"), self.mol_wt, self.idx_i_minus_1 + self.idx_i + self.idx_i_plus_1)
        return 

    def compute_pairs(self):
        self.pairs_wt = extract_pairs_from_topology(self.top_wt, self.mol_wt_name)
        self.pairs_mutant = map_exclusions_pairs(extract_pairs_from_topology(self.top_mutant, self.mol_mutant_name), self.mapping)
        self.pairs_mutant_minus_wt = self.pairs_mutant.difference(self.pairs_wt)
        return

    def _extract_exclusions_from_tpr(self, topology, structure, output_prefix, molecule_name, mdp : Optional[str] = None):
        tpr = TPR_Object(topology, structure, mdp)
        tpr.make_tpr_dump(output_prefix)
        exclusions = tpr.extract_exclusions_from_tpr_dump(molecule_name)
        return exclusions
    
    def _validate_pairs_subset(self, pairs, exclusions, topology_name):
        if not pairs.issubset(exclusions):
            breakpoint()
            raise ValueError(f"Error: Not all pairs in topology {topology_name} are present in the exclusions extracted from the tpr dump. This is not supposed to happen.")
        return


    def compute_exclusions(self):
        exclusions_wt = self._extract_exclusions_from_tpr(
            structure = self.paths.structure_PDBfile("wt/solvated"),
            topology = self.paths.topology_file("wt/solvated"),
            output_prefix="./tests/data/system_stage1",
            molecule_name=self.mol_wt_name,
            mdp="./tests/MDP/em.mdp"
        )

        self._validate_pairs_subset(self.pairs_wt, exclusions_wt, "wt")

        add_fc0_bond(self.mol_wt, self.idx_i_minus_1_C, self.idx_i_plus_1_N)
            
        topology = "./tests/data/test.top"
        topology_writer = Writer(self.top_wt, topology, formatter=GromacsFormatter())
        topology_writer.write_topology()

        exclusions_temp = self._extract_exclusions_from_tpr(
            structure = self.paths.structure_PDBfile("wt/solvated"),
            topology=topology,
            output_prefix="./tests/data/test",
            molecule_name=self.mol_wt_name,
            mdp="./tests/MDP/em_test.mdp"
        )
        
        exclusions_mutant = self._extract_exclusions_from_tpr(
            structure=self.paths.structure_PDBfile("mutant/solvated"),
            topology=self.paths.topology_file("mutant/solvated"),
            output_prefix="./tests/data/system_stage5",
            molecule_name=self.mol_mutant_name,
            mdp="./tests/MDP/em.mdp"
        )

        exclusions_mutant = map_exclusions_pairs(exclusions_mutant, self.mapping)

        self._validate_pairs_subset(self.pairs_mutant, exclusions_mutant, "mutant")

        return exclusions_wt, exclusions_temp, exclusions_mutant
    
    def compute_pairs_nb_exclusions_to_add(self):
        self.compute_pairs()
        exclusions_wt, exclusions_temp, exclusions_mutant = self.compute_exclusions()
        exclusions_mutant_minus_wt = exclusions_mutant.difference(exclusions_wt).difference(self.pairs_mutant_minus_wt)
        exclusions_temp_minus_wt_mutant = exclusions_temp.difference(exclusions_wt).difference(exclusions_mutant)
        return exclusions_mutant_minus_wt, exclusions_temp_minus_wt_mutant
    
    def add_pairs_nb_exclusions_to_topology(self):
        exclusions_mutant_minus_wt, exclusions_temp_minus_wt_mutant = self.compute_pairs_nb_exclusions_to_add()
        self.mol_wt.add_section(add_exclusions_section_to_topology(sorted(self.pairs_mutant_minus_wt.union(exclusions_mutant_minus_wt).union(exclusions_temp_minus_wt_mutant))))
        self.mol_wt.add_section(add_pairs_nb_section_to_topology(self.pairs_mutant_minus_wt, exclusions_mutant_minus_wt.union(exclusions_temp_minus_wt_mutant), self.edge1_steps, self.sigma_epsilon_charges, self.comb_rule, self.fudge_QQ))
        return
    
    def edit_header_sections(self):
        self.top_wt.replace_header_section_by_name("defaults", edit_defaults_section(self.top_wt.get_header_section_by_name("defaults")))
        self.top_wt.replace_header_section_by_name("atomtypes", add_dummy_atomtypes_to_topology(self.top_wt.get_header_section_by_name("atomtypes")))
        return
    
    def compute_dual_atoms(self):
        dual_state_atoms_to_add, dummy_atoms_to_add = get_dual_state_atoms(self.mol_wt, self.idx_i)
        return dual_state_atoms_to_add, dummy_atoms_to_add

    def compute_pairs_to_add(self):
        top_mutant_pairs_to_add = get_top_mutant_pairs_parameters(self.top_mutant, self.pairs_mutant_minus_wt, self.mapping)
        return top_mutant_pairs_to_add
    
    def compute_bonds_to_transform(self):
        harmonic_bonds_mutant_minus_wt = map_bonds(extract_harmonic_bonds_from_topology(self.mol_mutant), self.mapping).difference(extract_harmonic_bonds_from_topology(self.mol_wt))
        assert(len(harmonic_bonds_mutant_minus_wt) == 1), "Error: There should be only one harmonic bond in topology mutant that is not present in topology wt. "
        top_mutant_bond_to_add = get_top_mutant_bond_parameters(self.mol_mutant, harmonic_bonds_mutant_minus_wt, self.mapping)
        return top_mutant_bond_to_add

    def compute_angles_to_add(self):
        angles_mutant_minus_wt = extract_angles_from_topology(self.mol_mutant, self.mapping).difference(extract_angles_from_topology(self.mol_wt))
        top_mutant_angles_to_add = get_top_mutant_angles_parameters(self.mol_mutant, angles_mutant_minus_wt, self.mapping)
        return top_mutant_angles_to_add

    def compute_dihedrals_to_add(self):
        dihedrals_mutant_minus_wt = extract_dihedrals_from_topology(self.mol_mutant, self.mapping).difference(extract_dihedrals_from_topology(self.mol_wt))
        top_mutant_dihedrals_to_add = get_top_mutant_dihedrals_parameters(self.mol_mutant, dihedrals_mutant_minus_wt, self.mapping)
        return top_mutant_dihedrals_to_add

    def write_topology(self):
        topology_writer = Writer(self.top_wt, "./tests/data/test.top", formatter=GromacsFormatter())
        topology_writer.write_topology()
        return
    
    def replace_topology_sections(self):
        dual_state_atoms_to_add, dummy_atoms_to_add = self.compute_dual_atoms()
        self.mol_wt.replace_section("atoms", updated_atoms_section(self.mol_wt.get_section("atoms"), self.idx_i, dual_state_atoms_to_add, dummy_atoms_to_add))
        self.mol_wt.replace_section("pairs", updated_pairs_section(self.mol_wt.get_section("pairs"), self.idx_i, self.compute_pairs_to_add()))

        top_mutant_bond_to_add = self.compute_bonds_to_transform()
        self.mol_wt.replace_section("bonds", updated_bonds_section(self.mol_wt.get_section("bonds"), top_mutant_bond_to_add, self.idx_i, self.idx_i_minus_1_C, self.idx_i_N, self.idx_i_C, self.idx_i_plus_1_N))
        add_new_bond(self.mol_wt, top_mutant_bond_to_add)

        self.mol_wt.replace_section("angles", updated_angles_section(self.mol_wt.get_section("angles"), self.idx_i, self.compute_angles_to_add()))
        self.mol_wt.replace_section("dihedrals", updated_dihedrals_section(self.mol_wt.get_section("dihedrals"), self.idx_i, self.compute_dihedrals_to_add()))
        return
    
    def write_topology_output(self, output_file):
        writer = Writer(topology=self.top_wt, file_path=output_file, formatter=GromacsFormatter())
        writer.write_topology()
        return
    
    def generate_resdel_topology(self):
        self.get_residue_atom_info()
        self.add_pairs_nb_exclusions_to_topology()
        self.edit_header_sections()
        self.replace_topology_sections()
        #self.write_topology_output()
        return


    def __repr__(self):
        return (f"Edge1_Topologies with topology wt (top_wt), topology mutant (top_mutant), residue to delete: {self.residue_to_delete}, mol_wt name: {self.mol_wt_name}, mol_mutant name: {self.mol_mutant_name}, edge1 steps: {self.edge1_steps}")
        

    


    
    


    
    

    
        
    
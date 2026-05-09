from resdel.topology import *
from resdel.edge1.pair_nb import Pair_nb_object

def add_exclusions_section_to_topology(exclusions_list):
    exclusions_dict = {}
    for x in exclusions_list:
        exclusions_dict[x[0]] = exclusions_dict[x[0]] + " " + str(x[1]) if x[0] in exclusions_dict else str(x[1])
        exclusions_dict[x[1]] = exclusions_dict[x[1]] + " " + str(x[0]) if x[1] in exclusions_dict else str(x[0])
    
    exclusions_section = Section("exclusions")
    exclusions_section.add_line(Line("#ifdef EXCLS_ON"))
    for key, value in sorted(exclusions_dict.items()):
        exclusions_section.add_line(Line(f"{key}  {value}"))
    exclusions_section.add_line(Line("#endif"))
    exclusions_section.add_line(Line(""))
    exclusions_section.add_line(Line(""))
    return exclusions_section           


def get_pairs_nb_objects(pairs_list, exclusions_list, edge1_steps, sigma_epsilon_charges, comb_rule, fudge_QQ):
    pairs_nb_list = []
    for x in sorted(pairs_list):
        pairs_nb_list.append(Pair_nb_object(x, sigma_epsilon_charges[str(x[0])], sigma_epsilon_charges[str(x[1])], stateA="full_interactions", stateB="scaled_interactions", comb_rule=comb_rule, fudge_QQ=fudge_QQ))
    for x in sorted(exclusions_list):
        pairs_nb_list.append(Pair_nb_object(x, sigma_epsilon_charges[str(x[0])], sigma_epsilon_charges[str(x[1])], stateA="full_interactions", stateB="no_interactions", comb_rule=comb_rule, fudge_QQ=fudge_QQ))
    return pairs_nb_list
       

def add_pairs_nb_section_to_topology(pairs_list, exclusions_list, edge1_steps, sigma_epsilon_charges, comb_rule, fudge_QQ):
    pairs_nb_list = get_pairs_nb_objects(pairs_list, exclusions_list, edge1_steps, sigma_epsilon_charges, comb_rule, fudge_QQ)
    pairs_nb_section = Section("pairs_nb")
    for step in range(0, edge1_steps):
        pairs_nb_section.add_line(Line(f"#ifdef EDGE1_STEP{step}"))
        for pair_nb in pairs_nb_list:
            pairs_nb_section.add_line(pair_nb.get_pair_nb_line(step / edge1_steps))
        pairs_nb_section.add_line(Line("#endif"))
    pairs_nb_section.add_line(Line(f"#ifdef STAGE2"))
    for pair_nb in pairs_nb_list:
        pairs_nb_section.add_line(pair_nb.get_pair_nb_line(step / edge1_steps))
    pairs_nb_section.add_line(Line("#endif"))
    pairs_nb_section.add_line(Line(""))
    pairs_nb_section.add_line(Line(""))
    return pairs_nb_section
    

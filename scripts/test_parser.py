from resdel.topology.parser import Parser
from resdel.topology.writer import Writer
from resdel.tranformations import Transformation_Bonded_Params
from typing import Optional


def main():
    
    input_fileA = "tests/data/system_stage1.top"
    topA = Parser(file_path=input_fileA)
    topA.parse_topology()

    input_fileB = "tests/data/system_stage5.top"
    topB = Parser(file_path=input_fileB)
    topB.parse_topology()
    
    output_file = "tests/data/system_stage5_output.top"
    writer = Writer(topology=topB.top, file_path=output_file)
    writer.write_topology()
    print ("Topology parsed and written successfully.")

    transformed_topologies = Transformation_Bonded_Params(topA=topA.top, topB=topB.top, res_to_delete=3)
    transformed_topologies.build_atom_mapping()
    transformed_topologies.tranform_bonds()
    #print(topB.top.molecules[0].get_section("bonds").lines)
if __name__ == "__main__":
    main()
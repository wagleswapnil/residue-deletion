from resdel.topology.parser import Parser
from resdel.topology import write_topology

def main():
    
    input_file = "tests/data/system_stage1.top"
    parser = Parser(file_path=input_file)
    parser.parse_topology()
    breakpoint()  # Debugging point to inspect the parser state after parsing
    output_file = "tests/data/system_stage1_output.top"

    top = parser.parse_topology()
    #breakpoint()  # Debugging point to inspect the parsed topology
    write_topology(top, output_file)

    print ("Topology parsed and written successfully.")

if __name__ == "__main__":
    main()

from resdel.topology import parse_topology, write_topology

def main():
    input_file = "tests/data/system_stage1.top"
    output_file = "tests/data/system_stage1_output.top"

    top = parse_topology(input_file)
    breakpoint()  # Debugging point to inspect the parsed topology
    write_topology(top, output_file)

    print ("Topology parsed and written successfully.")

if __name__ == "__main__":
    main()
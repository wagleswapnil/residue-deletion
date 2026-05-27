rule pdb_to_gro:
    input:
        pdb=config["system"]["pdb"]

    output:
        gro="output/system.gro"

    shell:
        """
        {config[gromacs][gmx]} editconf \
        -f {input.pdb} -o {output.gro} \
        -d 0 -bt cubic
        """

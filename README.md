```
                               ██╗          ██╗  
                               ██║          ██║  
██████╗  █████╗   ██████╗  ██████║  █████╗  ██║  
██╔═══╝ ███████║ ╚█████╗  ██╔══██║ ███████║ ██║  
██║     ╚██████╗ ██████╔╝ ╚██████║ ╚██████╗ ██║  
╚═╝      ╚═════╝ ╚═════╝   ╚═════╝  ╚═════╝ ╚═╝  
░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
                ✦  r e s d e l  ✦
```

# resdel

Prepare GROMACS topology files and analyze free energy calculations for protein residue deletion.

`resdel` is a Python package that prepares GROMACS topology files to run alchemical free energy
calculations for protein residue deletion. The free energy calculation of residue deletion is a
non-trivial process that requires forming a ring structure, and breaking and forming bonds from
that ring structure. The free energy of residue deletion helps protein designers estimate whether
deleting a particular residue is likely to stabilize a protein.

`resdel` takes a protein sequence and the index of the residue to be deleted. It builds the GROMACS
topologies of the wild type protein, the mutant (the protein with the deleted residue), and a topology
to run alchemical free energy simulation for the residue deletion.

> **Status:** This package supports a manuscript in preparation. It is shared for
> review and discussion; please contact the author before using `resdel` in
> published work or building on the method.

---

## Requirements

- Python ≥ 3.10
- GROMACS 2022.1 or later, available on `PATH` as `gmx`
- Network access to `github.com` at install time (see [Dependencies](#dependencies))

On HPC systems you will usually need to load a Python module first:

```bash
module load python/3.10
```

## Installation

```bash
python -m venv .venv
source .venv/bin/activate

pip install git+https://github.com/wagleswapnil/residue-deletion
```

Or from a local clone:

```bash
git clone https://github.com/wagleswapnil/residue-deletion
cd resdel
pip install .
```

Verify the installation:

```bash
resdel --help
```

### Dependencies

`resdel` depends on
[bio2byte-peptidebuilder](https://github.com/Bio2Byte/bio2byte-peptidebuilder),
a fork of PeptideBuilder maintained by the Bio2Byte group that adds geometries
for phosphorylated residues, capping groups and alternative protonation states.
The fork publishes no tagged releases, so it is pinned to a specific commit and
installed from GitHub automatically — installation therefore requires network
access to `github.com`.

## Quick start

All three commands read the same configuration file:

```bash
# 1. Build wild type and mutant topologies and structures
resdel prepare config.yaml

# 2. Build the alchemical topology for the transformation
resdel transform config.yaml

# 3. Estimate the free energy difference from the completed runs
resdel analyze config.yaml
```

A minimal configuration — deleting the central alanine from a capped tripeptide:

```yaml
# config.yaml
system:
  name: AAA_AA
  sequence: "ACE-ALA-ALA-ALA-NME"
  residue_to_delete: 3

io:
  output_dir: resdel_output

prepare:
  GMX_executable: gmx
  forcefield: amber14sb      
  water_model: tip3p             

transform:
  molecule_name: protein
  number_of_lambdas: 20
```

## The three stages

### `prepare`

Creates `output_dir/wt` and `output_dir/mutant`, containing the GROMACS topologies and
PDB structures of the wild type and the mutant proteins, respectively.

### `transform`

Creates `output_dir/transformation`, containing the GROMACS topology and the PDB structure
required to run the alchemical free energy calculation of residue deletion. The free energy 
simulation can be run using script `scripts/run_transformations_solvated.sh`. The GROMACS 
MDP files required to run the free energy calculation will be shared soon. 

### `analyze`

Calculates the free energy difference between the end states (the wild type and mutant proteins)
using the Bennett acceptance ratio (BAR) method, and plots the phase-space overlap between
adjacent intermediate states.

## Configuration

Keys a new user has to set:

| Key | Description | Default |
|-----|-------------|---------|
| `system.sequence` | Amino acid sequence of the wild type protein | — |
| `system.residue_to_delete` | Index of the amino acid to be deleted | — |
| `io.output_dir` | Parent directory for the output | `resdel_output` |
| `prepare.GMX_executable` | GROMACS executable path | `gmx` |
| `prepare.forcefield` | Force field for protein | amber14sb |
| `prepare.water_model` | Water model | tip3p |
| `transform.number_of_lambdas` | Number of intermediate alchemical states | 20 |

## Development

```bash
git clone https://github.com/wagleswapnil/residue-deletion
cd resdel
pip install -e .
```

Run the test suite:

```bash
pytest
```

Tests that require a GROMACS installation are marked and skipped automatically
when `gmx` is not on `PATH`.

## Known limitations

     - Tested on amber14sb/tip3p so far.
     - Deletion of charged residue is not supported.
     - The test suite (pytest) is still being built out

## License

Copyright © 2026 Swapnil Wagle. All rights reserved.

This repository is shared for review purposes. No license to use, copy, modify,
or distribute the code is granted at this time. An open source license will be
applied on publication of the associated manuscript.

## Citation

A manuscript describing this work is in preparation. Please contact the author
before using `resdel` in published work.

The bundled bio2byte-based PeptideBuilder fork is based on:

> M. Z. Tien, D. K. Sydykova, A. G. Meyer, C. O. Wilke (2013). PeptideBuilder:
> A simple Python library to generate model peptides. *PeerJ* 1:e80.

## Acknowledgements

This work was carried out with the support of the
[Mobley Lab](https://mobleylab.org/) at the University of California, Irvine.

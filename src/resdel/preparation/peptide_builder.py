# We are going to build model peptide systems using PeptideBuilder.
# Originally, PeptideBuilder was written by Wilke Lab at UT Austin.  
# But we have found a fork (Bio2byte :: PeptideBuilder) of the original GitHub repo,
# which we will be using here. It supports some additional functionalities from the
# original PeptideBuilder, such as, terminal residues and three letter amino acid codes.

class ModelPeptide:
    def __init__(self, sequence, structure_path, topology_path):
        self.sequence = sequence
        self.structure_path = structure_path
        self.topology_path = topology_path

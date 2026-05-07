from resdel.topology import *
from typing import Optional


class Edge3_Topologies:
    def __init__(self, topA, topB, residue_to_delete : str):
        self.topA = topA
        self.topB = topB
        self.residue_to_delete = residue_to_delete
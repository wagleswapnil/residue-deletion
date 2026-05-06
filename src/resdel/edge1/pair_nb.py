from resdel.topology import Line
from typing import List, Optional

class Pair_nb():
    def __init__(self, idx1, idx2, ftype, q1, q2, sigma, epsilon):
        self.idx1 = idx1
        self.idx2 = idx2
        self.ftype = ftype
        self.q1 = q1
        self.q2 = q2
        self.sigma = sigma
        self.epsilon = epsilon
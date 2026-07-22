import math

from resdel.topology import Line
from typing import Optional

class Pair_nb_object():
    def __init__(self, x, value1, value2, stateA : Optional[str] = "full_interactions", stateB : Optional[str] = "no_interactions", comb_rule : Optional[str] = 2, fudge_QQ : Optional[str] = 1):
        self.idx1 = x[0]
        self.idx2 = x[1]
        self.sigma1, self.epsilon1, self.q1 = float(value1[0]), float(value1[1]), float(value1[2])
        self.sigma2, self.epsilon2, self.q2 = float(value2[0]), float(value2[1]), float(value2[2])
        self.stateA = stateA
        self.stateB = stateB
        self.comb_rule = comb_rule
        self.fudge_QQ = float(fudge_QQ)
        
        self.stateA_q1, self.stateA_q2 = self.calculate_state_charge(self.q1, self.q2, self.stateA)
        self.stateB_q1, self.stateB_q2 = self.calculate_state_charge(self.q1, self.q2, self.stateB)
        self.stateA_epsilon1 = self.calculate_state_epsilon(self.epsilon1, self.epsilon2, self.stateA)
        self.stateB_epsilon2 = self.calculate_state_epsilon(self.epsilon1, self.epsilon2, self.stateB)
        self.stateA_sigma1 = self.calculate_state_sigma(self.sigma1, self.sigma2, self.stateA)
        self.stateB_sigma2 = self.calculate_state_sigma(self.sigma1, self.sigma2, self.stateB)


    def calculate_state_charge(self, q1, q2, state):
        if state == "full_interactions":
            return q1, q2
        elif state == "no_interactions":
            return 0.0, 0.0
        elif state == "scaled_interactions":
            return q1 * math.sqrt(self.fudge_QQ), q2 * math.sqrt(self.fudge_QQ)
        else:
            raise ValueError(f"Unknown state: {state}")
        
    def calculate_state_epsilon(self, epsilon1, epsilon2, state):
        if state == "full_interactions":
            return math.sqrt(epsilon1 * epsilon2)
        elif state == "no_interactions":
            return 0.0
        elif state == "scaled_interactions":
            if self.comb_rule == "2":
                return math.sqrt(epsilon1 * epsilon2) * 0.5
            elif self.comb_rule == "3":
                return math.sqrt(epsilon1 * epsilon2) * 0.5
            else:
                raise ValueError(f"Unknown comb rule: {self.comb_rule}")
        else:
            raise ValueError(f"Unknown state: {state}")
        
    def calculate_state_sigma(self, sigma1, sigma2, state):
        if state == "full_interactions":
            return (sigma1 + sigma2) / 2.0
        elif state == "no_interactions":
            return 0.0
        elif state == "scaled_interactions":
            if self.comb_rule == "2":
                return (sigma1 + sigma2) * 0.5
            elif self.comb_rule == "3":
                return math.sqrt(sigma1 * sigma2)
        else:
            raise ValueError(f"Unknown state: {state}")
        

    def get_pair_nb_line(self, multiplication_factor : Optional[float] = 0):
        q1 = self.stateA_q1 * (1 - multiplication_factor) + self.stateB_q1 * multiplication_factor
        q2 = self.stateA_q2 * (1 - multiplication_factor) + self.stateB_q2 * multiplication_factor
        ftype = "1"
        sigma = self.stateA_sigma1 * (1 - multiplication_factor) + self.stateB_sigma2 * multiplication_factor
        epsilon = self.stateA_epsilon1 * (1 - multiplication_factor) + self.stateB_epsilon2 * multiplication_factor
        return Line(f"{self.idx1}\t{self.idx2}\t{ftype}\t{q1:.6f} \t{q2:.6f}\t {sigma:.6f}\t {epsilon:.6f}")

    def __repr__(self):
        return f"Pair_nb_object(idx1={self.idx1}, idx2={self.idx2}, stateA={self.stateA}, stateB={self.stateB}, stateA_q1={self.stateA_q1}, stateA_q2={self.stateA_q2}, stateB_q1={self.stateB_q1}, stateB_q2={self.stateB_q2}, stateA_epsilon1={self.stateA_epsilon1}, stateB_epsilon2={self.stateB_epsilon2}, stateA_sigma1={self.stateA_sigma1}, stateB_sigma2={self.stateB_sigma2})"
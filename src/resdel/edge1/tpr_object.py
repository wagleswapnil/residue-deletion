from typing import Optional

class TPR_Object:
    def __init__(self, topology, structure, mdp : Optional[str] = "./tests/MDP/em.mdp"):
        self.topology = topology
        self.structure = structure
        self.mdp = mdp

    def get_tpr_dump(self):
        return
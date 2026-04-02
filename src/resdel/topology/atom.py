from typing import Optional, List

class Atom:
    def __init__(self):
        self.index = None
        self.type = None
        self.resnr = None
        self.residue = None
        self.atom = None
        self.cgnr = None
        self.charge = None
        self.mass = None
        self.typeB = None
        self.chargeB = None
        self.massB = None

    def get_atom_info(self, tokens : List[str]):
        if len(tokens) < 8:
            raise ValueError(f"Expected at least 8 tokens for an atom line, got {len(tokens)}: {tokens}")
        elif len(tokens) == 8:
            self.nr, self.type, self.resnr, self.residue, self.atom, self.cgnr, self.charge, self.mass = tokens
        elif len(tokens) == 11:
            self.nr, self.type, self.resnr, self.residue, self.atom, self.cgnr, self.charge, self.mass, self.typeB, self.chargeB, self.massB = tokens

    def __repr__(self):
        return f"Atom(nr={self.nr}, type='{self.type}', resnr={self.resnr}, residue='{self.residue}', atom='{self.atom}', cgnr={self.cgnr}, charge={self.charge}, mass={self.mass}, typeB='{self.typeB}', chargeB={self.chargeB}, massB={self.massB})"
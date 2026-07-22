from .topology import Line
from typing import Optional

class Formatter:
    def format_line(self, line: Line, section_name: Optional[str]=None) -> str:
        return line.raw
        
class GromacsFormatter(Formatter):
    def format_line(self, line: Line, section_name: Optional[str]=None) -> str:
        section_name = section_name.lower() if section_name else None
        
        formatter_map = {
            "defaults": self._format_defaults_line,
            "atomtypes": self._format_atomtypes_line,
            "atoms": self._format_atoms_line,
            "bonds": self._format_bonds_line,
            "pairs": self._format_pairs_line,
            "angles": self._format_angles_line,
            "dihedrals": self._format_dihedrals_line,
            "exclusions": self._format_exclusions_line,
            "pairs_nb": self._format_pairs_nb_line
        }

        if not line.tokens:
            return super().format_line(line, section_name)

        formatter = formatter_map.get(section_name)
        if formatter:
            return formatter(line)
        return super().format_line(line, section_name)
    
    def _format_defaults_line(self, line: Line) -> str:
        if len(line.tokens) == 5:
            nbfunc, comb_rule, gen_pairs, fudgeLJ, fudgeQQ = line.tokens
            return self._append_comment_to_line(f"{nbfunc:<4s}            {comb_rule:<4s}            {gen_pairs:<4s}            {fudgeLJ:<4s} {fudgeQQ:>12s}", line.comment)
        else:
            return super().format_line(line)
            
    def _format_atomtypes_line(self, line: Line) -> str:    
        if len(line.tokens) == 7:
            atype, at_num, mass, charge, ptype, sigma, epsilon = line.tokens
            return self._append_comment_to_line(f"{atype:<8s} {at_num:>6s} {float(mass):10.6f} {float(charge):12.8f} {ptype:>4s} {float(sigma):12.8f} {float(epsilon):12.8f}", line.comment)
        elif len(line.tokens) == 6:
            atype, at_num, mass, ptype, sigma, epsilon = line.tokens
            return self._append_comment_to_line(f"{atype:<8s} {float(at_num):12.8f} {float(mass):10.6f} {ptype:>4s} {float(sigma):12.8f} {float(epsilon):12.8f}", line.comment)
        else:
            return super().format_line(line)

    def _format_atoms_line(self, line: Line) -> str:
        if len(line.tokens) == 8:
            nr, atype, resnr, resname, atomname, cgnr, charge, mass = line.tokens
            return self._append_comment_to_line(f"{nr:>6s} {atype:>10s} {int(resnr):5d} {resname:>6s} {atomname:>6s} {int(cgnr):5d} {float(charge):12.8f} {float(mass):10.6f}", line.comment)
        elif len(line.tokens) == 11:
            nr, atype, resnr, resname, atomname, cgnr, charge, mass, atypeB, chargeB, massB = line.tokens
            return self._append_comment_to_line(f"{nr:>6s} {atype:>10s} {int(resnr):5d} {resname:>6s} {atomname:>6s} {int(cgnr):5d} {float(charge):12.8f} {float(mass):10.6f}  {atypeB:<8s} {float(chargeB):12.8f} {float(massB):10.6f}", line.comment)
        else:
            return super().format_line(line)

    def _format_bonds_line(self, line: Line) -> str:
        if len(line.tokens) == 5:
            ai, aj, funct, r0, fc = line.tokens
            return self._append_comment_to_line(f"{ai:>6s} {aj:>6s} {funct:>4s}  {float(r0):10.6f} {float(fc):15.6f}", line.comment)
        elif len(line.tokens) == 9:
            ai, aj, funct, r0, depth, beta, r0B, depthB, betaB = line.tokens
            return self._append_comment_to_line(f"{ai:>6s} {aj:>6s} {funct:>4s}  {float(r0):10.6f} {float(depth):15.6f}  {float(beta):8.2f}  {float(r0B):10.6f} {float(depthB):15.6f}  {float(betaB):8.2f}", line.comment)
        elif len(line.tokens) == 3:
            ai, aj, funct = line.tokens
            return self._append_comment_to_line(f"{ai:>6s} {aj:>6s} {funct:>4s}", line.comment)
        elif len(line.tokens) == 11:
            ai, aj, funct, r1, r2, rmax, fc, r1B, r2B, rmaxB, fcB = line.tokens
            return self._append_comment_to_line(f"{ai:>6s} {aj:>6s} {funct:>4s}  {float(r1):10.6f} {float(r2):10.6f} {float(rmax):10.6f} {float(fc):15.6f}  {float(r1B):10.6f} {float(r2B):10.6f} {float(rmaxB):10.6f} {float(fcB):15.6f}", line.comment)
        else:
            return super().format_line(line)

    def _format_pairs_line(self, line: Line) -> str:
        if len(line.tokens) == 3:
            ai, aj, funct = line.tokens
            return self._append_comment_to_line(f"{ai:>6s} {aj:>6s} {funct:>4s}", line.comment)
        elif len(line.tokens) == 5:
            ai, aj, funct, sigma, epsilon = line.tokens
            return self._append_comment_to_line(f"{ai:>6s} {aj:>6s} {funct:>4s}  {float(sigma):12.9f} {float(epsilon):12.9f}", line.comment)
        else:
            return super().format_line(line)

    def _format_angles_line(self, line: Line) -> str:
        if len(line.tokens) == 6:
            ai, aj, ak, funct, theta0, fc = line.tokens
            return self._append_comment_to_line(f"{ai:>6s} {aj:>6s} {ak:>6s} {funct:>4s}  {float(theta0):12.7f}  {float(fc):12.6f}", line.comment)
        elif len(line.tokens) == 8:
            ai, aj, ak, funct, theta0, fc, theta0B, fcB = line.tokens
            return self._append_comment_to_line(f"{ai:>6s} {aj:>6s} {ak:>6s} {funct:>4s}  {float(theta0):12.7f}  {float(fc):12.6f}  {float(theta0B):12.7f}  {float(fcB):12.6f}", line.comment)
        else:
            return super().format_line(line)

    def _format_dihedrals_line(self, line: Line) -> str:
        if len(line.tokens) == 8:
            ai, aj, ak, al, funct, phi0, fc, multiplicity = line.tokens
            return self._append_comment_to_line(f"{ai:>6s} {aj:>6s} {ak:>6s} {al:>6s} {funct:>4s}  {float(phi0):12.7f}  {float(fc):12.7f}  {int(multiplicity):2d}", line.comment)
        elif len(line.tokens) == 11:
            ai, aj, ak, al, funct, phi0, fc, multiplicity, phi0B, fcB, multiplicityB = line.tokens
            return self._append_comment_to_line(f"{ai:>6s} {aj:>6s} {ak:>6s} {al:>6s} {funct:>4s}  {float(phi0):12.7f}  {float(fc):12.7f}  {int(multiplicity):2d}  {float(phi0B):12.7f}  {float(fcB):12.7f}  {int(multiplicityB):2d}", line.comment)
        else:
            return super().format_line(line)
        
    def _format_exclusions_line(self, line: Line) -> str:
        if line.tokens:
            return " ".join([f"{token:>4s}" for token in line.tokens])
        else:
            return super().format_line(line)

    def _format_pairs_nb_line(self, line: Line) -> str:
        if len(line.tokens) == 7:
            ai, aj, ftype, q1, q2, sigma, epsilon = line.tokens
            return self._append_comment_to_line(f"{ai:>6s} {aj:>6s} {ftype:>4s} {float(q1):12.8f} {float(q2):12.8f} {float(sigma):12.8f} {float(epsilon):12.8f}", line.comment)
        else:
            return super().format_line(line)


    def _append_comment_to_line(self, line: Line, comment: str) -> str:
        if comment:
            return f"{line}; {comment}" 
        else:
            return line
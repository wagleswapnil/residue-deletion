from openmm import app
import parmed as pmd
from openmm.unit import *

# This function takes the 3D structure of a protein/peptide as input and uses OpenMM to generate its topology in the solvated phase.
# Defult forcefield for protein: amberff14sb, for water: tip3p, default ions: Na+ and Cl-, default padding of the box: 2 nm.
# Defaults pH for side chain protonation: 7.
def create_receptor_system(input_structure):
    ff14sb = app.ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')
    pdbfile = app.PDBFile(input_structure)
    modeller = app.Modeller(pdbfile.topology, pdbfile.positions)
    modeller.addHydrogens(ff14sb, pH=7.0)
    modeller.addSolvent(ff14sb, padding=3*nanometers, model='tip3p', boxShape='cube', ionicStrength=0.2*molar, positiveIon='Na+', negativeIon='Cl-')
    peptide_system = ff14sb.createSystem(modeller.topology, nonbondedMethod=app.PME, rigidWater=False)
    pmd_receptor_struct = pmd.openmm.load_topology(modeller.topology, peptide_system, modeller.positions)
    return pmd_receptor_struct
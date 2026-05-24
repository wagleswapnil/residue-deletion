from openmm import app
import parmed as pmd
from openmm.unit import *

def create_receptor_system(input_structure):
    ff14sb = app.ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')
    pdbfile = app.PDBFile(input_structure)
    modeller = app.Modeller(pdbfile.topology, pdbfile.positions)
    modeller.addHydrogens(ff14sb, pH=7.0)
    modeller.addSolvent(ff14sb, padding=1*nanometers, model='tip3p', ionicStrength=0.2*molar, positiveIon='Na+', negativeIon='Cl-')
    peptide_system = ff14sb.createSystem(modeller.topology, nonbondedMethod=app.PME, rigidWater=False)
    pmd_receptor_struct = pmd.openmm.load_topology(modeller.topology, peptide_system, modeller.positions)
    return pmd_receptor_struct
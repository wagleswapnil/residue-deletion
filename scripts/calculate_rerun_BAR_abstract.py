import os, sys
import numpy as np
import pandas as pd
import pmx
import subprocess
from resdel.analysis.plot_work import plot_work_dist

p1 = os.path.abspath(sys.argv[1])
t1 = os.path.abspath(sys.argv[2])
p2 = os.path.abspath(sys.argv[3])
t2 = os.path.abspath(sys.argv[4])
recalc = int(sys.argv[5])

out_folder = "temp"
input_option = "Potential"

out_folder= os.path.join(os.getcwd(), out_folder)
if not os.path.exists(out_folder):
    os.makedirs(out_folder)

if recalc:
    cmd11 = f"gmx mdrun -s {p1} -rerun {t1} -deffnm p1_t1"
    subprocess.run(cmd11, check=True, shell=True, cwd=out_folder)

    cmd11e = f"gmx energy -f p1_t1.edr -o p1_t1_energy.xvg"
    subprocess.run(cmd11e, check=True, shell=True, text=True, input=input_option, cwd=out_folder)

    cmd12 = f"gmx mdrun -s {p1} -rerun {t2} -deffnm p1_t2"
    subprocess.run(cmd12, check=True, shell=True, cwd=out_folder)

    cmd12e = f"gmx energy -f p1_t2.edr -o p1_t2_energy.xvg"
    subprocess.run(cmd12e, check=True, shell=True, text=True, input=input_option, cwd=out_folder)

    cmd21 = f"gmx mdrun -s {p2} -rerun {t1} -deffnm p2_t1"
    subprocess.run(cmd21, check=True, shell=True, cwd=out_folder)

    cmd21e = f"gmx energy -f p2_t1.edr -o p2_t1_energy.xvg"
    subprocess.run(cmd21e, check=True, shell=True, text=True, input=input_option, cwd=out_folder)

    cmd22 = f"gmx mdrun -s {p2} -rerun {t2} -deffnm p2_t2"
    subprocess.run(cmd22, check=True, shell=True, cwd=out_folder)

    cmd22e = f"gmx energy -f p2_t2.edr -o p2_t2_energy.xvg"
    subprocess.run(cmd22e, check=True, shell=True, text=True, input=input_option, cwd=out_folder)

p1t1 = os.path.join(out_folder, f"p1_t1_energy.xvg")
p1t2 = os.path.join(out_folder, f"p1_t2_energy.xvg")
p2t1 = os.path.join(out_folder, f"p2_t1_energy.xvg")
p2t2 = os.path.join(out_folder, f"p2_t2_energy.xvg")

mp1t1 = pd.read_csv(p1t1, skiprows=23, sep=r'\s+')
mp2t1 = pd.read_csv(p2t1, skiprows=23, sep=r'\s+')

mp1t2 = pd.read_csv(p1t2, skiprows=23, sep=r'\s+')
mp2t2 = pd.read_csv(p2t2, skiprows=23, sep=r'\s+')

rpf = mp2t1['s0']-mp1t1['s0']    # forward work
rpb = mp2t2['s0']-mp1t2['s0']    # backward work

rpf = rpf.dropna()
rpb = rpb.dropna()

rpf = np.array(rpf)
rpb = np.array(rpb)

rpf = rpf[rpf <= 100]
rpb = rpb[rpb <= 100]

energy = pmx.estimators.BAR(rpf, rpb, 298.15)
fe_value = energy.dg * 0.239006 
fe_error = energy.err * 0.239006
print(f"Free energy difference: {fe_value:.3f} ± {fe_error:.3f} kcal/mol")
plot_work_dist(wf=rpf, wr=rpb, nbins=10, dG=energy.dg, dGerr=energy.err, fname=os.path.join(out_folder, "Wdist.png"))

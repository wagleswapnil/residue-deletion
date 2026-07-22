import os

import numpy as np
import pandas as pd
import pmx
from resdel.analysis.plot_work import plot_work_dist

def calculate_fe_from_reduced_potentials(state1, state2, out_folder):
    p1t1 = os.path.join(out_folder, f"{state1}_{state1}_energy.xvg")
    p1t2 = os.path.join(out_folder, f"{state1}_{state2}_energy.xvg")
    p2t1 = os.path.join(out_folder, f"{state2}_{state1}_energy.xvg")
    p2t2 = os.path.join(out_folder, f"{state2}_{state2}_energy.xvg")

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

    energy = pmx.estimators.BAR(rpf, rpb, 298.15)
    fe_value = energy.dg * 0.239006 
    fe_error = energy.err * 0.239006
    print(f"Free energy difference between {state1} and {state2}: {fe_value:.3f} ± {fe_error:.3f} kcal/mol")
    #breakpoint()
    #plot_work_dist(wf=rpf, wr=rpb, nbins=10, dG=energy.dg, dGerr=energy.err, fname=os.path.join(out_folder, f"{state1}_{state2}_Wdist.png"))
    return fe_value, fe_error
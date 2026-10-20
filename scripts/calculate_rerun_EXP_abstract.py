import os, sys
import numpy as np
import pandas as pd
from pymbar.other_estimators import exp
import subprocess
import matplotlib.pyplot as plt

p1 = os.path.abspath(sys.argv[1])
p2 = os.path.abspath(sys.argv[2])
t1 = os.path.abspath(sys.argv[3])
recalc = int(sys.argv[4])

out_folder = "./temp"
input_option = "Potential"

out_folder= os.path.join(os.getcwd(), out_folder)
if not os.path.exists(out_folder):
    os.makedirs(out_folder)

if recalc:
    cmd11 = f"gmx mdrun -s {p1} -rerun {t1} -deffnm p1_t1"
    subprocess.run(cmd11, check=True, shell=True, cwd=out_folder)

    cmd11e = f"gmx energy -f p1_t1.edr -o p1_t1_energy.xvg"
    subprocess.run(cmd11e, check=True, shell=True, text=True, input=input_option, cwd=out_folder)

    cmd21 = f"gmx mdrun -s {p2} -rerun {t1} -deffnm p2_t1"
    subprocess.run(cmd21, check=True, shell=True, cwd=out_folder)

    cmd21e = f"gmx energy -f p2_t1.edr -o p2_t1_energy.xvg"
    subprocess.run(cmd21e, check=True, shell=True, text=True, input=input_option, cwd=out_folder)


p1t1 = os.path.join(out_folder, f"p1_t1_energy.xvg")
p2t1 = os.path.join(out_folder, f"p2_t1_energy.xvg")

mp1t1 = pd.read_csv(p1t1, skiprows=23, sep=r'\s+')
mp2t1 = pd.read_csv(p2t1, skiprows=23, sep=r'\s+')


rpf = mp2t1['s0']-mp1t1['s0']    # forward work


rpf = rpf.dropna()


rpf = np.array(rpf)


energy = exp(rpf, compute_uncertainty=True, is_timeseries=False)
#fe_value = energy.Delta_f * 0.239006 
#fe_error = energy.dDelta_f * 0.239006
print(f"Free energy difference in kcal/mol: {energy}")

temperature = 300.0          # K

R = 0.008314462618           # kJ/mol/K
beta = 1.0 / (R * temperature)

running_dG = np.zeros(len(rpf))

for i in range(1, len(rpf)+1):
    running_dG[i-1] = -(1.0/beta) * np.log(
        np.mean(np.exp(-beta*rpf[:i]))
    )

# Final estimate
final_dG = running_dG[-1]

print(f"Zwanzig free energy = {final_dG:.3f} kJ/mol")

# =====================================
# Plot
# =====================================

plt.figure(figsize=(7,5))

plt.plot(np.arange(1,len(rpf)+1), running_dG,
         lw=2)

plt.axhline(final_dG,
            ls='--',
            color='red',
            label=f'Final = {final_dG:.2f} kJ/mol')

plt.xlabel("Number of frames")
plt.ylabel("ΔG (kJ/mol)")
plt.title("Running exponential averaging (Zwanzig)")
plt.legend()

plt.tight_layout()
plt.savefig("./temp/running_exp_average.png", dpi=300, bbox_inches="tight")

plt.figure(figsize=(6,4))
plt.hist(rpf, bins=50)
plt.xlabel(r'$\Delta U$ (kJ/mol)')
plt.ylabel('Count')
plt.tight_layout()
plt.savefig("./temp/deltaU_histogram.png", dpi=300)

weights = np.exp(-beta*rpf)
weights /= weights.sum()

plt.figure(figsize=(7,4))
plt.plot(weights)
plt.xlabel("Frame")
plt.ylabel("Normalized Boltzmann weight")
plt.tight_layout()
plt.savefig("./temp/weights.png", dpi=300)

ESS = 1.0 / np.sum(weights**2)

print("Effective sample size =", ESS)
print("Fraction of frames =", ESS/len(weights))

np.savetxt(
    "running_dG.dat",
    np.column_stack((np.arange(1, len(rpf)+1), running_dG)),
    header="Frame  Running_dG(kJ/mol)",
    fmt=["%8d", "%12.6f"]
)
from plot_work import plot_work_dist
import numpy as np

rng = np.random.default_rng()
x = rng.integers(low=10, high=50, size=50)
y = rng.integers(low=10, high=50, size=50)

plot_work_dist(wf=x, wr=y, nbins=10, dG=3, dGerr=0.1, fname="/dfs9/dmobley-lab/swapnilw/pdgfra_project/residue_deletion/src/resdel/analysis/test1.png")

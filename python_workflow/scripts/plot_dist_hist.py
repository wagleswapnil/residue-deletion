import glob
import os, sys, re
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

in_path = sys.argv[1]

from scipy.stats import gaussian_kde
import numpy as np

def bhattacharyya_coefficient(x1, x2, npoints=1000):
    """
    Calculate the Bhattacharyya coefficient between two 1D datasets.
    """

    kde1 = gaussian_kde(x1)
    kde2 = gaussian_kde(x2)

    xmin = min(np.min(x1), np.min(x2))
    xmax = max(np.max(x1), np.max(x2))

    grid = np.linspace(xmin, xmax, npoints)

    p = kde1(grid)
    q = kde2(grid)

    bc = np.trapz(np.sqrt(p*q), grid)

    return bc

# Find all dist.xvg files
files = sorted(glob.glob(
    f"{in_path}/edge3_*/production/dist.xvg"
))
files = sorted(
    files,
    key=lambda x: int(re.search(r'edge3_(\d+)', x).group(1))
)
print(f"Found {len(files)} files.")

# ---------- Read all data ----------
col2_data = []
col3_data = []
labels = []

for file in files:
    data = np.loadtxt(file, comments=('#', '@'))

    col2_data.append(data[:, 1])
    col3_data.append(data[:, 2])

    # Extract edge number from path
    edge = os.path.basename(os.path.dirname(os.path.dirname(file)))
    labels.append(edge)

# Common bins
bins2 = np.linspace(
    min(np.min(x) for x in col2_data),
    max(np.max(x) for x in col2_data),
    40
)

bins3 = np.linspace(
    min(np.min(x) for x in col3_data),
    max(np.max(x) for x in col3_data),
    40
)

# ---------- Plot column 2 ----------
plt.figure(figsize=(8,6))

colors = plt.cm.tab20.colors 

plt.figure(figsize=(8,6))

for values, label, color in zip(col2_data, labels, colors):

    kde = gaussian_kde(values)
    x = np.linspace(min(values), max(values), 300)
    y = kde(x)

    plt.fill_between(x, y, alpha=0.3, color=color)
    plt.plot(x, y, color=color, lw=1.5)
    imax = np.argmax(y)

    plt.text(x[imax],
            y[imax],
            label.replace("edge3_", ""),
            fontsize=8,
            ha='center',
            va='bottom',
            color=color)

plt.xlabel("Distance 1 (nm)")
plt.ylabel("Probability density")
plt.tight_layout()
plt.savefig("/dfs9/dmobley-lab/swapnilw/pdgfra_project/residue_deletion/temp/first_third_histograms.png", dpi=300, bbox_inches="tight")

plt.show()
# ---------- Plot column 3 ----------
plt.figure(figsize=(8,6))

colors = plt.cm.Dark2.colors

plt.figure(figsize=(8,6))

for values, label, color in zip(col3_data, labels, colors):

    kde = gaussian_kde(values)
    x = np.linspace(min(values), max(values), 300)
    y = kde(x)

    plt.fill_between(x, y, alpha=0.3, color=color)
    plt.plot(x, y, color=color, lw=1.5)
    imax = np.argmax(y)

    plt.text(x[imax],
            y[imax],
            label.replace("edge3_", ""),
            fontsize=8,
            ha='center',
            va='bottom',
            color=color)

plt.xlabel("Distance 1 (nm)")
plt.ylabel("Probability density")
plt.tight_layout()

plt.savefig("/dfs9/dmobley-lab/swapnilw/pdgfra_project/residue_deletion/temp/second_third_histograms.png", dpi=300, bbox_inches="tight")

plt.show()
plt.show()


fig, axes = plt.subplots(4, 5,
                         figsize=(16, 12),
                         sharex=True,
                         sharey=True)

colors = plt.cm.tab20.colors

for i, (values, label, color) in enumerate(zip(col2_data, labels, colors)):

    ax = axes.flat[i]

    ax.hist(values,
            bins=bins2,
            density=True,
            color=color,
            edgecolor='black',
            linewidth=0.5)

    ax.set_title(label.replace("edge3_", "λ = "), fontsize=10)

    mean = np.mean(values)
    std = np.std(values)

    ax.text(0.05, 0.95,
            f"μ = {mean:.3f}\nσ = {std:.3f}",
            transform=ax.transAxes,
            fontsize=8,
            verticalalignment='top',
            bbox=dict(boxstyle='round',
                    facecolor='white',
                    alpha=0.8,
                    edgecolor='gray'))

# Common labels
fig.supxlabel("Distance 1 (nm)", fontsize=14)
fig.supylabel("Probability density", fontsize=14)

plt.tight_layout()

plt.savefig("/dfs9/dmobley-lab/swapnilw/pdgfra_project/residue_deletion/temp/column2_subplots.png",
            dpi=300,
            bbox_inches="tight")

plt.show()


fig, axes = plt.subplots(4, 5,
                         figsize=(16, 12),
                         sharex=True,
                         sharey=True)

colors = plt.cm.tab20.colors

for i, (values, label, color) in enumerate(zip(col3_data, labels, colors)):

    ax = axes.flat[i]

    ax.hist(values,
            bins=bins3,
            density=True,
            color=color,
            edgecolor='black',
            linewidth=0.5)

    ax.set_title(label.replace("edge3_", "λ = "), fontsize=10)

    mean = np.mean(values)
    std = np.std(values)

    ax.text(0.05, 0.95,
            f"μ = {mean:.3f}\nσ = {std:.3f}",
            transform=ax.transAxes,
            fontsize=8,
            verticalalignment='top',
            bbox=dict(boxstyle='round',
                    facecolor='white',
                    alpha=0.8,
                    edgecolor='gray'))

# Common labels
fig.supxlabel("Distance 1 (nm)", fontsize=14)
fig.supylabel("Probability density", fontsize=14)

plt.tight_layout()

plt.savefig("/dfs9/dmobley-lab/swapnilw/pdgfra_project/residue_deletion/temp/column3_subplots.png",
            dpi=300,
            bbox_inches="tight")


plt.show()


n = len(col2_data)

BC = np.zeros((n,n))

for i in range(n):
    for j in range(n):
        BC[i,j] = bhattacharyya_coefficient(col2_data[i], col2_data[j])

plt.figure(figsize=(7,6))

plt.imshow(BC,
           origin='lower',
           cmap='viridis',
           vmin=0,
           vmax=1)

plt.colorbar(label='Bhattacharyya coefficient')

plt.xlabel("Lambda window")
plt.ylabel("Lambda window")

plt.xticks(range(n))
plt.yticks(range(n))

plt.tight_layout()
plt.savefig("/dfs9/dmobley-lab/swapnilw/pdgfra_project/residue_deletion/temp/BC_matrix_column2.png", dpi=300)
plt.show()

n = len(col3_data)

BC = np.zeros((n,n))

for i in range(n):
    for j in range(n):
        BC[i,j] = bhattacharyya_coefficient(col3_data[i], col3_data[j])

plt.figure(figsize=(7,6))

plt.imshow(BC,
           origin='lower',
           cmap='viridis',
           vmin=0,
           vmax=1)

plt.colorbar(label='Bhattacharyya coefficient')

plt.xlabel("Lambda window")
plt.ylabel("Lambda window")

plt.xticks(range(n))
plt.yticks(range(n))

plt.tight_layout()
plt.savefig("/dfs9/dmobley-lab/swapnilw/pdgfra_project/residue_deletion/temp/BC_matrix_column3.png", dpi=300)
plt.show()


print("Adjacent window overlaps Column 2")
print("------------------------")

for i in range(len(col2_data)-1):
    bc = bhattacharyya_coefficient(col2_data[i], col2_data[i+1])
    print(f"edge3_{i:2d} ↔ edge3_{i+1:2d}: BC = {bc:.3f}")


print("Adjacent window overlaps Column 3")
print("------------------------")

for i in range(len(col3_data)-1):
    bc = bhattacharyya_coefficient(col3_data[i], col3_data[i+1])
    print(f"edge3_{i:2d} ↔ edge3_{i+1:2d}: BC = {bc:.3f}")
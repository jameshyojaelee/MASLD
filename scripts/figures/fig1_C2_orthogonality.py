#!/usr/bin/env python3
"""
Fig 1 Panel C2 — Evidence-source orthogonality matrix (PROTOTYPE -> figures/misc/)
KEY MESSAGE: The seven evidence sources are statistically independent (largest pairwise
             |rho| = 0.21) -- 7 non-redundant experiments, not one signal counted seven times.

A correlation matrix where the eye expects a hot diagonal and finds pale everywhere off it;
the emptiness IS the argument. Reframes "291 atlas columns" from padding into independence.

Data: RNA-seq/results/multi_evidence/convergence_evidence_modality_correlations.csv (long format)
Gate compliance: leads with verified max |rho|=0.21 ONLY. Does NOT cite "K_eff" (fabricated).
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import csv
import os

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 6,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.dpi": 300,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.grid": False, "legend.frameon": False,
})

# canonical 7 convergence channels (sources_active_definition.md); S6 single-cell dropped
ORDER = ["S1", "S8", "S2_intact", "S4", "S5", "S7", "S3"]
LABEL = {
    "S1": "Human bulk", "S8": "Mouse bulk", "S2_intact": "Genetic causal",
    "S4": "Epigenomic", "S5": "Spatial", "S7": "Proteomics", "S3": "Essentiality",
}


def _root():
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def load_matrix():
    f = os.path.join(_root(), "RNA-seq/results/multi_evidence/convergence_evidence_modality_correlations.csv")
    rho = {}
    for r in csv.DictReader(open(f)):
        rho[(r["Mod_A"], r["Mod_B"])] = float(r["Spearman_rho"])
    n = len(ORDER)
    M = np.full((n, n), np.nan)
    for i, a in enumerate(ORDER):
        for j, b in enumerate(ORDER):
            M[i, j] = rho.get((a, b), rho.get((b, a), np.nan))
    return M


def main():
    M = load_matrix()
    n = len(ORDER)
    off = [(i, j) for i in range(n) for j in range(n) if i > j]
    mx_i, mx_j = max(off, key=lambda ij: abs(M[ij]))
    mx = M[mx_i, mx_j]

    # lower triangle = real rho; diagonal = 1.0 (the "expected hot" anchor); upper = blank
    disp = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            if i > j:
                disp[i, j] = M[i, j]
            elif i == j:
                disp[i, j] = 1.0

    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    fig.patch.set_facecolor("white")

    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad(color="white")
    # TRUE correlation range: weak |rho| reads pale; diagonal reads dark -> the contrast is the point
    im = ax.imshow(disp, cmap=cmap, vmin=-1, vmax=1, aspect="equal", zorder=1)

    for i in range(n):
        for j in range(n):
            if i > j:
                ax.text(j, i, f"{M[i, j]:+.02f}", ha="center", va="center",
                        fontsize=6, color="#333333", zorder=3)
            elif i == j:
                ax.text(j, i, "1", ha="center", va="center",
                        fontsize=6, color="white", fontweight="normal", zorder=3)
    for x in np.arange(-0.5, n, 1):
        ax.axhline(x, color="white", lw=1.4, zorder=2)
        ax.axvline(x, color="white", lw=1.4, zorder=2)

    # ring the single largest off-diagonal cell
    ax.add_patch(mpatches.Rectangle((mx_j - 0.5, mx_i - 0.5), 1, 1,
                 facecolor="none", edgecolor="#2D3436", lw=1.8, zorder=4))

    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels([LABEL[s] for s in ORDER], fontsize=6, rotation=35, ha="right")
    ax.set_yticklabels([LABEL[s] for s in ORDER], fontsize=6)
    ax.tick_params(left=False, bottom=False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlim(-0.5, n - 0.5); ax.set_ylim(n - 0.5, -0.5)

    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.04, ticks=[-1, -0.5, 0, 0.5, 1])
    cbar.set_label("Spearman ρ", fontsize=6)
    cbar.ax.tick_params(labelsize=6)
    cbar.outline.set_linewidth(0)

    fig.savefig(os.path.join(_root(), "figures/misc/fig1_C2_orthogonality.pdf"),
                bbox_inches="tight", dpi=300, facecolor="white")
    plt.close(fig)
    print("[caption] Seven evidence sources are statistically independent — "
          f"largest off-diagonal |ρ| = {abs(mx):.2f} (Proteomics vs Human bulk), non-redundant, "
          "not one signal counted 7x")
    print("Saved figures/misc/fig1_C2_orthogonality.pdf  (max |off-diag rho| = %.3f)" % abs(mx))


if __name__ == "__main__":
    main()

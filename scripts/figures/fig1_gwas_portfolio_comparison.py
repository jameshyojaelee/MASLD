#!/usr/bin/env python3
"""
Fig 1 GWAS Portfolio — Options 1 and 2 (for visual comparison)
KEY MESSAGE: 28 GWAS spanning 4 ancestries and 8 liver-disease traits enable
             multi-ancestry causal inference at scale.

Produces two PDFs in figures/misc/ for side-by-side review:
  gwas_option1_tile_matrix.pdf   — ancestry × trait tiles, color = log10(N), count annotation
  gwas_option2_bubble_matrix.pdf — one bubble per GWAS, size = sqrt(N), color = ancestry
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import numpy as np
import os
from collections import defaultdict

# ── Figure defaults (FIGURE_GUIDELINES.md) ───────────────────────────────────
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":         9,
    "pdf.fonttype":      42,
    "ps.fonttype":       42,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
    "savefig.dpi":       300,
    "figure.facecolor":  "white",
    "axes.facecolor":    "white",
    "axes.grid":         False,
    "legend.frameon":    False,
})

# ── GWAS inventory (28 canonical; 17 EUR + 5 EAS + 3 AFR + 3 SAS) ──────────
# (ancestry, trait, N_total, short_label)
GWAS_DATA = [
    # European (17)
    ("EUR", "NAFLD",       8_434, "Namjou 2019"),
    ("EUR", "NAFLD",       9_491, "Anstee 2020"),
    ("EUR", "NAFLD",     778_614, "Ghodsian 2021"),
    ("EUR", "NAFLD",     370_000, "Sveinbj. deCODE"),
    ("EUR", "NAFLD",     111_000, "Sveinbj. Intrm."),
    ("EUR", "NAFLD",     400_000, "Sveinbj. UKBB"),
    ("EUR", "NAFLD",     438_857, "FinnGen R12"),
    ("EUR", "NASH",      435_000, "FinnGen R12"),
    ("EUR", "HCC",       435_000, "FinnGen R12"),
    ("EUR", "ALT",       343_850, "UKBB"),
    ("EUR", "AST",       343_850, "UKBB"),
    ("EUR", "GGT",       343_850, "UKBB"),
    ("EUR", "Cirrhosis", 431_122, "Ghouse 2024"),
    ("EUR", "HCC",       310_000, "Ghouse 2025"),
    ("EUR", "PDFF",       36_116, "Liu 2021"),
    ("EUR", "PDFF",       32_858, "Haas 2021"),
    ("EUR", "PDFF",       44_867, "van der Meer 2022"),
    # East Asian (5)
    ("EAS", "Cirrhosis", 376_326, "Ishigaki 2020"),
    ("EAS", "HCC",       376_326, "Ishigaki 2020"),
    ("EAS", "ALT",       160_000, "BBJ"),
    ("EAS", "AST",       160_000, "BBJ"),
    ("EAS", "GGT",       160_000, "BBJ"),
    # African (3)
    ("AFR", "ALT",         6_636, "Pan-UKBB"),
    ("AFR", "AST",         6_636, "Pan-UKBB"),
    ("AFR", "GGT",         6_636, "Pan-UKBB"),
    # South Asian (3)
    ("SAS", "ALT",         8_876, "Pan-UKBB"),
    ("SAS", "AST",         8_876, "Pan-UKBB"),
    ("SAS", "GGT",         8_876, "Pan-UKBB"),
]

ANCESTRIES = ["EUR", "EAS", "AFR", "SAS"]
ANCESTRY_LABELS = {
    "EUR": "European",
    "EAS": "East Asian",
    "AFR": "African",
    "SAS": "South Asian",
}
ANCESTRY_N = {"EUR": 17, "EAS": 5, "AFR": 3, "SAS": 3}

# Okabe-Ito palette (colorblind-safe)
ANCESTRY_COLORS = {
    "EUR": "#0072B2",
    "EAS": "#E69F00",
    "AFR": "#009E73",
    "SAS": "#CC79A7",
}

TRAITS = ["ALT", "AST", "GGT", "NAFLD", "NASH", "Cirrhosis", "HCC", "PDFF"]

# Trait groups for background bands — (col_start_inclusive, col_end_exclusive)
TRAIT_GROUPS = [
    ("Liver enzymes",       0, 3, "#0072B2"),
    ("Disease endpoints",   3, 7, "#CC79A7"),
    ("Imaging",             7, 8, "#009E73"),
]


def _project_root():
    """scripts/figures/ → scripts/ → project root"""
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _out_path(filename):
    misc = os.path.join(_project_root(), "figures", "misc")
    os.makedirs(misc, exist_ok=True)
    return os.path.join(misc, filename)


def _build_matrix():
    """cell_data[anc][trait] = list of (N, label)"""
    mat = defaultdict(lambda: defaultdict(list))
    for anc, trait, N, label in GWAS_DATA:
        mat[anc][trait].append((N, label))
    return mat


def _add_group_bands(ax, n_trait, alpha=0.07):
    for (_, c0, c1, color) in TRAIT_GROUPS:
        ax.axvspan(c0 - 0.5, c1 - 0.5, color=color, alpha=alpha, zorder=0, linewidth=0)


def _add_group_labels_top(ax):
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    centers = [(c0 + c1 - 1) / 2 for (_, c0, c1, _) in TRAIT_GROUPS]
    labels  = [label for (label, _, _, _) in TRAIT_GROUPS]
    ax2.set_xticks(centers)
    ax2.set_xticklabels(labels, fontsize=7, color="#666666", style="italic")
    ax2.tick_params(top=False, labeltop=True, left=False, right=False, bottom=False)
    for spine in ax2.spines.values():
        spine.set_visible(False)
    return ax2


def _add_ancestry_strips(ax, width=0.07, x_offset=-0.55):
    """Colored strips to the left of ancestry row labels."""
    for i, anc in enumerate(ANCESTRIES):
        ax.add_patch(mpatches.Rectangle(
            (x_offset, i - 0.45), width, 0.90,
            transform=ax.transData,
            color=ANCESTRY_COLORS[anc],
            clip_on=False,
            zorder=5,
        ))


def _format_n(N):
    if N >= 1_000_000:
        return f"{N/1e6:.1f}M"
    if N >= 1_000:
        return f"{N/1e3:.0f}K"
    return str(N)


# ─────────────────────────────────────────────────────────────────────────────
# OPTION 1 — TILE MATRIX
# ─────────────────────────────────────────────────────────────────────────────
def build_option1(out_path):
    """Tile matrix: color = log10(max N per cell), annotate study count + max N."""
    cell_data = _build_matrix()
    n_anc, n_trait = len(ANCESTRIES), len(TRAITS)

    # Build value matrix (log10 max N) and count matrix
    val_mat   = np.full((n_anc, n_trait), np.nan)
    count_mat = np.zeros((n_anc, n_trait), dtype=int)
    for i, anc in enumerate(ANCESTRIES):
        for j, trait in enumerate(TRAITS):
            items = cell_data[anc][trait]
            if items:
                val_mat[i, j]   = np.log10(max(N for N, _ in items))
                count_mat[i, j] = len(items)

    fig, ax = plt.subplots(figsize=(8.5, 3.0))
    fig.patch.set_facecolor("white")

    _add_group_bands(ax, n_trait, alpha=0.06)

    cmap = plt.get_cmap("YlOrRd").copy()
    cmap.set_bad(color="#EBEBEB")

    vmin = np.log10(6_000)    # ~AFR/SAS smallest
    vmax = np.log10(800_000)  # ~largest EUR
    im = ax.imshow(val_mat, cmap=cmap, vmin=vmin, vmax=vmax,
                   aspect="auto", zorder=1, interpolation="nearest")

    # Cell text: study count (top) + max N (center)
    for i in range(n_anc):
        for j in range(n_trait):
            if count_mat[i, j] == 0:
                continue
            # Background brightness determines text color
            norm_val = (val_mat[i, j] - vmin) / (vmax - vmin)
            txt_color = "#222222" if norm_val < 0.65 else "white"

            ax.text(j, i + 0.33, f"n={count_mat[i,j]}",
                    ha="center", va="bottom", fontsize=6.0,
                    color=txt_color, alpha=0.85, zorder=3)

            max_n = max(N for N, _ in cell_data[ANCESTRIES[i]][TRAITS[j]])
            ax.text(j, i - 0.08, _format_n(max_n),
                    ha="center", va="center", fontsize=7.5,
                    color=txt_color, fontweight="bold", zorder=3)

    # White cell-dividing lines
    for x in np.arange(-0.5, n_trait, 1):
        ax.axvline(x, color="white", linewidth=1.2, zorder=2)
    for y in np.arange(-0.5, n_anc, 1):
        ax.axhline(y, color="white", linewidth=1.2, zorder=2)

    # Axes
    ax.set_xticks(range(n_trait))
    ax.set_xticklabels(TRAITS, fontsize=8.5)
    ax.set_yticks(range(n_anc))
    ax.set_yticklabels(
        [f"{ANCESTRY_LABELS[a]} ({ANCESTRY_N[a]})" for a in ANCESTRIES],
        fontsize=8.5,
    )
    ax.tick_params(left=False, bottom=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    _add_ancestry_strips(ax)
    _add_group_labels_top(ax)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, pad=0.015, fraction=0.022, aspect=20)
    cbar.set_label("log₁₀(max N)", fontsize=7.5)
    cbar.ax.tick_params(labelsize=7)
    cbar.outline.set_linewidth(0)
    # Add reference ticks
    ref_ns = [10_000, 100_000, 500_000]
    cbar.set_ticks([np.log10(n) for n in ref_ns])
    cbar.set_ticklabels([_format_n(n) for n in ref_ns])

    ax.set_title(
        "Option 1 — Tile matrix  |  color = max N per cell  |  n = # studies  |  bold = max N",
        fontsize=8, pad=16, loc="left", color="#555555",
    )

    fig.tight_layout(rect=[0, 0, 1, 1])
    fig.savefig(out_path, bbox_inches="tight", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# OPTION 2 — BUBBLE MATRIX
# ─────────────────────────────────────────────────────────────────────────────
def build_option2(out_path):
    """Bubble matrix: one bubble per GWAS, size ∝ sqrt(N), color = ancestry."""
    cell_data = _build_matrix()
    n_anc, n_trait = len(ANCESTRIES), len(TRAITS)

    all_N   = [N for _, _, N, _ in GWAS_DATA]
    max_N   = max(all_N)
    MAX_PTS = 450   # scatter s-units for largest bubble

    fig, ax = plt.subplots(figsize=(8.5, 3.6))
    fig.patch.set_facecolor("white")

    _add_group_bands(ax, n_trait, alpha=0.07)

    # Light grid
    for x in np.arange(-0.5, n_trait, 1):
        ax.axvline(x, color="#E8E8E8", linewidth=0.8, zorder=1)
    for y in np.arange(-0.5, n_anc, 1):
        ax.axhline(y, color="#E8E8E8", linewidth=0.8, zorder=1)

    # Draw bubbles
    xs, ys, ss, cs = [], [], [], []
    for i, anc in enumerate(ANCESTRIES):
        for j, trait in enumerate(TRAITS):
            items = cell_data[anc][trait]
            if not items:
                continue
            n_items = len(items)
            # Horizontal spread within cell; cap so bubbles don't leave the cell
            if n_items == 1:
                offsets = [0.0]
            else:
                spread = min(0.38 * (n_items - 1) / max(n_items - 1, 1), 0.40)
                offsets = np.linspace(-spread, spread, n_items)

            for (N, _), dx in zip(items, offsets):
                xs.append(j + dx)
                ys.append(i)
                ss.append(MAX_PTS * (N / max_N))   # area ∝ N (not sqrt, more legible)
                cs.append(ANCESTRY_COLORS[anc])

    ax.scatter(xs, ys, s=ss, c=cs, alpha=0.72, zorder=3, linewidths=0.3,
               edgecolors="white")

    # Axes
    ax.set_xlim(-0.5, n_trait - 0.5)
    ax.set_ylim(-0.5, n_anc - 0.5)
    ax.set_xticks(range(n_trait))
    ax.set_xticklabels(TRAITS, fontsize=8.5)
    ax.set_yticks(range(n_anc))
    ax.set_yticklabels(
        [f"{ANCESTRY_LABELS[a]} ({ANCESTRY_N[a]})" for a in ANCESTRIES],
        fontsize=8.5,
    )
    ax.tick_params(left=False, bottom=False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    _add_ancestry_strips(ax, width=0.06)
    _add_group_labels_top(ax)

    # Size legend (bottom-right inset)
    ref_ns  = [10_000, 100_000, 400_000]
    ref_s   = [MAX_PTS * (n / max_N) for n in ref_ns]
    leg_h   = [ax.scatter([], [], s=s, c="#888888", alpha=0.72, linewidths=0,
                          label=_format_n(n))
               for s, n in zip(ref_s, ref_ns)]
    size_leg = ax.legend(
        handles=leg_h, title="N", title_fontsize=7,
        fontsize=6.5, loc="lower right",
        handletextpad=0.4, labelspacing=0.5,
        bbox_to_anchor=(1.0, -0.05),
        frameon=False,
    )
    ax.add_artist(size_leg)

    # Ancestry color legend
    anc_patches = [
        mpatches.Patch(color=ANCESTRY_COLORS[a], alpha=0.75,
                       label=f"{ANCESTRY_LABELS[a]} ({ANCESTRY_N[a]})")
        for a in ANCESTRIES
    ]
    ax.legend(handles=anc_patches, fontsize=7, ncol=2,
              loc="lower left", bbox_to_anchor=(0.0, -0.22),
              frameon=False, handlelength=1.2)

    ax.set_title(
        "Option 2 — Bubble matrix  |  each bubble = one GWAS  |  area ~ N  |  color = ancestry",
        fontsize=8, pad=16, loc="left", color="#555555",
    )

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    build_option1(_out_path("gwas_option1_tile_matrix.pdf"))
    build_option2(_out_path("gwas_option2_bubble_matrix.pdf"))
    print("Done.")

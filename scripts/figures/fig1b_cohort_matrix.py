#!/usr/bin/env python3
"""
Figure 1b — Cohort × Modality dot/tile matrix.

Rows   = 9 human bulk RNA-seq cohorts, ordered largest → smallest.
Columns = 5 modalities (GWAS/eQTL, Bulk RNA-seq, scRNA-seq, Spatial, Proteomics).
Filled circle  = cohort has that data.
Open  circle   = cohort does not have that data.
★ beside cohort label = control-bearing (contributes to mega-analysis, n=846 total).

Output: figures/misc/fig1b_cohort_matrix.pdf
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os, sys

# ── Canonical palette ────────────────────────────────────────────────────────
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig1_palette import MODALITY_COLORS, CONTROL_GRAY

plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":        8,
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
    "figure.dpi":       150,
    "axes.facecolor":   "white",
    "figure.facecolor": "white",
})

# ── Output helpers ────────────────────────────────────────────────────────────
def _root():
    return os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    )

def _out(fn):
    d = os.path.join(_root(), "figures", "misc")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, fn)

# ── Data definitions ──────────────────────────────────────────────────────────

# Cohorts: (pi_label, geo_id, n_samples, is_control_bearing)
# Ordered largest → smallest
COHORTS = [
    ("Chen",        "GSE213621", 358, True),
    ("Govaere",     "GSE135251", 215, True),
    ("Hoshida",     "GSE193066", 160, False),
    ("Bril",        "GSE162694", 142, True),
    ("Kozumi",      "GSE167523",  97, False),
    ("Kawamura",    "GSE174478",  93, False),
    ("Hoang",       "GSE130970",  76, True),
    ("Verschuren",  "GSE240729",  64, False),
    ("Suppli",      "GSE126848",  55, True),
]

# Modalities: (key, label, color)
MODALITIES = [
    ("gwas",    "GWAS/eQTL",    MODALITY_COLORS["gwas"]),
    ("bulk",    "Bulk RNA-seq", MODALITY_COLORS["bulk"]),
    ("scrna",   "scRNA-seq",    MODALITY_COLORS["scrna"]),
    ("spatial", "Spatial",      MODALITY_COLORS["spatial"]),
    ("proteo",  "Proteomics",   MODALITY_COLORS["proteo"]),
]

# Filled-circle matrix: (cohort_index, modality_key)
# GWAS/eQTL: all 9 (population-level)
# Bulk RNA-seq: all 9
# Spatial (Visium GSE192741): Govaere only
# Proteomics (Olink, 177 staged subjects): Govaere only
# scRNA: none of the 9 main cohorts directly (shown via footer note)
COHORT_NAMES = [c[0] for c in COHORTS]
MOD_KEYS     = [m[0] for m in MODALITIES]

FILLED = {(name, mkey): False for name in COHORT_NAMES for mkey in MOD_KEYS}

# GWAS + Bulk: all 9
for name in COHORT_NAMES:
    FILLED[(name, "gwas")] = True
    FILLED[(name, "bulk")] = True

# Spatial + Proteomics: Govaere only
FILLED[("Govaere", "spatial")] = True
FILLED[("Govaere", "proteo")]  = True


# ── Layout constants ──────────────────────────────────────────────────────────
N_ROWS    = len(COHORTS)     # 9
N_COLS    = len(MODALITIES)  # 5

FIG_W     = 6.0   # inches
FIG_H     = 5.0   # inches

DOT_R     = 120   # scatter marker size (s=, in points²)
OPEN_LW   = 1.0   # linewidth for open circles
OPEN_ALPHA = 0.30  # alpha for open circles

COL_LABEL_Y_OFFSET = 0.55   # fraction of a row above top row
ROW_GRID_COLOR     = "#F5F5F5"
STAR_COLOR         = "#C62828"  # marker color for control-bearing cohorts


def main():
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    fig.subplots_adjust(left=0.28, right=0.88, top=0.88, bottom=0.10)

    # ── Draw horizontal background stripes (every other row) ─────────────────
    for ri in range(N_ROWS):
        if ri % 2 == 1:
            ax.axhspan(ri - 0.5, ri + 0.5,
                       facecolor="#FAFAFA", edgecolor="none", zorder=0)

    # ── Draw thin horizontal grid lines between rows ──────────────────────────
    for ri in range(N_ROWS + 1):
        ax.axhline(ri - 0.5, color=ROW_GRID_COLOR, linewidth=0.5, zorder=1)

    # ── Plot dots ─────────────────────────────────────────────────────────────
    for ci, (mkey, mlabel, mcolor) in enumerate(MODALITIES):
        for ri, (pi_name, geo_id, n, is_ctrl) in enumerate(COHORTS):
            # y = 0 is top row (Chen), y = N_ROWS-1 is bottom (Suppli)
            row_y = ri

            if FILLED[(pi_name, mkey)]:
                ax.scatter(ci, row_y, s=DOT_R, color=mcolor,
                           edgecolors="white", linewidths=0.5,
                           zorder=5)
            else:
                ax.scatter(ci, row_y, s=DOT_R,
                           facecolors="none",
                           edgecolors=mcolor,
                           linewidths=OPEN_LW,
                           alpha=OPEN_ALPHA,
                           zorder=4)

    # ── Row labels (left side) ────────────────────────────────────────────────
    for ri, (pi_name, geo_id, n, is_ctrl) in enumerate(COHORTS):
        # Cohort PI name (italic)
        ax.text(-0.55, ri, pi_name,
                ha="right", va="center",
                fontsize=8.5, fontstyle="italic",
                fontweight="semibold",
                color="#212121")
        # GEO ID + n (smaller, gray)
        ax.text(-0.55, ri + 0.28,
                f"{geo_id}  n={n}",
                ha="right", va="center",
                fontsize=6.0, color="#757575")
        # Marker for control-bearing cohorts (bold asterisk, Helvetica-safe)
        if is_ctrl:
            ax.text(N_COLS - 0.40, ri, "*",
                    ha="left", va="center",
                    fontsize=10, fontweight="bold",
                    color=STAR_COLOR, zorder=6)

    # ── Column headers ────────────────────────────────────────────────────────
    for ci, (mkey, mlabel, mcolor) in enumerate(MODALITIES):
        ax.text(ci, -0.78, mlabel,
                ha="center", va="bottom",
                fontsize=8, fontweight="bold",
                color=mcolor, rotation=0)

    # ── Bracket / annotation for mega-analysis cohorts ───────────────────────
    # Control-bearing rows: Chen(0), Govaere(1), Bril(3), Hoang(6), Suppli(8)
    ctrl_rows = [ri for ri, (_, _, _, is_ctrl) in enumerate(COHORTS) if is_ctrl]
    r_min, r_max = min(ctrl_rows), max(ctrl_rows)
    bx = N_COLS - 0.10   # x position just right of dots

    # Draw bracket
    bracket_x = bx + 0.18
    bracket_color = STAR_COLOR
    ax.annotate("",
                xy=(bracket_x, r_max + 0.45),
                xytext=(bracket_x, r_min - 0.45),
                arrowprops=dict(arrowstyle="-",
                                color=bracket_color,
                                lw=1.2),
                annotation_clip=False)
    # Tick lines at top and bottom of bracket
    for row_val in [r_min - 0.45, r_max + 0.45]:
        ax.plot([bracket_x - 0.05, bracket_x + 0.05], [row_val, row_val],
                color=bracket_color, lw=1.2, clip_on=False)

    # Label
    mid_y = (r_min + r_max) / 2
    ax.text(bracket_x + 0.16, mid_y,
            "* control-bearing\n(n=846 mega-analysis)",
            ha="left", va="center",
            fontsize=6.0, color=STAR_COLOR,
            linespacing=1.4,
            clip_on=False)

    # ── Vertical separator after GWAS column ─────────────────────────────────
    ax.axvline(0.55, color="#E0E0E0", linewidth=0.8, linestyle="--",
               zorder=2, ymin=0.02, ymax=0.98)

    # ── Footer note ───────────────────────────────────────────────────────────
    fig.text(0.28, 0.025,
             "† scRNA-seq: 7 separate datasets, 895,542 cells (269 donors); "
             "not directly linked to the 9 bulk cohorts above.",
             ha="left", va="bottom",
             fontsize=5.5, color="#757575",
             style="italic")

    # ── Axes cosmetics ────────────────────────────────────────────────────────
    ax.set_xlim(-0.55, N_COLS - 0.5 + 1.6)
    ax.set_ylim(N_ROWS - 0.5, -1.0)   # invert: row 0 (Chen) at top

    ax.axis("off")

    # ── Legend (bottom-right, inside the panel) ───────────────────────────────
    legend_x = N_COLS * 0.10 + 0.5
    legend_y = N_ROWS + 0.25
    # Filled dot
    ax.scatter(legend_x, legend_y, s=60, color="#555555",
               edgecolors="white", linewidths=0.4, zorder=10, clip_on=False)
    ax.text(legend_x + 0.18, legend_y, "Data available",
            ha="left", va="center", fontsize=6.5, color="#555555")
    # Open dot
    ax.scatter(legend_x + 1.5, legend_y, s=60,
               facecolors="none", edgecolors="#555555",
               linewidths=1.0, alpha=0.5, zorder=10, clip_on=False)
    ax.text(legend_x + 1.68, legend_y, "Not available",
            ha="left", va="center", fontsize=6.5, color="#555555")

    # ── Save ─────────────────────────────────────────────────────────────────
    out_path = _out("fig1b_cohort_matrix.pdf")
    fig.savefig(out_path, bbox_inches="tight", dpi=300, facecolor="white")
    print(f"Saved: {out_path}")
    plt.close(fig)


if __name__ == "__main__":
    main()

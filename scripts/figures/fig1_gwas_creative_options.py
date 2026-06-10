#!/usr/bin/env python3
"""
Fig 1 GWAS Portfolio — Creative Options A–D (breadth-focused, no N encoding)
KEY MESSAGE: 23 GWAS spanning 4 ancestries x liver-disease traits — show coverage breadth.

Outputs to figures/misc/:
  gwas_creative_A_binary_grid.pdf
  gwas_creative_B_alluvial.pdf
  gwas_creative_C_radial.pdf
  gwas_creative_D_unit_chart.pdf
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
import numpy as np
import os
from collections import defaultdict

# ── rcParams ──────────────────────────────────────────────────────────────────
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

# ── Data ──────────────────────────────────────────────────────────────────────
GWAS_DATA = [
    ("EUR", "NAFLD",     8_434),
    ("EUR", "NAFLD",     9_491),
    ("EUR", "NAFLD",   778_614),
    ("EUR", "NAFLD",   370_000),
    ("EUR", "NAFLD",   111_000),
    ("EUR", "NAFLD",   400_000),
    ("EUR", "NAFLD",   438_857),
    ("EUR", "NASH",    435_000),
    ("EUR", "ALT",     343_850),
    ("EUR", "AST",     343_850),
    ("EUR", "GGT",     343_850),
    ("EUR", "PDFF",     36_116),
    ("EUR", "PDFF",     32_858),
    ("EUR", "PDFF",     44_867),
    ("EAS", "ALT",     160_000),
    ("EAS", "AST",     160_000),
    ("EAS", "GGT",     160_000),
    ("AFR", "ALT",       6_636),
    ("AFR", "AST",       6_636),
    ("AFR", "GGT",       6_636),
    ("SAS", "ALT",       8_876),
    ("SAS", "AST",       8_876),
    ("SAS", "GGT",       8_876),
]

ANCESTRIES     = ["EUR", "EAS", "AFR", "SAS"]
ANCESTRY_LABEL = {"EUR": "European", "EAS": "East Asian", "AFR": "African", "SAS": "South Asian"}
ANCESTRY_N     = {"EUR": 14, "EAS": 3, "AFR": 3, "SAS": 3}

# Okabe-Ito palette (colorblind-safe)
ANCESTRY_COLOR = {
    "EUR": "#0072B2", "EAS": "#E69F00", "AFR": "#009E73", "SAS": "#CC79A7",
}

TRAITS = ["ALT", "AST", "GGT", "NAFLD", "NASH", "PDFF"]
# Trait group: (label, col_start, col_end_excl, band_color)
TRAIT_GROUPS = [
    ("Liver enzymes",     0, 3, "#9ECAE1"),
    ("Disease diagnoses", 3, 5, "#FCBBA1"),
    ("Imaging",           5, 6, "#A1D99B"),
]

ABSENT_COLOR = "#E8E8E8"


def _cell_data():
    mat = defaultdict(lambda: defaultdict(int))
    for anc, trait, _ in GWAS_DATA:
        mat[anc][trait] += 1
    return mat


def _project_root():
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _out_path(fname):
    d = os.path.join(_project_root(), "figures", "misc")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, fname)


def _group_bands(ax, alpha=0.07):
    for _, c0, c1, color in TRAIT_GROUPS:
        ax.axvspan(c0 - 0.5, c1 - 0.5, color=color, alpha=alpha, zorder=0, linewidth=0)


def _group_top_labels(ax):
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks([(c0 + c1 - 1) / 2 for _, c0, c1, _ in TRAIT_GROUPS])
    ax2.set_xticklabels([lbl for lbl, *_ in TRAIT_GROUPS],
                        fontsize=7, color="#666666", style="italic")
    ax2.tick_params(top=False, labeltop=True, left=False, right=False, bottom=False)
    for sp in ax2.spines.values():
        sp.set_visible(False)


def _anc_strips(ax, xoff=-0.55, w=0.07):
    for i, anc in enumerate(ANCESTRIES):
        ax.add_patch(mpatches.Rectangle(
            (xoff, i - 0.45), w, 0.9,
            transform=ax.transData, color=ANCESTRY_COLOR[anc],
            clip_on=False, zorder=6,
        ))


def _grid_axis(ax):
    n_trait = len(TRAITS)
    ax.set_xlim(-0.5, n_trait - 0.5)
    ax.set_ylim(-0.5, len(ANCESTRIES) - 0.5)
    ax.set_xticks(range(n_trait))
    ax.set_xticklabels(TRAITS, fontsize=8.5)
    ax.set_yticks(range(len(ANCESTRIES)))
    ax.set_yticklabels([ANCESTRY_LABEL[a] for a in ANCESTRIES], fontsize=8.5)
    ax.tick_params(left=False, bottom=False)
    for sp in ax.spines.values():
        sp.set_visible(False)


# ─────────────────────────────────────────────────────────────────────────────
# A: Binary presence grid
# ─────────────────────────────────────────────────────────────────────────────
def build_A(path):
    """Binary grid: cell = ancestry color if covered, gray if absent. Count annotated."""
    cd = _cell_data()
    fig, ax = plt.subplots(figsize=(8.5, 3.0))
    fig.patch.set_facecolor("white")
    _group_bands(ax, alpha=0.05)

    for i, anc in enumerate(ANCESTRIES):
        for j, trait in enumerate(TRAITS):
            n = cd[anc][trait]
            color = ANCESTRY_COLOR[anc] if n > 0 else ABSENT_COLOR
            ax.add_patch(mpatches.FancyBboxPatch(
                (j - 0.43, i - 0.43), 0.86, 0.86,
                boxstyle="round,pad=0.04",
                facecolor=color, edgecolor="white", linewidth=1.0,
                alpha=0.80 if n > 0 else 1.0, zorder=2,
            ))
            if n > 1:
                ax.text(j, i, str(n), ha="center", va="center",
                        fontsize=10, color="white", fontweight="bold", zorder=3)

    for x in np.arange(-0.5, len(TRAITS), 1):
        ax.axvline(x, color="white", linewidth=1.5, zorder=1)
    for y in np.arange(-0.5, len(ANCESTRIES), 1):
        ax.axhline(y, color="white", linewidth=1.5, zorder=1)

    _grid_axis(ax)
    _anc_strips(ax)
    _group_top_labels(ax)

    leg = [mpatches.Patch(facecolor=ANCESTRY_COLOR[a], alpha=0.80,
                          label=f"{ANCESTRY_LABEL[a]} ({ANCESTRY_N[a]} studies)")
           for a in ANCESTRIES]
    leg.append(mpatches.Patch(facecolor=ABSENT_COLOR, label="Not covered"))
    ax.legend(handles=leg, fontsize=7, loc="lower right",
              bbox_to_anchor=(1.0, -0.26), ncol=3, frameon=False)

    ax.set_title("Option A — Binary coverage grid  |  number = study count  |  color = ancestry",
                 fontsize=8, pad=16, loc="left", color="#555555")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Saved: {path}")


# ─────────────────────────────────────────────────────────────────────────────
# B: Alluvial flow  (1 unit = 1 GWAS study)
# ─────────────────────────────────────────────────────────────────────────────
def _ribbon(ax, x0, y0b, y0t, x1, y1b, y1t, color, alpha=0.40):
    cx = (x0 + x1) / 2
    verts = [
        (x0, y0b), (cx, y0b), (cx, y1b), (x1, y1b),
        (x1, y1t), (cx, y1t), (cx, y0t), (x0, y0t),
        (x0, y0b),
    ]
    codes = [Path.MOVETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.CLOSEPOLY]
    ax.add_patch(mpatches.PathPatch(
        Path(verts, codes),
        facecolor=color, edgecolor="none", alpha=alpha, zorder=2,
    ))


def build_B(path, extra_paths=None):
    """Alluvial: unit width = 1 GWAS study. Left = ancestries, Right = traits."""
    cd = _cell_data()
    GAP = 0.55
    BW  = 0.055   # block width

    # Trait order: most GWAS first (descending), keeps largest ribbons on top → less crossing
    trait_total = {t: sum(cd[a][t] for a in ANCESTRIES) for t in TRAITS}
    right_order = sorted(TRAITS, key=lambda t: -trait_total[t])

    # Ancestry blocks (left, top-to-bottom, y increasing downward after invert_yaxis)
    left_start = {}
    y = 0
    for anc in ANCESTRIES:
        left_start[anc] = y
        y += ANCESTRY_N[anc] + GAP
    total_left = y - GAP

    # Trait blocks (right, same order as right_order)
    right_start = {}
    y = 0
    for t in right_order:
        right_start[t] = y
        y += trait_total[t] + GAP
    total_right = y - GAP

    total_h = max(total_left, total_right) + 1.0
    # Extra bottom margin to accommodate the annotation tag
    TAG_MARGIN = 1.5

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    fig.patch.set_facecolor("white")
    ax.set_xlim(-0.28, 1.28)
    ax.set_ylim(0, total_h + TAG_MARGIN)
    ax.invert_yaxis()
    ax.axis("off")

    # Left ancestry blocks + labels
    for anc in ANCESTRIES:
        y0 = left_start[anc]
        h  = ANCESTRY_N[anc]
        ax.add_patch(mpatches.Rectangle(
            (-BW / 2, y0), BW, h,
            facecolor=ANCESTRY_COLOR[anc], edgecolor="none", alpha=0.92, zorder=4,
        ))
        ax.text(-BW / 2 - 0.02, y0 + h / 2,
                f"{ANCESTRY_LABEL[anc]}\nn={ANCESTRY_N[anc]}",
                ha="right", va="center", fontsize=8.5,
                color=ANCESTRY_COLOR[anc], fontweight="bold")

    # Right trait blocks + labels  (colored by trait group)
    grp_color_of = {}
    for (_, c0, c1, gc) in TRAIT_GROUPS:
        for t in TRAITS[c0:c1]:
            grp_color_of[t] = gc

    for t in right_order:
        y0 = right_start[t]
        h  = trait_total[t]
        ax.add_patch(mpatches.Rectangle(
            (1.0 - BW / 2, y0), BW, h,
            facecolor=grp_color_of[t], edgecolor="none", alpha=0.92, zorder=4,
        ))
        ax.text(1.0 + BW / 2 + 0.02, y0 + h / 2,
                f"{t}  n={h}", ha="left", va="center",
                fontsize=8.5, color="#444444")

    # Ribbons: iterate traits in right_order within each ancestry block
    l_off = {anc: left_start[anc] for anc in ANCESTRIES}
    r_off = {t: right_start[t] for t in right_order}

    for anc in ANCESTRIES:
        for t in right_order:
            n = cd[anc][t]
            if n == 0:
                continue
            _ribbon(ax, 0.0, l_off[anc], l_off[anc] + n,
                        1.0, r_off[t],   r_off[t] + n,
                    ANCESTRY_COLOR[anc], alpha=0.38)
            l_off[anc] += n
            r_off[t]   += n

    # Category legend (right side top)
    cat_patches = [mpatches.Patch(facecolor=gc, alpha=0.92, label=lbl)
                   for lbl, _, _, gc in TRAIT_GROUPS]
    ax.legend(handles=cat_patches, fontsize=7, loc="upper right",
              bbox_to_anchor=(1.26, 1.0), frameon=False)

    # Annotation tag: just below the bottom of the existing plot nodes
    # After invert_yaxis(), total_h is visually at the bottom; we place tag slightly below it.
    # Use DejaVu Sans for this text element to ensure the arrow glyph renders correctly.
    ax.text(0.5, total_h + 0.85,
            "23 GWAS → 368 colocalized genes",
            ha="center", va="top",
            fontsize=7.5, color="#2D3436", style="italic",
            fontfamily="DejaVu Sans",
            transform=ax.transData)

    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white")
    print(f"Saved: {path}")
    if extra_paths:
        for ep in extra_paths:
            fig.savefig(ep, bbox_inches="tight", dpi=300, facecolor="white")
            print(f"Saved: {ep}")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# C: Radial coverage portrait
# ─────────────────────────────────────────────────────────────────────────────
def build_C(path):
    """Radial: outer ring = each trait arc split into 4 ancestry segments (filled/absent).
       Inner ring = trait category colors. Center = study count."""
    cd = _cell_data()

    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    ax.set_aspect("equal")
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.55, 1.55)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    R_OUT     = 1.05   # outer ancestry ring outer radius
    R_IN      = 0.62   # outer ancestry ring inner radius
    R_CAT_OUT = 0.55   # inner category ring outer radius
    R_CAT_IN  = 0.30   # inner category ring inner radius
    R_CENTER  = R_CAT_IN - 0.01

    N_TRAITS      = len(TRAITS)
    DEG_PER_TRAIT = 360.0 / N_TRAITS
    TRAIT_GAP     = 2.0   # degrees between traits
    ANC_GAP       = 0.5   # degrees between ancestry segments within a trait

    # Inner category ring
    cat_angle = {"Liver enzymes": (0, 135), "Disease endpoints": (135, 315), "Imaging": (315, 360)}
    for lbl, c0, c1, gc in TRAIT_GROUPS:
        t1, t2 = c0 * DEG_PER_TRAIT, c1 * DEG_PER_TRAIT
        ax.add_patch(mpatches.Wedge((0, 0), R_CAT_OUT, t1, t2,
                                    width=R_CAT_OUT - R_CAT_IN,
                                    facecolor=gc, edgecolor="white",
                                    linewidth=0.8, alpha=0.80, zorder=1))
        # Category label
        mid_rad = np.radians((t1 + t2) / 2)
        r_mid = (R_CAT_OUT + R_CAT_IN) / 2
        lx, ly = r_mid * np.cos(mid_rad), r_mid * np.sin(mid_rad)
        short = lbl.replace("Disease endpoints", "Disease\nendpoints").replace("Liver enzymes", "Liver\nenzymes")
        ax.text(lx, ly, short, ha="center", va="center",
                fontsize=5.5, color="#444444", style="italic", zorder=5,
                multialignment="center")

    # Outer ancestry coverage ring
    for i, trait in enumerate(TRAITS):
        t_start = i * DEG_PER_TRAIT + TRAIT_GAP / 2
        t_end   = (i + 1) * DEG_PER_TRAIT - TRAIT_GAP / 2
        eff     = t_end - t_start
        n_anc   = len(ANCESTRIES)
        anc_span = (eff - (n_anc - 1) * ANC_GAP) / n_anc

        for k, anc in enumerate(ANCESTRIES):
            a1 = t_start + k * (anc_span + ANC_GAP)
            a2 = a1 + anc_span
            covered = cd[anc][trait] > 0
            color = ANCESTRY_COLOR[anc] if covered else ABSENT_COLOR
            alpha = 0.85 if covered else 0.90
            ax.add_patch(mpatches.Wedge((0, 0), R_OUT, a1, a2,
                                        width=R_OUT - R_IN,
                                        facecolor=color, edgecolor="white",
                                        linewidth=0.7, alpha=alpha, zorder=2))

        # Trait label outside ring
        mid_deg = i * DEG_PER_TRAIT + DEG_PER_TRAIT / 2
        mid_rad = np.radians(mid_deg)
        lx = (R_OUT + 0.17) * np.cos(mid_rad)
        ly = (R_OUT + 0.17) * np.sin(mid_rad)
        # Text alignment based on position
        ha = "center"
        if np.cos(mid_rad) > 0.4:
            ha = "left"
        elif np.cos(mid_rad) < -0.4:
            ha = "right"
        va = "center"
        if np.sin(mid_rad) > 0.4:
            va = "bottom"
        elif np.sin(mid_rad) < -0.4:
            va = "top"
        ax.text(lx, ly, trait, ha=ha, va=va,
                fontsize=8.5, color="#222222", fontweight="bold", zorder=5)

    # White center circle
    ax.add_patch(mpatches.Circle((0, 0), R_CENTER,
                                  facecolor="white", edgecolor="#CCCCCC",
                                  linewidth=0.6, zorder=4))
    ax.text(0, 0.02, "23", ha="center", va="center",
            fontsize=14, color="#333333", fontweight="bold", zorder=6)
    ax.text(0, -0.13, "GWAS", ha="center", va="center",
            fontsize=7, color="#777777", zorder=6)

    # Ancestry legend
    anc_leg = [mpatches.Patch(facecolor=ANCESTRY_COLOR[a], alpha=0.85,
                               label=f"{ANCESTRY_LABEL[a]} ({ANCESTRY_N[a]})")
               for a in ANCESTRIES]
    absent_leg = mpatches.Patch(facecolor=ABSENT_COLOR, label="Not covered")
    ax.legend(handles=anc_leg + [absent_leg], fontsize=7.5,
              loc="lower center", bbox_to_anchor=(0.5, -0.07),
              ncol=3, frameon=False, handlelength=1.2)

    ax.set_title("Option C — Radial coverage  |  outer ring: ancestry × trait  |  filled = covered",
                 fontsize=8, pad=8, loc="center", color="#555555")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Saved: {path}")


# ─────────────────────────────────────────────────────────────────────────────
# D: Unit chart  (stacked squares)
# ─────────────────────────────────────────────────────────────────────────────
def build_D(path):
    """Unit chart: each small square = 1 GWAS study, stacked bottom-up per cell."""
    cd = _cell_data()
    n_anc, n_trait = len(ANCESTRIES), len(TRAITS)

    SQ  = 0.115   # square side in data units
    GAP = 0.018   # gap between stacked squares

    fig, ax = plt.subplots(figsize=(8.5, 3.8))
    fig.patch.set_facecolor("white")
    _group_bands(ax, alpha=0.05)

    for i, anc in enumerate(ANCESTRIES):
        for j, trait in enumerate(TRAITS):
            n = cd[anc][trait]
            if n == 0:
                # Absent placeholder
                ax.add_patch(mpatches.FancyBboxPatch(
                    (j - 0.43, i - 0.43), 0.86, 0.86,
                    boxstyle="round,pad=0.04",
                    facecolor=ABSENT_COLOR, edgecolor="white",
                    linewidth=0.8, alpha=0.6, zorder=1,
                ))
                continue
            cell_bot = i - 0.44
            for k in range(n):
                sq_x = j - SQ / 2
                sq_y = cell_bot + k * (SQ + GAP)
                ax.add_patch(mpatches.FancyBboxPatch(
                    (sq_x, sq_y), SQ, SQ,
                    boxstyle="round,pad=0.012",
                    facecolor=ANCESTRY_COLOR[anc],
                    edgecolor="white", linewidth=0.6,
                    alpha=0.82, zorder=2,
                ))

    for x in np.arange(-0.5, n_trait, 1):
        ax.axvline(x, color="white", linewidth=1.8, zorder=0)
    for y in np.arange(-0.5, n_anc, 1):
        ax.axhline(y, color="white", linewidth=1.8, zorder=0)

    _grid_axis(ax)
    _anc_strips(ax)
    _group_top_labels(ax)

    # Legend
    leg = [mpatches.Patch(facecolor=ANCESTRY_COLOR[a], alpha=0.82,
                          label=f"{ANCESTRY_LABEL[a]} ({ANCESTRY_N[a]})")
           for a in ANCESTRIES]
    leg.append(mpatches.Patch(facecolor=ABSENT_COLOR, label="Not covered"))
    ax.legend(handles=leg, fontsize=7, loc="lower right",
              bbox_to_anchor=(1.0, -0.28), ncol=3, frameon=False)

    ax.set_title("Option D — Unit chart  |  each square = 1 GWAS study  |  color = ancestry",
                 fontsize=8, pad=16, loc="left", color="#555555")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Saved: {path}")


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    build_A(_out_path("gwas_creative_A_binary_grid.pdf"))
    build_B(_out_path("gwas_creative_B_alluvial.pdf"),
            extra_paths=[_out_path("fig1c_gwas_alluvial.pdf")])
    print("All done.")

#!/usr/bin/env python3
"""
Fig 1 GWAS Portfolio — v3 (beauty pass: no absence drawn)
KEY MESSAGE: 28 GWAS across 4 ancestries x 8 liver traits — a fully-colored
             coverage portrait where richness of color = cross-ancestry breadth.

Outputs to figures/misc/:
  gwas_v3_donut.pdf     — fully-filled coverage donut (no gray); color richness = breadth
  gwas_v3_backbone.pdf  — cross-ancestry backbone vs European depth (clean rebuild)
  gwas_v3_skyline.pdf   — stacked ancestry "skyline" per trait (no gray)
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
from collections import defaultdict

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

GWAS_DATA = [
    ("EUR", "NAFLD",     8_434), ("EUR", "NAFLD",     9_491),
    ("EUR", "NAFLD",   778_614), ("EUR", "NAFLD",   370_000),
    ("EUR", "NAFLD",   111_000), ("EUR", "NAFLD",   400_000),
    ("EUR", "NAFLD",   438_857), ("EUR", "NASH",    435_000),
    ("EUR", "HCC",     435_000), ("EUR", "ALT",     343_850),
    ("EUR", "AST",     343_850), ("EUR", "GGT",     343_850),
    ("EUR", "Cirrhosis", 431_122), ("EUR", "HCC",   310_000),
    ("EUR", "PDFF",     36_116), ("EUR", "PDFF",     32_858),
    ("EUR", "PDFF",     44_867),
    ("EAS", "Cirrhosis", 376_326), ("EAS", "HCC",   376_326),
    ("EAS", "ALT",     160_000), ("EAS", "AST",     160_000),
    ("EAS", "GGT",     160_000),
    ("AFR", "ALT",       6_636), ("AFR", "AST",       6_636),
    ("AFR", "GGT",       6_636),
    ("SAS", "ALT",       8_876), ("SAS", "AST",       8_876),
    ("SAS", "GGT",       8_876),
]

ANCESTRIES     = ["EUR", "EAS", "AFR", "SAS"]
ANCESTRY_LABEL = {"EUR": "European", "EAS": "East Asian", "AFR": "African", "SAS": "South Asian"}
ANCESTRY_N     = {"EUR": 17, "EAS": 5, "AFR": 3, "SAS": 3}
ANCESTRY_COLOR = {"EUR": "#0072B2", "EAS": "#E69F00", "AFR": "#009E73", "SAS": "#CC79A7"}
ANC_ALPHA = 0.80   # softer than the old crisp 0.92 (matches the alluvial's gentler feel)

TRAIT_CAT = {
    "ALT": "Liver enzymes", "AST": "Liver enzymes", "GGT": "Liver enzymes",
    "NAFLD": "Disease diagnoses", "NASH": "Disease diagnoses",
    "Cirrhosis": "Disease diagnoses", "HCC": "Disease diagnoses",
    "PDFF": "Imaging",
}
# pastel trait-category accents, matched to the alluvial (gwas_creative_B)
CATEGORY_COLOR = {
    "Liver enzymes":     "#9ECAE1",
    "Disease diagnoses": "#FCBBA1",
    "Imaging":           "#A1D99B",
}


def _cd():
    mat = defaultdict(lambda: defaultdict(int))
    for anc, trait, _ in GWAS_DATA:
        mat[anc][trait] += 1
    return mat


def _project_root():
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _out(fname):
    d = os.path.join(_project_root(), "figures", "misc")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, fname)


def _covs(cd, t):
    """Covering ancestries in canonical order."""
    return [a for a in ANCESTRIES if cd[a][t] > 0]


def _anc_count(cd, t):
    return len(_covs(cd, t))


def _study_count(cd, t):
    return sum(cd[a][t] for a in ANCESTRIES)


def _ranked_traits(cd):
    return sorted(TRAIT_CAT.keys(),
                  key=lambda t: (-_anc_count(cd, t), -_study_count(cd, t)))


# ─────────────────────────────────────────────────────────────────────────────
# DONUT — fully-filled coverage portrait (no gray)
# ─────────────────────────────────────────────────────────────────────────────
def build_donut(path):
    cd = _cd()
    traits = _ranked_traits(cd)   # colorful (4-anc) → solid (1-anc)
    N = len(traits)

    fig, ax = plt.subplots(figsize=(6.6, 6.6))
    ax.set_aspect("equal")
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.5, 1.42)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    R_OUT, R_IN = 1.12, 0.66
    R_CAT_OUT, R_CAT_IN = 0.64, 0.46   # inner trait-category ring (alluvial pastels)
    R_CENTER = 0.44
    GAP = 1.6   # degrees between trait wedges
    START = 90  # start at top, go clockwise
    # Wedge angle ∝ # ancestries covering the trait → every ancestry segment is one equal
    # "coverage unit"; cross-ancestry traits get bigger arcs, European-only ones become slivers.
    U = sum(_anc_count(cd, t) for t in traits)
    unit = (360.0 - N * GAP) / U
    bounds, _ang = {}, START
    for _t in traits:
        _span = _anc_count(cd, _t) * unit
        bounds[_t] = (_ang - _span, _ang)
        _ang -= _span + GAP

    for trait in traits:
        t_lo, t_hi = bounds[trait]
        covs = _covs(cd, trait)
        m = len(covs)
        seg = (t_hi - t_lo) / m
        for k, anc in enumerate(covs):
            a1 = t_lo + k * seg
            ax.add_patch(mpatches.Wedge(
                (0, 0), R_OUT, a1, a1 + seg, width=R_OUT - R_IN,
                facecolor=ANCESTRY_COLOR[anc], edgecolor="white",
                linewidth=1.0, alpha=ANC_ALPHA, zorder=2))
        # inner trait-category ring (pastel), spanning the whole wedge
        ax.add_patch(mpatches.Wedge(
            (0, 0), R_CAT_OUT, t_lo, t_hi, width=R_CAT_OUT - R_CAT_IN,
            facecolor=CATEGORY_COLOR[TRAIT_CAT[trait]], edgecolor="white",
            linewidth=0.8, alpha=0.90, zorder=2))

        # trait label at one consistent outer radius; # GWAS sits inside the ancestry wedge
        mid = np.radians((t_lo + t_hi) / 2)
        lx, ly = (R_OUT + 0.20) * np.cos(mid), (R_OUT + 0.20) * np.sin(mid)
        ha = "left" if np.cos(mid) > 0.35 else ("right" if np.cos(mid) < -0.35 else "center")
        va = "bottom" if np.sin(mid) > 0.35 else ("top" if np.sin(mid) < -0.35 else "center")
        ax.text(lx, ly, trait, ha=ha, va=va, fontsize=9.5,
                color="#222222", fontweight="bold", zorder=5)
        r_mid = (R_OUT + R_IN) / 2
        ax.text(r_mid * np.cos(mid), r_mid * np.sin(mid), str(_study_count(cd, trait)),
                ha="center", va="center", fontsize=7, color="white",
                fontweight="bold", zorder=5)

    # center hub (hero number only)
    ax.add_patch(mpatches.Circle((0, 0), R_CENTER, facecolor="white",
                                 edgecolor="#E0E0E0", linewidth=0.8, zorder=4))
    ax.text(0, 0.07, "28", ha="center", va="center", fontsize=21,
            color="#2D3436", fontweight="bold", zorder=6)
    ax.text(0, -0.15, "GWAS", ha="center", va="center", fontsize=9,
            color="#636E72", zorder=6, fontweight="bold")

    # ancestry legend + pastel category legend (keys for the two color dimensions)
    anc_leg = [mpatches.Patch(facecolor=ANCESTRY_COLOR[a], alpha=ANC_ALPHA,
                              label=f"{ANCESTRY_LABEL[a]} ({ANCESTRY_N[a]})")
               for a in ANCESTRIES]
    cat_leg = [mpatches.Patch(facecolor=CATEGORY_COLOR[c], alpha=0.90, label=c)
               for c in ["Liver enzymes", "Disease diagnoses", "Imaging"]]
    leg1 = ax.legend(handles=anc_leg, fontsize=8, loc="lower center",
                     bbox_to_anchor=(0.5, -0.01), ncol=4, frameon=False,
                     handlelength=1.1, columnspacing=1.2)
    ax.add_artist(leg1)
    ax.legend(handles=cat_leg, fontsize=7.5, loc="lower center",
              bbox_to_anchor=(0.5, -0.085), ncol=3, frameon=False,
              handlelength=1.1, columnspacing=1.2)
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Saved: {path}")


# ─────────────────────────────────────────────────────────────────────────────
# BACKBONE — clean rebuild (no label collisions)
# ─────────────────────────────────────────────────────────────────────────────
def build_backbone(path):
    cd = _cd()
    traits = _ranked_traits(cd)
    n_cross = sum(1 for t in traits if _anc_count(cd, t) >= 2)
    n_rows = len(traits)

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    fig.patch.set_facecolor("white")

    LAB_X   = -0.10                     # trait label column (right-aligned here)
    SLOT_X0 = 0.50                      # first ancestry slot
    SLOT_DX = 0.85
    RIGHT_X = SLOT_X0 + 3 * SLOT_DX     # last ancestry slot
    SQ = 0.50

    y_of = {t: (n_rows - 1 - r) for r, t in enumerate(traits)}
    y_div = (n_rows - 1 - n_cross) + 0.5

    # per-row background tinted by trait category (alluvial pastels)
    for t in traits:
        ax.axhspan(y_of[t] - 0.5, y_of[t] + 0.5,
                   color=CATEGORY_COLOR[TRAIT_CAT[t]], alpha=0.30, zorder=0, linewidth=0)
    # divider between cross-ancestry and European-only zones
    ax.axhline(y_div, color="#9AA0A6", lw=0.8, ls=(0, (4, 3)), zorder=1)

    # ancestry column headers
    for j, anc in enumerate(ANCESTRIES):
        ax.text(SLOT_X0 + j * SLOT_DX, n_rows - 0.30, anc, ha="center", va="bottom",
                fontsize=7.5, color=ANCESTRY_COLOR[anc], fontweight="bold")
    ax.text(RIGHT_X + 0.5, n_rows - 0.30, "anc·GWAS", ha="left", va="bottom",
            fontsize=6.0, color="#9AA0A6", style="italic")

    for t in traits:
        y = y_of[t]
        for j, anc in enumerate(ANCESTRIES):
            x = SLOT_X0 + j * SLOT_DX
            if cd[anc][t] > 0:
                ax.add_patch(mpatches.FancyBboxPatch(
                    (x - SQ / 2, y - SQ / 2), SQ, SQ, boxstyle="round,pad=0.02",
                    facecolor=ANCESTRY_COLOR[anc], edgecolor="white",
                    linewidth=0.8, alpha=ANC_ALPHA, zorder=3))
                if cd[anc][t] > 1:
                    ax.text(x, y, str(cd[anc][t]), ha="center", va="center",
                            fontsize=7, color="white", fontweight="bold", zorder=4)
            else:
                ax.add_patch(mpatches.Circle((x, y), 0.05, facecolor="none",
                             edgecolor="#BBBBBB", lw=1.0, zorder=3))
        ax.text(LAB_X, y, t, ha="right", va="center", fontsize=8.5,
                color="#222222", fontweight="bold")
        ax.text(RIGHT_X + 0.5, y, f"{_anc_count(cd, t)}·{_study_count(cd, t)}",
                ha="left", va="center", fontsize=7, color="#777777")

    # zone labels (rotated, left gutter)
    ax.text(-0.80, (y_div + n_rows - 0.5) / 2, "CROSS-ANCESTRY",
            ha="center", va="center", fontsize=7, color="#0072B2",
            fontweight="bold", rotation=90)
    ax.text(-0.80, (y_div - 0.6) / 2, "European only",
            ha="center", va="center", fontsize=7, color="#666666",
            fontweight="bold", rotation=90)

    # pastel category legend
    cat_leg = [mpatches.Patch(facecolor=CATEGORY_COLOR[c], alpha=0.55, label=c)
               for c in ["Liver enzymes", "Disease diagnoses", "Imaging"]]
    ax.legend(handles=cat_leg, fontsize=6.8, loc="lower center",
              bbox_to_anchor=(0.5, -0.11), ncol=3, frameon=False, handlelength=1.1)

    ax.set_xlim(-1.05, RIGHT_X + 1.2)
    ax.set_ylim(-0.95, n_rows + 0.05)
    ax.axis("off")
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white", pad_inches=0.03)
    plt.close(fig)
    print(f"Saved: {path}")


# ─────────────────────────────────────────────────────────────────────────────
# SKYLINE — stacked ancestry blocks per trait (no gray)
# ─────────────────────────────────────────────────────────────────────────────
def build_skyline(path):
    cd = _cd()
    traits = _ranked_traits(cd)
    n_trait = len(traits)

    fig, ax = plt.subplots(figsize=(8.0, 3.3))
    fig.patch.set_facecolor("white")

    BW = 0.70
    for j, t in enumerate(traits):
        covs = _covs(cd, t)
        for k, anc in enumerate(covs):
            ax.add_patch(mpatches.FancyBboxPatch(
                (j - BW / 2, k + 0.06), BW, 0.88, boxstyle="round,pad=0.01",
                facecolor=ANCESTRY_COLOR[anc], edgecolor="white",
                linewidth=1.0, alpha=ANC_ALPHA, zorder=2))
        # pastel category strip beneath each bar
        ax.add_patch(mpatches.FancyBboxPatch(
            (j - BW / 2, -0.40), BW, 0.24, boxstyle="round,pad=0.01",
            facecolor=CATEGORY_COLOR[TRAIT_CAT[t]], edgecolor="white",
            linewidth=0.6, alpha=0.85, zorder=2))
        ax.text(j, len(covs) + 0.16, f"{_study_count(cd, t)} GWAS",
                ha="center", va="bottom", fontsize=6.8, color="#888888")

    ax.set_xlim(-0.6, n_trait - 0.4)
    ax.set_ylim(-0.7, 4.7)
    ax.set_xticks(range(n_trait))
    ax.set_xticklabels(traits, fontsize=9, fontweight="bold")
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_ylabel("Ancestries with GWAS", fontsize=9)
    ax.tick_params(bottom=False, labelsize=8)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_linewidth(0.4)
    ax.spines["bottom"].set_visible(False)

    # both legends tucked into the empty upper-right (above the short right-hand bars)
    anc_leg = [mpatches.Patch(facecolor=ANCESTRY_COLOR[a], alpha=ANC_ALPHA,
                              label=f"{ANCESTRY_LABEL[a]} ({ANCESTRY_N[a]})")
               for a in ANCESTRIES]
    cat_leg = [mpatches.Patch(facecolor=CATEGORY_COLOR[c], alpha=0.85, label=c)
               for c in ["Liver enzymes", "Disease diagnoses", "Imaging"]]
    leg1 = ax.legend(handles=anc_leg, fontsize=7.5, loc="upper right",
                     bbox_to_anchor=(1.0, 1.0), ncol=2, frameon=False, handlelength=1.1)
    ax.add_artist(leg1)
    ax.legend(handles=cat_leg, fontsize=7, loc="upper right",
              bbox_to_anchor=(1.0, 0.76), ncol=1, frameon=False, handlelength=1.1)

    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white", pad_inches=0.03)
    plt.close(fig)
    print(f"Saved: {path}")


if __name__ == "__main__":
    build_donut(_out("gwas_v3_donut.pdf"))
    build_backbone(_out("gwas_v3_backbone.pdf"))
    build_skyline(_out("gwas_v3_skyline.pdf"))
    print("All done.")

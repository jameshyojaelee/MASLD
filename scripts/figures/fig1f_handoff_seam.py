#!/usr/bin/env python3
"""
Fig 1f — Handoff Seam: genetics → bulk RNA-seq atlas
KEY MESSAGE: 368 SuSiE-COLOC effector genes (PP4>0.5) are 1.28x enriched among
             1,885 Tier-1 DEGs (continuous Wilcoxon p=7.2e-51).
             Modest overlap is expected: coding-effect variants (PNPLA3) are
             excluded by construction (no cis-eQTL).

Output: figures/misc/fig1f_handoff_seam.pdf

All numbers are HARDCODED from verified facts pack:
  - 368  SuSiE-COLOC genes (PP4 > 0.5, 23 GWAS)
  - 67   EUR + EAS replicated
  - 1885 Tier-1 DEGs (padj < 0.05, |LFC| > 0.5)
  - 1.28x enrichment (continuous Wilcoxon, p = 7.2e-51)
  - 27187 genes x 315 columns in atlas
  - 9 cohorts, n = 846
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import numpy as np
import os

# ── rcParams ─────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":        "sans-serif",
    "font.sans-serif":    ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":          9,
    "pdf.fonttype":       42,
    "ps.fonttype":        42,
    "figure.dpi":         150,
    "savefig.dpi":        300,
    "figure.facecolor":   "white",
    "axes.facecolor":     "white",
    "axes.grid":          False,
    "legend.frameon":     False,
})

# ── Canonical palette (single source of truth) ────────────────────────────────
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig1_palette import MODALITY_COLORS, CONTROL_GRAY

GWAS_COL  = MODALITY_COLORS["gwas"]   # #000000 genetics (GWAS/eQTL) - black
BULK_COL  = MODALITY_COLORS["bulk"]   # #0072B2 bulk RNA-seq           - Okabe-Ito blue
ATAC_COL  = MODALITY_COLORS["atac"]   # #E69F00 ATAC / epigenomic       - orange
C_GRAY    = "#636E72"
C_LGRAY   = "#BDBDBD"
C_TEXT    = "#171A1C"


def _root():
    return os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    )


def _blend(c1, c2, t=0.5):
    """Linear blend of two hex colors at fraction t toward c2."""
    import matplotlib.colors as mc
    rgb1 = np.array(mc.to_rgb(c1))
    rgb2 = np.array(mc.to_rgb(c2))
    return mc.to_hex(np.clip(rgb1 * (1 - t) + rgb2 * t, 0, 1))


def main():
    # ── Canvas ───────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(9, 4.5))
    fig.patch.set_facecolor("white")

    # Single axes spanning the whole figure; we draw everything in data coords
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    # ── Hardcoded facts ───────────────────────────────────────────────────────
    N_COLOC      = 368
    N_REPLICATED = 67
    N_TIER1      = 1885
    FOLD_ENRICH  = 1.28
    PVAL_STR     = "7.2×10⁻⁵¹"   # 7.2×10⁻⁵¹  (unicode superscript)
    N_ATLAS_GENES = 27_187
    N_ATLAS_COLS  = 315
    N_COHORTS     = 9
    N_SAMPLES     = 846

    # ── Layout constants ──────────────────────────────────────────────────────
    #   Left bar   : x-center  = 18
    #   Venn left  : x-center  = 40
    #   Venn right : x-center  = 60
    #   Right bar  : x-center  = 82
    BAR_Y_TOP    = 80     # top of bars
    BAR_Y_BOT    = 28     # bottom of bars
    BAR_H        = BAR_Y_TOP - BAR_Y_BOT
    BAR_W        = 10     # half-width of bar → full = 20 data units
    VENN_CY      = 54     # vertical center of Venn circles
    VENN_R       = 13     # radius of each Venn circle
    VENN_OVERLAP = 5.5    # half-distance between circle centers → ~17% area overlap

    LC_X = 18   # left bar x-center
    RC_X = 82   # right bar x-center

    # Left Venn circle center
    VL_X = 44
    # Right Venn circle center
    VR_X = 56

    # ── Helper: filled rounded rect (bar) ────────────────────────────────────
    def _bar(cx, color, alpha=0.88):
        rect = mpatches.FancyBboxPatch(
            (cx - BAR_W / 2, BAR_Y_BOT),
            BAR_W, BAR_H,
            boxstyle="round,pad=0.8",
            facecolor=color, edgecolor=color,
            linewidth=1.2, alpha=alpha, zorder=4,
        )
        ax.add_patch(rect)

    # ── Draw left bar (COLOC / GWAS) ─────────────────────────────────────────
    _bar(LC_X, GWAS_COL)

    # Label inside bar
    ax.text(LC_X, (BAR_Y_TOP + BAR_Y_BOT) / 2 + 4,
            f"{N_COLOC:,}",
            ha="center", va="center", fontsize=18, fontweight="bold",
            color="white", zorder=6)
    ax.text(LC_X, (BAR_Y_TOP + BAR_Y_BOT) / 2 - 6,
            "colocalized\neffector genes",
            ha="center", va="center", fontsize=7.5, color="white",
            linespacing=1.45, zorder=6)

    # Sub-label below bar
    ax.text(LC_X, BAR_Y_BOT - 4.5,
            f"{N_REPLICATED} EUR+EAS\nreplicated",
            ha="center", va="top", fontsize=6.5, color=C_GRAY,
            linespacing=1.45, zorder=5)

    # Modality label above bar
    ax.text(LC_X, BAR_Y_TOP + 4,
            "Genetic effector genes",
            ha="center", va="bottom", fontsize=7.5, fontweight="bold",
            color=GWAS_COL, zorder=5)

    # ── Draw right bar (Bulk DEGs) ────────────────────────────────────────────
    _bar(RC_X, BULK_COL)

    ax.text(RC_X, (BAR_Y_TOP + BAR_Y_BOT) / 2 + 4,
            f"{N_TIER1:,}",
            ha="center", va="center", fontsize=18, fontweight="bold",
            color="white", zorder=6)
    ax.text(RC_X, (BAR_Y_TOP + BAR_Y_BOT) / 2 - 6,
            "disease-\nassociated genes",
            ha="center", va="center", fontsize=7.5, color="white",
            linespacing=1.45, zorder=6)

    ax.text(RC_X, BAR_Y_BOT - 4.5,
            f"{N_ATLAS_GENES:,} genes × {N_ATLAS_COLS} columns\n"
            f"({N_COHORTS} cohorts, n = {N_SAMPLES:,})",
            ha="center", va="top", fontsize=6.5, color=C_GRAY,
            linespacing=1.45, zorder=5)

    ax.text(RC_X, BAR_Y_TOP + 4,
            "Bulk RNA-seq atlas",
            ha="center", va="bottom", fontsize=7.5, fontweight="bold",
            color=BULK_COL, zorder=5)

    # ── Flow arrows: bar → Venn ───────────────────────────────────────────────
    arrow_kw = dict(arrowstyle="-|>", color=C_LGRAY, lw=1.4, mutation_scale=12)

    # Left bar right edge → Venn left circle left edge
    ax.annotate("",
                xy=(VL_X - VENN_R + 0.5, VENN_CY),
                xytext=(LC_X + BAR_W / 2 + 0.8, VENN_CY),
                arrowprops=arrow_kw, zorder=3)

    # Venn right circle right edge → Right bar left edge
    ax.annotate("",
                xy=(RC_X - BAR_W / 2 - 0.8, VENN_CY),
                xytext=(VR_X + VENN_R - 0.5, VENN_CY),
                arrowprops=arrow_kw, zorder=3)

    # ── Venn diagram ──────────────────────────────────────────────────────────
    blend_col = _blend(GWAS_COL, BULK_COL, t=0.5)

    left_circle = mpatches.Circle(
        (VL_X, VENN_CY), VENN_R,
        facecolor=GWAS_COL, edgecolor=GWAS_COL,
        linewidth=0, alpha=0.40, zorder=5,
    )
    right_circle = mpatches.Circle(
        (VR_X, VENN_CY), VENN_R,
        facecolor=BULK_COL, edgecolor=BULK_COL,
        linewidth=0, alpha=0.40, zorder=5,
    )

    # Draw lighter overlap patch: a filled ellipse approximating the lens
    # True lens position: x center = midpoint of the two circle centers = 50
    overlap_cx = (VL_X + VR_X) / 2
    # Approximate overlap region width (chord length) based on geometry
    d = VR_X - VL_X                              # distance between centers = 12
    # Half-chord of intersection at center line
    # For two circles of radius R with center distance d:
    # half_chord = sqrt(R^2 - (d/2)^2)
    half_chord = np.sqrt(max(0, VENN_R**2 - (d / 2)**2))  # ~3.2

    overlap_ellipse = mpatches.Ellipse(
        (overlap_cx, VENN_CY),
        width=d * 0.70,        # horizontal extent of lens approximation
        height=half_chord * 2,
        facecolor=blend_col,
        edgecolor="none",
        alpha=0.55, zorder=6,
    )

    ax.add_patch(left_circle)
    ax.add_patch(right_circle)
    ax.add_patch(overlap_ellipse)

    # Draw thin outlines on top for definition
    for cx_circ, col in [(VL_X, GWAS_COL), (VR_X, BULK_COL)]:
        outline = mpatches.Circle(
            (cx_circ, VENN_CY), VENN_R,
            facecolor="none", edgecolor=col,
            linewidth=1.0, alpha=0.55, zorder=7,
        )
        ax.add_patch(outline)

    # ── Key statistic at overlap center ───────────────────────────────────────
    stat_x = overlap_cx
    stat_y = VENN_CY + 4.5

    # Fold enrichment — large bold
    ax.text(stat_x, stat_y,
            f"{FOLD_ENRICH}×",
            ha="center", va="center",
            fontsize=13, fontweight="bold",
            color=C_TEXT, zorder=10)

    # "enriched for colocalization"
    ax.text(stat_x, stat_y - 5.5,
            "enriched for\ncolocalization",
            ha="center", va="center",
            fontsize=6.8, color=C_TEXT,
            linespacing=1.35, zorder=10)

    # p-value
    ax.text(stat_x, stat_y - 12.0,
            f"continuous Wilcoxon\np = {PVAL_STR}",
            ha="center", va="center",
            fontsize=6.0, color=C_GRAY,
            linespacing=1.35, zorder=10)

    # ── Caveat note below Venn ────────────────────────────────────────────────
    caveat_y = BAR_Y_BOT - 4.5
    ax.text(
        overlap_cx, caveat_y,
        "Modest overlap expected: coding-effect variants\n"
        "(e.g. PNPLA3 p.I148M) excluded by construction",
        ha="center", va="top",
        fontsize=6.0, color=C_GRAY,
        linespacing=1.4, style="italic", zorder=5,
    )

    # ── Section labels above Venn circles ────────────────────────────────────
    ax.text(VL_X - VENN_R / 2, VENN_CY + VENN_R + 3.5,
            "GWAS / eQTL",
            ha="center", va="bottom",
            fontsize=6.5, color=GWAS_COL, fontweight="bold", zorder=5)

    ax.text(VR_X + VENN_R / 2, VENN_CY + VENN_R + 3.5,
            "Bulk RNA-seq",
            ha="center", va="bottom",
            fontsize=6.5, color=BULK_COL, fontweight="bold", zorder=5)

    # ── Save ──────────────────────────────────────────────────────────────────
    misc_dir = os.path.join(_root(), "figures", "misc")
    os.makedirs(misc_dir, exist_ok=True)
    out_path = os.path.join(misc_dir, "fig1f_handoff_seam.pdf")
    fig.savefig(out_path, bbox_inches="tight", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()

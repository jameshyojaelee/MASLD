#!/usr/bin/env python3
"""
Figure 1a — Multi-Modal Data Landscape (Modality Cards)

Publication-quality overview of ALL datasets in the MASLD Transcriptomic Atlas.
Six modality tiles arranged in a 3x2 grid with a dark summary ribbon.
Each tile: modality name, dataset count, scale, species, technology.

Output: PDF (600 DPI) to figures/
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as path_effects
import numpy as np
import os

# ─── Sanjana Lab Publication Standards ───────────────────────────────────────
plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":        10,
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
    "axes.spines.top":  False,
    "axes.spines.right": False,
    "figure.dpi":       150,
})

# ─── Colors ──────────────────────────────────────────────────────────────────
TEXT_DARK   = "#2D3436"
LABEL_GRAY  = "#636E72"
BG_PAGE     = "#FAFAFA"
RIBBON_BG   = "#2D3436"

MODALITY_COLORS = {
    "bulk":    "#2A9D8F",
    "scrna":   "#E76F51",
    "spatial": "#6C5CE7",
    "atac":    "#00B894",
    "gwas":    "#F4A261",
    "pharma":  "#D4A017",
}

SPECIES_COLORS = {
    "Human": ("#E8F5E9", "#2E7D32"),
    "Mouse": ("#E3F2FD", "#1565C0"),
}

# ─── Modality Data ───────────────────────────────────────────────────────────
# (title, color_key, species_list, headline_number, headline_unit,
#  detail_line, tech_line)

MODALITIES = [
    ("Bulk RNA-seq",          "bulk",    ["Human", "Mouse"],
     "17", "cohorts",
     "1,277 human \u00b7 443 mouse samples",
     "Illumina HiSeq / NextSeq / NovaSeq"),

    ("Single-Cell RNA-seq",   "scrna",   ["Human"],
     "7", "datasets",
     "524,701+ cells \u00b7 17 cell types",
     "10x Chromium scRNA / snRNA"),

    ("Spatial Transcriptomics", "spatial", ["Human", "Mouse"],
     "2", "datasets",
     "50 sections \u00b7 6,546+ spots",
     "10x Visium (55 \u00b5m)"),

    ("ATAC-seq",              "atac",    ["Human", "Mouse"],
     "4", "datasets",
     "12 mouse + 73 human samples",
     "Bulk ATAC / snATAC"),

    ("GWAS & eQTL",           "gwas",    ["Human"],
     "12", "resources",
     "778K GWAS \u00b7 1,183 eQTL donors",
     "MR \u00b7 TWAS \u00b7 COLOC \u00b7 sc-eQTL \u00b7 ieQTL"),

    ("Pharmacogenomics",      "pharma",  ["Human"],
     "4", "databases",
     "1,107 LINCS compounds \u00b7 625 PPI drugs",
     "LINCS L1000 \u00b7 CGP \u00b7 STRING \u00b7 SomaScan"),
]

SUMMARY_ITEMS = [
    ("46", "datasets"),
    ("1,912", "bulk samples"),
    ("524K+", "cells"),
    ("778K+", "GWAS subjects"),
    ("7", "evidence layers"),
]


def draw_tile(ax, x, y, w, h, mod_data):
    """Draw a single modality tile at (x, y) with width w, height h."""
    title, color_key, species, num, unit, detail, tech = mod_data
    color = MODALITY_COLORS[color_key]

    # Card background
    card = mpatches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.015",
        facecolor="white",
        edgecolor="#E0E0E0",
        linewidth=0.8,
    )
    ax.add_patch(card)

    # Left color accent bar
    bar_w = w * 0.025
    bar_pad = w * 0.008
    accent = mpatches.FancyBboxPatch(
        (x + bar_pad, y + h * 0.08), bar_w, h * 0.84,
        boxstyle="round,pad=0.003",
        facecolor=color,
        edgecolor="none",
    )
    ax.add_patch(accent)

    # Text x offset (after accent bar)
    tx = x + bar_w + w * 0.06

    # Title
    ax.text(tx, y + h * 0.88, title,
            fontsize=9, fontweight="bold", color=color,
            va="top", ha="left")

    # Species badges (top right)
    badge_x = x + w * 0.97
    for i, sp in enumerate(reversed(species)):
        bg_col, fg_col = SPECIES_COLORS[sp]
        bw, bh = w * 0.16, h * 0.11
        bx = badge_x - (i + 1) * (bw + w * 0.015)
        by = y + h * 0.83
        badge = mpatches.FancyBboxPatch(
            (bx, by), bw, bh,
            boxstyle="round,pad=0.008",
            facecolor=bg_col,
            edgecolor="none",
        )
        ax.add_patch(badge)
        ax.text(bx + bw / 2, by + bh / 2, sp,
                fontsize=5.5, color=fg_col, fontweight="bold",
                ha="center", va="center")

    # Big number + unit
    ax.text(tx, y + h * 0.55, num,
            fontsize=18, fontweight="bold", color=TEXT_DARK,
            va="center", ha="left")
    # Unit next to number
    num_width = len(num) * 0.035 * w + w * 0.02
    ax.text(tx + num_width, y + h * 0.55, unit,
            fontsize=10, color=LABEL_GRAY,
            va="center", ha="left")

    # Detail line
    ax.text(tx, y + h * 0.30, detail,
            fontsize=7, color="#888888",
            va="center", ha="left")

    # Tech line
    ax.text(tx, y + h * 0.13, tech,
            fontsize=6, color="#AAAAAA",
            va="center", ha="left")


def main():
    fig_w, fig_h = 11, 5.5
    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("auto")
    ax.axis("off")
    fig.patch.set_facecolor("white")

    # ── Title ────────────────────────────────────────────────────────────────
    ax.text(0.5, 0.97,
            "Multi-Modal Data Landscape \u2014 MASLD Transcriptomic Atlas",
            fontsize=14, fontweight="bold", color=TEXT_DARK,
            ha="center", va="top",
            transform=ax.transAxes)

    # ── Grid layout (3 cols x 2 rows) ────────────────────────────────────────
    n_cols, n_rows = 3, 2
    margin_x = 0.04
    margin_top = 0.10
    gap_x = 0.025
    gap_y = 0.035
    tile_w = (1 - 2 * margin_x - (n_cols - 1) * gap_x) / n_cols
    tile_h = 0.30

    grid_top = 1.0 - margin_top

    for idx, mod in enumerate(MODALITIES):
        col = idx % n_cols
        row = idx // n_cols
        x = margin_x + col * (tile_w + gap_x)
        y = grid_top - (row + 1) * tile_h - row * gap_y
        draw_tile(ax, x, y, tile_w, tile_h, mod)

    # ── Summary Ribbon ───────────────────────────────────────────────────────
    ribbon_y = 0.05
    ribbon_h = 0.08
    ribbon = mpatches.FancyBboxPatch(
        (margin_x, ribbon_y), 1 - 2 * margin_x, ribbon_h,
        boxstyle="round,pad=0.012",
        facecolor=RIBBON_BG,
        edgecolor="none",
    )
    ax.add_patch(ribbon)

    n_items = len(SUMMARY_ITEMS)
    usable_w = 1 - 2 * margin_x
    item_w = usable_w / n_items

    for i, (val, label) in enumerate(SUMMARY_ITEMS):
        cx = margin_x + (i + 0.5) * item_w
        cy = ribbon_y + ribbon_h / 2

        ax.text(cx, cy + 0.012, val,
                fontsize=12, fontweight="bold", color="white",
                ha="center", va="center")
        ax.text(cx, cy - 0.018, label,
                fontsize=7, color="#BBBBBB",
                ha="center", va="center")

        # Separator line (skip last)
        if i < n_items - 1:
            sep_x = margin_x + (i + 1) * item_w
            ax.plot([sep_x, sep_x],
                    [ribbon_y + ribbon_h * 0.2, ribbon_y + ribbon_h * 0.8],
                    color="#555555", linewidth=0.6)

    # ── Save ─────────────────────────────────────────────────────────────────
    out_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__)))),
        "figures",
    )
    os.makedirs(out_dir, exist_ok=True)
    pdf_path = os.path.join(out_dir, "fig1a_data_landscape.pdf")
    fig.savefig(pdf_path, bbox_inches="tight", dpi=600, facecolor="white")
    plt.close(fig)
    print(f"Saved: {pdf_path}")


if __name__ == "__main__":
    main()

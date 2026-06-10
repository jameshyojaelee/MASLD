#!/usr/bin/env python3
"""
Fig 1D — Mechanism class: coding/protein-function vs regulatory/expression-mediated.

KEY MESSAGE: MASLD's biggest GWAS genes (PNPLA3, TM6SF2…) are genetically real but
invisible to expression methods — because they act through protein function, not
transcription.  The panel shows (LEFT) the GWAS-ATAC filtering cascade and the top
disrupted TFs, and (RIGHT) a two-column class map separating coding/protein-function
genes (COLOC-null) from regulatory/expression-mediated genes (SuSiE PP4 ≥ 0.5).

Data: hard-coded canonical numbers (verified facts pack).
Output: figures/misc/fig1d_mechanism_class.pdf
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np
import os

# ── shared style ────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":         8,
    "pdf.fonttype":      42,
    "ps.fonttype":       42,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.spines.left":  False,
    "axes.spines.bottom": False,
    "figure.dpi":        150,
    "savefig.dpi":       300,
    "figure.facecolor":  "white",
    "axes.facecolor":    "white",
    "axes.grid":         False,
    "legend.frameon":    False,
})

# ── palette ──────────────────────────────────────────────────────────────────
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig1_palette import MODALITY_COLORS, CONTROL_GRAY

GWAS_BLUE = MODALITY_COLORS["gwas"]   # "#000000" — used for funnel fill tinted below
CODING_COLOR  = "#D9D9D9"             # absent / coding-mechanism genes
REG_COLOR     = "#0072B2"             # regulatory / expression-mediated (Okabe-Ito blue)

# Tinted funnel steps: three shades from light→dark of Okabe-Ito blue
FUNNEL_COLORS = ["#B8D6EE", "#5BA3D6", "#1565A8"]

# ── output path ──────────────────────────────────────────────────────────────
def _root():
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

OUT_DIR = os.path.join(_root(), "figures", "misc")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_PDF = os.path.join(OUT_DIR, "fig1d_mechanism_class.pdf")

# ── canonical hard-coded data ─────────────────────────────────────────────────
FUNNEL_LEVELS = [
    (8329, "credible-set variants"),
    (617,  "in open chromatin\n(scATAC peaks)"),
    (434,  "disrupt predicted\nTF motifs"),
]

TOP_TFS = [
    ("HNF4A", 12, "#E69F00"),   # atac color (orange) — master hepatic TF
    ("RORA",  11, "#56B4E9"),   # sky blue
    ("THRB",  10, "#009E73"),   # bluish green
]

# LEFT column: coding/protein-function genes (COLOC-null or marginal)
CODING_GENES = [
    ("PNPLA3",   0.008, ""),
    ("HSD17B13", 0.031, ""),
    ("MBOAT7",   0.044, ""),
    ("TM6SF2",   0.52,  "marginal"),
]

# RIGHT column: regulatory/expression-mediated genes (SuSiE PP4 ≥ 0.5)
REG_GENES = [
    ("RORA",   0.999),
    ("HKDC1",  0.992),
    ("GCKR",   0.950),
    ("MTTP",   0.835),
]

# ── figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(8.5, 5.5))

# Two main axes: left (funnel) and right (class map)
ax_funnel = fig.add_axes([0.04, 0.05, 0.38, 0.88])   # [left, bottom, width, height]
ax_class  = fig.add_axes([0.50, 0.05, 0.48, 0.88])

for ax in (ax_funnel, ax_class):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

# ═══════════════════════════════════════════════════════════════════════════
# LEFT PANEL — vertical variant funnel
# ═══════════════════════════════════════════════════════════════════════════

# Funnel geometry: three rectangles, centred, decreasing width
CENTER = 0.50
WIDTHS  = [0.80, 0.54, 0.42]   # full / medium / narrow
HEIGHTS = [0.14, 0.14, 0.14]
Y_TOPS  = [0.90, 0.68, 0.46]   # y of top edge of each box

# Panel section label
ax_funnel.text(CENTER, 0.97, "GWAS-ATAC\nFiltering Cascade",
               ha="center", va="top", fontsize=9, fontweight="bold", color="#333333")

for i, ((n, label), w, h, y_top, col) in enumerate(
        zip(FUNNEL_LEVELS, WIDTHS, HEIGHTS, Y_TOPS, FUNNEL_COLORS)):

    x0 = CENTER - w / 2
    rect = FancyBboxPatch(
        (x0, y_top - h), w, h,
        boxstyle="round,pad=0.01",
        facecolor=col, edgecolor="white", linewidth=1.2,
        zorder=3
    )
    ax_funnel.add_patch(rect)

    # Bold number
    ax_funnel.text(CENTER, y_top - h / 2 + 0.022, f"{n:,}",
                   ha="center", va="center",
                   fontsize=14, fontweight="bold", color="white", zorder=4)
    # Small label below number
    ax_funnel.text(CENTER, y_top - h / 2 - 0.025, label,
                   ha="center", va="center",
                   fontsize=6.5, color="white", zorder=4,
                   linespacing=1.3)

    # Arrow connector to next level
    if i < len(FUNNEL_LEVELS) - 1:
        arr_x   = CENTER
        arr_y1  = y_top - h - 0.003
        arr_y2  = Y_TOPS[i + 1] + 0.004
        ax_funnel.annotate(
            "", xy=(arr_x, arr_y2), xytext=(arr_x, arr_y1),
            arrowprops=dict(arrowstyle="-|>", color="#555555",
                            lw=1.4, mutation_scale=10),
            zorder=2
        )

# ── TF chips below the last funnel box ───────────────────────────────────────
chip_y      = 0.26
chip_gap    = 0.30
chip_w      = 0.22
chip_h      = 0.07
chip_x_vals = [CENTER - chip_gap, CENTER, CENTER + chip_gap]

ax_funnel.text(CENTER, chip_y + chip_h + 0.05,
               "Top disrupted TFs:",
               ha="center", va="bottom", fontsize=7, color="#555555")

for (tf_name, tf_n, tf_col), cx in zip(TOP_TFS, chip_x_vals):
    chip_rect = FancyBboxPatch(
        (cx - chip_w / 2, chip_y), chip_w, chip_h,
        boxstyle="round,pad=0.015",
        facecolor=tf_col, edgecolor="none", alpha=0.85, zorder=3
    )
    ax_funnel.add_patch(chip_rect)
    ax_funnel.text(cx, chip_y + chip_h / 2 + 0.010, tf_name,
                   ha="center", va="center",
                   fontsize=7, fontweight="bold", color="white", zorder=4)
    ax_funnel.text(cx, chip_y + chip_h / 2 - 0.018, f"n={tf_n}",
                   ha="center", va="center",
                   fontsize=6.2, color="white", zorder=4)

# small note at bottom of funnel
ax_funnel.text(CENTER, 0.04,
               "motif disruptions in\nscATAC-accessible chromatin",
               ha="center", va="bottom", fontsize=6.5, color="#888888",
               linespacing=1.3)

# ═══════════════════════════════════════════════════════════════════════════
# RIGHT PANEL — two-column class map
# ═══════════════════════════════════════════════════════════════════════════

COL_LEFT  = 0.24   # x center of left column
COL_RIGHT = 0.76   # x center of right column
ROW_H     = 0.10   # row height per gene chip
ROW_GAP   = 0.015  # gap between rows
START_Y   = 0.74   # y of top of first gene row

# ── column headers ────────────────────────────────────────────────────────
for cx, title, subtitle, col in [
    (COL_LEFT,  "Coding /",        "protein-function",  "#666666"),
    (COL_RIGHT, "Regulatory /",    "expression-mediated", REG_COLOR),
]:
    ax_class.text(cx, 0.94, title,
                  ha="center", va="bottom", fontsize=9,
                  fontweight="bold", color=col)
    ax_class.text(cx, 0.90, subtitle,
                  ha="center", va="bottom", fontsize=8, color=col)

# dividing vertical line
ax_class.axvline(0.50, ymin=0.05, ymax=0.88,
                 color="#CCCCCC", lw=1.0, zorder=1)

# ── gene row helper ──────────────────────────────────────────────────────
CHIP_W   = 0.34   # width of the name + bar block
CHIP_H   = 0.072
BAR_MAXW = 0.20   # full PP4=1.0 maps to this width
BAR_H    = 0.024
NAME_X_OFFSET = -0.13   # name left edge relative to chip center

def draw_gene_chip(ax, cx, cy, gene, pp4, is_coding, note=""):
    """Draw a rounded chip: gene name on left, PP4 bar on right."""
    bg_col    = CODING_COLOR if is_coding else "#D6E8F5"
    bar_col   = "#AAAAAA"    if is_coding else REG_COLOR
    name_col  = "#555555"    if is_coding else "#1A3D5C"
    half_w    = CHIP_W / 2

    # chip background
    chip = FancyBboxPatch(
        (cx - half_w, cy - CHIP_H / 2), CHIP_W, CHIP_H,
        boxstyle="round,pad=0.012",
        facecolor=bg_col, edgecolor="white", linewidth=0.8, zorder=3
    )
    ax.add_patch(chip)

    # gene name
    ax.text(cx - half_w + 0.02, cy + 0.008, gene,
            ha="left", va="center",
            fontsize=8, fontweight="bold", color=name_col, zorder=4)

    # PP4 bar (right side of chip)
    bar_x0  = cx + 0.01
    bar_y0  = cy - BAR_H / 2 - 0.006
    bar_full_w = half_w - 0.03
    # background track
    track = FancyBboxPatch(
        (bar_x0, bar_y0), bar_full_w, BAR_H,
        boxstyle="round,pad=0.005",
        facecolor="#E8E8E8", edgecolor="none", zorder=3
    )
    ax.add_patch(track)
    # filled bar proportional to pp4
    bar_w = max(bar_full_w * pp4, 0.004)
    filled = FancyBboxPatch(
        (bar_x0, bar_y0), bar_w, BAR_H,
        boxstyle="round,pad=0.005",
        facecolor=bar_col, edgecolor="none", zorder=4
    )
    ax.add_patch(filled)

    # PP4 value
    pp4_label = f"{pp4:.3f}" if pp4 < 0.1 else f"{pp4:.2f}"
    ax.text(bar_x0 + bar_full_w + 0.010, cy - 0.006,
            pp4_label,
            ha="left", va="center", fontsize=6.5, color=name_col, zorder=4)

    # optional marginal note
    if note:
        ax.text(cx - half_w + 0.02, cy - 0.022, note,
                ha="left", va="center",
                fontsize=5.5, color="#999999", style="italic", zorder=4)


# ── "PP4" axis label above bar area ──────────────────────────────────────
for cx in (COL_LEFT, COL_RIGHT):
    ax_class.text(cx + 0.12, START_Y + ROW_H * 0.5 + 0.04,
                  "SuSiE PP4",
                  ha="center", va="bottom", fontsize=6.5, color="#888888")

# ── draw left (coding) genes ──────────────────────────────────────────────
for i, (gene, pp4, note) in enumerate(CODING_GENES):
    cy = START_Y - i * (ROW_H + ROW_GAP)
    draw_gene_chip(ax_class, COL_LEFT, cy, gene, pp4,
                   is_coding=True, note=note)

# ── draw right (regulatory) genes ────────────────────────────────────────
for i, (gene, pp4) in enumerate(REG_GENES):
    cy = START_Y - i * (ROW_H + ROW_GAP)
    draw_gene_chip(ax_class, COL_RIGHT, cy, gene, pp4,
                   is_coding=False)

# ── column footers ────────────────────────────────────────────────────────
last_left_y  = START_Y - (len(CODING_GENES) - 1) * (ROW_H + ROW_GAP) - ROW_H
last_right_y = START_Y - (len(REG_GENES)   - 1) * (ROW_H + ROW_GAP) - ROW_H
footer_y = min(last_left_y, last_right_y) - 0.06

ax_class.text(COL_LEFT, footer_y,
              "COLOC-null:\nprotein mediates risk",
              ha="center", va="top", fontsize=6.5,
              color="#666666", linespacing=1.4,
              style="italic")

ax_class.text(COL_RIGHT, footer_y,
              "n = 368 regulatory-effector genes\n(SuSiE PP4 > 0.5)",
              ha="center", va="top", fontsize=6.5,
              color=REG_COLOR, linespacing=1.4)

# ── x-axis bar legend ─────────────────────────────────────────────────────
legend_y  = 0.065
legend_x0 = 0.10
legend_w  = 0.36

# show a reference scale bar for PP4
ax_class.text(legend_x0, legend_y + 0.030, "PP4 scale:",
              ha="left", va="bottom", fontsize=6, color="#888888")
for val, label in [(0.0, "0"), (0.5, "0.5"), (1.0, "1.0")]:
    tick_x = legend_x0 + val * legend_w
    ax_class.plot([tick_x, tick_x], [legend_y, legend_y + 0.010],
                  color="#AAAAAA", lw=0.8, zorder=3)
    ax_class.text(tick_x, legend_y - 0.004, label,
                  ha="center", va="top", fontsize=5.5, color="#AAAAAA")

ax_class.plot([legend_x0, legend_x0 + legend_w], [legend_y, legend_y],
              color="#AAAAAA", lw=0.8, zorder=3)

# ── save ──────────────────────────────────────────────────────────────────
fig.savefig(OUT_PDF, bbox_inches="tight", dpi=300)
print(f"Saved: {OUT_PDF}")

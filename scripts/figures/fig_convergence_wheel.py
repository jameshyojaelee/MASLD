#!/usr/bin/env python3
"""
Figure 1a (alt) — Multi-Modal Evidence Convergence Wheel

Radial diagram with the MASLD Multi-Evidence Atlas at center.
Six modality sectors arranged as colored arcs, connected to the hub
by subtle lines.  Dataset counts and scale placed alongside each arc.

Usage:
  python fig_convergence_wheel.py                  # standalone with title
  python fig_convergence_wheel.py --no-title       # panel version (no title, saves to fig1b_wheel_notitle.pdf)

Output: PDF (600 DPI) to figures/
"""

import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, Arc, Wedge
import numpy as np
import os

# ─── Sanjana Lab Publication Standards ───────────────────────────────────────
plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":        6,
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
    "figure.dpi":       150,
})

# ─── Colors ──────────────────────────────────────────────────────────────────
TEXT_DARK  = "#2D3436"
LABEL_GRAY = "#636E72"
HUB_COLOR  = "#2D3436"

MODALITY_COLORS = {
    "bulk":    "#2A9D8F",
    "scrna":   "#E76F51",
    "spatial": "#6C5CE7",
    "atac":    "#00B894",
    "gwas":    "#F4A261",
    "pharma":  "#D4A017",
}

# ─── Modality definitions ───────────────────────────────────────────────────
# (key, title, line1, line2, arc_weight)
# arc_weight controls visual thickness (proportional emphasis, not data)

SECTORS = [
    ("bulk",    "Bulk RNA-seq",
     "17 datasets \u00b7 1,912 samples", "Human + Mouse",    5.0),
    ("scrna",   "Single-Cell RNA-seq",
     "8 datasets \u00b7 2.29M+ cells",  "Human + Mouse",    3.5),
    ("spatial", "Spatial Transcriptomics",
     "2 datasets \u00b7 Visium",        "Human + Mouse",    2.0),
    ("atac",    "ATAC-seq",
     "4 datasets \u00b7 85 samples",    "Human + Mouse",    2.0),
    ("pharma",  "Pharmacogenomics",
     "4 datasets \u00b7 1,107 cmpds",  "Cell lines",            3.0),
    ("gwas",    "GWAS & eQTLs",
     "12 GWAS \u00b7 778K+ subj.", "Liver bulk & sc-eQTLs", 3.5),
]


def main(no_title=False):
    fig, ax = plt.subplots(1, 1, figsize=(6.5, 6.5))
    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.35, 1.35)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("white")

    # ── Title removed (publication figure) ─────────────────────────────────
    print("[caption] Multi-omic MASLD Atlas")

    # ── Central hub ──────────────────────────────────────────────────────────
    hub_r = 0.20
    hub = plt.Circle((0, 0), hub_r, facecolor=HUB_COLOR,
                      edgecolor="white", linewidth=1.5, zorder=10)
    ax.add_patch(hub)

    # Hub glow
    for r_off, alpha in [(0.02, 0.04), (0.04, 0.02)]:
        glow = plt.Circle((0, 0), hub_r + r_off,
                           facecolor="none", edgecolor=HUB_COLOR,
                           linewidth=1, alpha=alpha, zorder=9)
        ax.add_patch(glow)

    ax.text(0, 0.025, "MASLD", fontsize=6,
            color="white", ha="center", va="center", zorder=11)
    ax.text(0, -0.07, "Atlas", fontsize=6,
            color="white", ha="center", va="center", zorder=11)

    # ── Sector arcs ──────────────────────────────────────────────────────────
    n = len(SECTORS)
    gap_deg = 8  # gap between arcs in degrees
    total_gap = n * gap_deg
    total_arc = 360 - total_gap

    # Distribute arc span proportional to weight
    total_weight = sum(s[4] for s in SECTORS)
    arc_spans = [(s[4] / total_weight) * total_arc for s in SECTORS]

    arc_radius = 0.62
    arc_lw_base = 14

    angle = 164  # start point customized to center bulk/sc at top and pharma at bottom

    for i, (key, title, line1, line2, weight) in enumerate(SECTORS):
        span = arc_spans[i]
        start_angle = angle
        end_angle = angle - span
        mid_angle = (start_angle + end_angle) / 2
        color = MODALITY_COLORS[key]

        # Draw arc
        lw = arc_lw_base * (weight / max(s[4] for s in SECTORS))
        lw = max(lw, 6)
        arc = Arc((0, 0), 2 * arc_radius, 2 * arc_radius,
                  angle=0, theta1=end_angle, theta2=start_angle,
                  color=color, linewidth=lw, alpha=0.82,
                  capstyle="round", zorder=5)
        ax.add_patch(arc)

        # Connection line from arc midpoint to hub edge
        mid_rad = np.deg2rad(mid_angle)
        arc_x = arc_radius * np.cos(mid_rad)
        arc_y = arc_radius * np.sin(mid_rad)
        hub_x = (hub_r + 0.02) * np.cos(mid_rad)
        hub_y = (hub_r + 0.02) * np.sin(mid_rad)

        ax.plot([hub_x, arc_x], [hub_y, arc_y],
                color=color, linewidth=1.2, alpha=0.25,
                linestyle=(0, (4, 3)), zorder=4)

        # Small dot at arc midpoint
        ax.plot(arc_x, arc_y, 'o', color=color, markersize=4,
                alpha=0.6, zorder=6)

        # Label placement (outside arc)
        label_r = arc_radius + 0.18
        label_x = label_r * np.cos(mid_rad)
        label_y = label_r * np.sin(mid_rad)

        # Determine text alignment based on angle
        cos_val = np.cos(mid_rad)
        if abs(cos_val) < 0.15:
            ha = "center"
        elif cos_val > 0:
            ha = "left"
        else:
            ha = "right"

        ax.text(label_x, label_y + 0.06, title,
                fontsize=6, color=color,
                ha=ha, va="center", zorder=7)
        ax.text(label_x, label_y - 0.04, line1,
                fontsize=6, color="black",
                ha=ha, va="center", zorder=7)
        ax.text(label_x, label_y - 0.13, line2,
                fontsize=6, color="black",
                ha=ha, va="center", zorder=7)

        # Thin line from label to arc
        conn_r = arc_radius + 0.06
        conn_x = conn_r * np.cos(mid_rad)
        conn_y = conn_r * np.sin(mid_rad)
        ax.plot([conn_x, label_x], [conn_y, label_y + 0.00],
                color=color, linewidth=0.7, alpha=0.4, zorder=4)

        # Advance angle
        angle = end_angle - gap_deg

    # ── Summary stats (bottom) ───────────────────────────────────────────────
    summary_items = [
        ("46 datasets", None),
        ("1,912 bulk samples", None),
        ("524K+ cells", None),
        ("778K+ GWAS subjects", None),
    ]
    summary_text = "  \u2502  ".join(s[0] for s in summary_items)
    ax.text(0, -1.25, summary_text,
            fontsize=6, color="black",
            ha="center", va="center", zorder=7)

    # ── Save ─────────────────────────────────────────────────────────────────
    out_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__)))),
        "figures", "fig1", "panels",
    )
    os.makedirs(out_dir, exist_ok=True)

    if no_title:
        # fig1b_wheel_notitle.pdf removed — not used in publication
        plt.close(fig)
        print("--no-title variant disabled; skipping save.")
        return

    pdf_path = os.path.join(out_dir, "fig1a_convergence_wheel.pdf")
    fig.savefig(pdf_path, bbox_inches="tight", dpi=600, facecolor="white")
    plt.close(fig)
    print(f"Saved: {pdf_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-title", action="store_true",
                        help="Omit title (for embedding as a panel in Fig 1)")
    args = parser.parse_args()
    main(no_title=args.no_title)

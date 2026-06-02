#!/usr/bin/env python3
"""
Figure 1a — Study Design Overview (graphic, icon-driven)
Full-width (180mm), matches fig1b sunburst style.
Uses data coordinates (not axes fraction) for proper icon scaling.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, PathPatch, FancyArrowPatch
from matplotlib.path import Path
import matplotlib.colors as mc
import numpy as np
import os

plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":       7,
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
    "figure.dpi":      150,
})

# ── Canonical Okabe-Ito modality palette (single source of truth) ────────────
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig1_palette import MODALITY_COLORS as _MC
MOD_BULK    = _MC["bulk"]
MOD_SCRNA   = _MC["scrna"]
MOD_SPATIAL = _MC["spatial"]
MOD_ATAC    = _MC["atac"]
MOD_PROTEO  = _MC["proteo"]
MOD_GWAS    = _MC["gwas"]
MOD_PHARMA  = _MC["pharma"]

SP_HUMAN = "#6BAED6"
SP_MOUSE = "#74C476"

C_DISC    = "#1565C0"
C_BIO     = "#C2185B"
C_CAUS    = "#7B1FA2"
C_TRANS   = "#00695C"

C_HEALTHY  = "#9E9E9E"
C_MASL     = "#F4A674"
C_MASH     = "#C9265E"
C_FIBROSIS = "#880E4F"

C_TEXT  = "#171A1C"
C_GRAY  = "#636E72"
C_LGRAY = "#BDBDBD"
C_BG    = "white"

FIG_W = 180 / 25.4
FIG_H = 70 / 25.4  # wide + short = landscape banner


def lighten(c, amt=0.85):
    rgb = np.array(mc.to_rgb(c))
    return mc.to_hex(rgb * (1 - amt) + np.ones(3) * amt)


def main():
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    fig.patch.set_facecolor(C_BG)
    ax.set_xlim(0, 100)
    ax.set_ylim(16, 50)
    ax.set_aspect("equal")
    ax.axis("off")

    # =====================================================================
    # LEFT BLOCK: Disease spectrum + species (x: 2–32)
    # =====================================================================

    # ── Liver icons ──
    liver_x = [6, 14, 22, 30]
    liver_y_base = 42
    stages = ["Healthy", "MASL", "MASH", "Fibrosis"]
    stage_cols = [C_HEALTHY, C_MASL, C_MASH, C_FIBROSIS]
    liver_cfg = [
        (lighten(C_HEALTHY, 0.45), C_HEALTHY, 0, 0),
        (lighten(C_MASL, 0.25), C_MASL, 8, 0),
        (lighten(C_MASH, 0.25), C_MASH, 12, 3),
        (lighten(C_FIBROSIS, 0.30), C_FIBROSIS, 6, 8),
    ]

    for i, (stage, col) in enumerate(zip(stages, stage_cols)):
        cx, cy = liver_x[i], liver_y_base
        fc, ec, dots, lines = liver_cfg[i]
        _draw_liver(ax, cx, cy, 2.2, fc, ec, dots, lines)
        ax.text(cx, cy - 2.8, stage, ha="center", va="top",
                fontsize=5.5, fontweight="bold", color=col)
        if i < 3:
            ax.annotate("", xy=(liver_x[i+1]-2.5, cy), xytext=(cx+2.5, cy),
                        arrowprops=dict(arrowstyle="->", color=C_LGRAY, lw=0.7))

    # Gradient bar under livers
    _gradient_bar(ax, liver_x[0]-1.5, liver_x[-1]+1.5, liver_y_base - 3.8,
                  0.35, stage_cols)

    # ── Human sample count ──
    sp_y = 34.5
    _human(ax, 10, sp_y, SP_HUMAN)
    ax.text(14, sp_y + 0.8, "1,259", fontsize=9, fontweight="bold",
            color=SP_HUMAN, va="center")
    ax.text(14, sp_y - 0.6, "human samples (QC-passing)",
            fontsize=4.5, color=C_GRAY, va="center")

    # =====================================================================
    # CENTER BLOCK: Analytical pipeline (x: 36–90)
    # =====================================================================

    modules = [
        ("Genetics",    "Fig 2", C_CAUS,
         ["SuSiE-COLOC · GWAS-ATAC", "Multi-ancestry · sc-eQTL"],
         [MOD_GWAS, MOD_ATAC]),
        ("Discovery",   "Fig 3", C_DISC,
         ["Mega-analysis · LOO-CV", "1,885 DEGs · 27,187 genes"],
         [MOD_BULK]),
        ("Biology",     "Fig 4", C_BIO,
         ["Cascade · Sex", "Deconv. · Cross-species"],
         [MOD_BULK, MOD_SCRNA]),
        ("Translation", "Fig 5", C_TRANS,
         ["Convergence · Drugs", "Proteomics · Spatial"],
         [MOD_PHARMA, MOD_PROTEO, MOD_SPATIAL]),
    ]

    node_xs = [40, 52, 64, 76]
    node_y = 38
    node_r = 4.0

    for i, (title, flbl, color, kws, mod_cols) in enumerate(modules):
        nx = node_xs[i]

        # Outer circle
        circle = plt.Circle((nx, node_y), node_r,
                             facecolor=lighten(color, 0.88),
                             edgecolor=color, linewidth=1.5, zorder=4)
        ax.add_patch(circle)

        # Inner icon  (i=0 Genetics→DNA, i=1 Discovery→volcano, i=2 Biology→sex, i=3 Translation→pill)
        if i == 0:
            _icon_dna(ax, nx, node_y, color)
        elif i == 1:
            _icon_volcano(ax, nx, node_y, color)
        elif i == 2:
            _icon_sex(ax, nx, node_y, color)
        elif i == 3:
            _icon_pill(ax, nx, node_y, color)

        # Modality dots above circle
        _mod_dots(ax, nx, node_y + node_r + 1.0, mod_cols)

        # Title below
        ax.text(nx, node_y - node_r - 0.8, title,
                ha="center", va="top", fontsize=6.5, fontweight="bold", color=color)
        ax.text(nx, node_y - node_r - 2.2, flbl,
                ha="center", va="top", fontsize=5, color=C_GRAY)

        # Keywords (compact, 2 lines max)
        kw_y = node_y - node_r - 3.5
        for ki, kw in enumerate(kws):
            ax.text(nx, kw_y - ki * 1.2, kw,
                    ha="center", va="top", fontsize=3.5, color=C_GRAY)

        # Arrow chevron to next
        if i < 3:
            mid_x = (nx + node_xs[i+1]) / 2
            ax.annotate("", xy=(mid_x + 0.8, node_y),
                        xytext=(mid_x - 0.8, node_y),
                        arrowprops=dict(arrowstyle="-|>", color=C_LGRAY, lw=1.2,
                                        mutation_scale=10))

    # =====================================================================
    # RIGHT BLOCK: Atlas output (x: 80–98)
    # =====================================================================
    # Arrow from pipeline to atlas
    ax.annotate("", xy=(86, 38), xytext=(81, 38),
                arrowprops=dict(arrowstyle="-|>", color=C_LGRAY, lw=1.5,
                                mutation_scale=12))

    # Atlas badge
    atl_x, atl_y = 93, 38
    atl_r = 4.5
    circle = plt.Circle((atl_x, atl_y), atl_r,
                         facecolor="#0D47A1", edgecolor="#0D47A1",
                         linewidth=1.5, zorder=4, alpha=0.95)
    ax.add_patch(circle)
    ax.text(atl_x, atl_y + 1.0, "MASLD", ha="center", va="center",
            fontsize=7, fontweight="bold", color="white", zorder=10)
    ax.text(atl_x, atl_y - 0.5, "Atlas", ha="center", va="center",
            fontsize=7, fontweight="bold", color="white", zorder=10)

    # Modality ring around atlas
    all_mods = [MOD_BULK, MOD_SCRNA, MOD_SPATIAL, MOD_ATAC,
                MOD_PROTEO, MOD_GWAS, MOD_PHARMA]
    for j, mc_ in enumerate(all_mods):
        angle = np.pi/2 + j * 2*np.pi/7
        dx = (atl_r + 1.3) * np.cos(angle)
        dy = (atl_r + 1.3) * np.sin(angle)
        ax.plot(atl_x + dx, atl_y + dy, "o", color=mc_, markersize=5,
                markeredgecolor="white", markeredgewidth=0.4, zorder=10)

    # ── Modality legend (bottom strip — tightened) ──
    mod_names = ["RNA-seq", "scRNA", "Spatial", "ATAC", "Proteo", "Genomics", "Pharma"]
    mod_colors = [MOD_BULK, MOD_SCRNA, MOD_SPATIAL, MOD_ATAC,
                  MOD_PROTEO, MOD_GWAS, MOD_PHARMA]
    leg_y = 22
    leg_x0 = 40
    leg_sp = 8.5
    for j, (nm, cl) in enumerate(zip(mod_names, mod_colors)):
        lx = leg_x0 + j * leg_sp
        ax.plot(lx, leg_y, "s", color=cl, markersize=4.5,
                markeredgecolor="white", markeredgewidth=0.3, zorder=10)
        ax.text(lx + 1.0, leg_y, nm, fontsize=4, color=C_TEXT,
                va="center", fontweight="medium")

    # Species legend (inline, below modalities — human only)
    sp_y = leg_y - 2
    ax.plot(leg_x0, sp_y, "o", color=SP_HUMAN, markersize=4.5)
    ax.text(leg_x0 + 1.0, sp_y, "Human", fontsize=4, color=C_TEXT,
            va="center", fontweight="medium")

    # ── Save ─────────────────────────────────────────────────────────────
    out_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "figures",
    )
    os.makedirs(out_dir, exist_ok=True)
    # Primary output (canonical path)
    pdf_path = os.path.join(out_dir, "fig1a_study_design.pdf")
    fig.savefig(pdf_path, bbox_inches="tight", dpi=600, facecolor=C_BG)
    print(f"Saved: {pdf_path}")
    # Secondary output for genetics-first arc
    misc_dir = os.path.join(out_dir, "misc")
    os.makedirs(misc_dir, exist_ok=True)
    pdf_path2 = os.path.join(misc_dir, "fig1a_genetics_first.pdf")
    fig.savefig(pdf_path2, bbox_inches="tight", dpi=600, facecolor=C_BG)
    print(f"Saved: {pdf_path2}")
    plt.close(fig)


# ═════════════════════════════════════════════════════════════════════════════
# HELPER DRAWING FUNCTIONS (data coordinates, not axes fraction)
# ═════════════════════════════════════════════════════════════════════════════

def _draw_liver(ax, cx, cy, s, fc, ec, dots=0, lines=0):
    verts = [
        (cx-0.6*s, cy+0.05*s), (cx-0.4*s, cy+0.35*s),
        (cx-0.05*s, cy+0.38*s), (cx+0.05*s, cy+0.30*s),
        (cx+0.15*s, cy+0.40*s), (cx+0.55*s, cy+0.35*s),
        (cx+0.65*s, cy+0.05*s), (cx+0.55*s, cy-0.25*s),
        (cx+0.15*s, cy-0.35*s), (cx-0.3*s, cy-0.28*s),
        (cx-0.6*s, cy+0.05*s),
    ]
    codes = [Path.MOVETO] + [Path.CURVE3] * (len(verts) - 1)
    ax.add_patch(PathPatch(Path(verts, codes), fc=fc, ec=ec, lw=0.8,
                           alpha=0.9, zorder=6))
    rng = np.random.RandomState(int(cx*100)+42)
    for _ in range(dots):
        ax.plot(cx+rng.uniform(-0.22*s,0.32*s), cy+rng.uniform(-0.13*s,0.13*s),
                "o", color="white", alpha=0.5, markersize=rng.uniform(1.5,2.5), zorder=7)
    for _ in range(lines):
        lx = cx+rng.uniform(-0.18*s,0.28*s)
        ly = cy+rng.uniform(-0.10*s,0.04*s)
        ax.plot([lx,lx+rng.uniform(-0.03*s,0.03*s)],
                [ly,ly+rng.uniform(0.05*s,0.12*s)],
                color=ec, alpha=0.4, lw=0.5, zorder=7)


def _gradient_bar(ax, x0, x1, y, h, colors):
    n = 60
    w = x1 - x0
    for j in range(n):
        f = j / n
        x = x0 + f * w
        if f < 0.33:
            t = f/0.33
            c = np.array(mc.to_rgb(colors[0]))*(1-t) + np.array(mc.to_rgb(colors[1]))*t
        elif f < 0.66:
            t = (f-0.33)/0.33
            c = np.array(mc.to_rgb(colors[1]))*(1-t) + np.array(mc.to_rgb(colors[2]))*t
        else:
            t = (f-0.66)/0.34
            c = np.array(mc.to_rgb(colors[2]))*(1-t) + np.array(mc.to_rgb(colors[3]))*t
        ax.add_patch(plt.Rectangle((x, y), w/n+0.05, h,
                     facecolor=mc.to_hex(np.clip(c,0,1)), edgecolor="none", zorder=3))


def _human(ax, cx, cy, color):
    ax.plot(cx, cy+1.8, "o", color=color, markersize=6, zorder=8)
    ax.fill([cx-0.8, cx-0.5, cx+0.5, cx+0.8],
            [cy-1.5, cy+0.8, cy+0.8, cy-1.5],
            color=color, alpha=0.8, zorder=7)


def _mouse(ax, cx, cy, color):
    body = plt.Circle((cx, cy), 1.2, fc=color, ec=color, lw=0.5,
                       alpha=0.8, zorder=7)
    ax.add_patch(body)
    for dx in [-0.7, 0.7]:
        ear = plt.Circle((cx+dx, cy+1.0), 0.45, fc=color, ec=color,
                          lw=0.3, alpha=0.7, zorder=7)
        ax.add_patch(ear)
    ax.plot([cx+1.2, cx+2.0, cx+2.5], [cy-0.2, cy-0.7, cy-0.3],
            color=color, lw=1.2, solid_capstyle="round", zorder=7)


def _icon_volcano(ax, cx, cy, color):
    """Mini volcano scatter."""
    rng = np.random.RandomState(77)
    for _ in range(20):
        dx, dy = rng.normal(0, 1.2), rng.exponential(0.6)
        c = color if abs(dx) > 1.0 and dy > 0.8 else C_LGRAY
        ax.plot(cx+dx, cy+dy-1.0, "o", color=c, markersize=1.5,
                alpha=0.7, zorder=8)
    ax.plot([cx-2.5, cx+2.5], [cy, cy], color=C_LGRAY, lw=0.4,
            ls="--", zorder=7)


def _icon_sex(ax, cx, cy, color):
    """Female + Male symbols."""
    # Female
    fc = plt.Circle((cx-1.0, cy+0.5), 0.8, fc="none", ec=color, lw=1.0, zorder=8)
    ax.add_patch(fc)
    ax.plot([cx-1.0, cx-1.0], [cy-0.3, cy-1.5], color=color, lw=0.8, zorder=8)
    ax.plot([cx-1.5, cx-0.5], [cy-1.0, cy-1.0], color=color, lw=0.8, zorder=8)
    # Male
    mc_ = plt.Circle((cx+1.0, cy+0.5), 0.8, fc="none", ec="#1A237E", lw=1.0, zorder=8)
    ax.add_patch(mc_)
    ax.plot([cx+1.6, cx+2.3], [cy+1.1, cy+1.8], color="#1A237E", lw=0.8, zorder=8)
    ax.plot([cx+1.8, cx+2.3, cx+2.3], [cy+1.8, cy+1.8, cy+1.3],
            color="#1A237E", lw=0.7, zorder=8)


def _icon_dna(ax, cx, cy, color):
    """Double helix."""
    t = np.linspace(-2.2, 2.2, 40)
    x1 = cx + 0.8 * np.sin(t * 3)
    x2 = cx - 0.8 * np.sin(t * 3)
    y = cy + t
    ax.plot(x1, y, color=color, lw=1.0, alpha=0.9, zorder=8)
    ax.plot(x2, y, color=MOD_GWAS, lw=1.0, alpha=0.9, zorder=8)
    for k in range(0, len(t), 5):
        ax.plot([x1[k], x2[k]], [y[k], y[k]], color=C_LGRAY, lw=0.3, zorder=7)


def _icon_pill(ax, cx, cy, color):
    """Pill capsule (two half-circles + rectangle)."""
    # Left cap
    angles = np.linspace(np.pi/2, 3*np.pi/2, 30)
    cap_r = 1.2
    xl = cx - 0.8 + cap_r * np.cos(angles)
    yl = cy + cap_r * np.sin(angles)
    ax.fill(xl, yl, color=color, alpha=0.85, zorder=8)
    # Right cap
    angles2 = np.linspace(-np.pi/2, np.pi/2, 30)
    xr = cx + 0.8 + cap_r * np.cos(angles2)
    yr = cy + cap_r * np.sin(angles2)
    ax.fill(xr, yr, color=lighten(color, 0.4), alpha=0.85, zorder=8)
    # Center rect
    ax.fill([cx-0.8, cx+0.8, cx+0.8, cx-0.8],
            [cy-cap_r, cy-cap_r, cy+cap_r, cy+cap_r],
            color=lighten(color, 0.2), alpha=0.85, zorder=7)
    # Divider line
    ax.plot([cx, cx], [cy-cap_r, cy+cap_r], color="white", lw=0.8, zorder=9)


def _mod_dots(ax, cx, y, colors):
    """Row of modality dots centered at cx."""
    n = len(colors)
    spacing = 2.0
    x0 = cx - (n-1)*spacing/2
    for j, c in enumerate(colors):
        ax.plot(x0 + j*spacing, y, "o", color=c, markersize=5.5,
                markeredgecolor="white", markeredgewidth=0.4, zorder=10)


if __name__ == "__main__":
    main()

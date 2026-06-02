#!/usr/bin/env python3
"""
Publication-quality dataset overview figure for the MASLD Cas13 library design paper.

Generates a two-panel horizontal bar chart:
  Left:  Human patient cohorts (6 datasets) — Pink/Purple family
  Right: Mouse diet models (7+ datasets) — Cyan/Green family

Each bar is split into Disease vs Control segments, annotated with total N,
GEO accession, and year. Follows Sanjana Lab Publication Standards.

Output: PDF + PNG (600 DPI) to RNA-seq/figures/
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import numpy as np

# ─── Sanjana Lab Publication Standards ───────────────────────────────────────
plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":        10,
    "axes.titlesize":   13,
    "axes.labelsize":   11,
    "xtick.labelsize":  9,
    "ytick.labelsize":  10,
    "pdf.fonttype":     42,       # Type 42 (TrueType) for editable text
    "ps.fonttype":      42,
    "axes.spines.top":  False,
    "axes.spines.right": False,
    "figure.dpi":       150,
})

# ─── Color Definitions (from KI) ────────────────────────────────────────────
# Human: Pink/Purple family
HUMAN_DISEASE   = "#C9265E"   # Liang deep magenta
HUMAN_CONTROL   = "#9E9E9E"   # Liang neutral gray (project invariant: control = gray)

# Mouse: Cyan/Green family
MOUSE_DISEASE   = "#3a7b36"   # Medium green (species-categorical encoding retained)
MOUSE_CONTROL   = "#9E9E9E"   # Liang neutral gray (project invariant: control = gray)

# Accents
ACCENT_BORDER   = "#333333"
LABEL_GRAY      = "#555555"
BG_HUMAN        = "#fdf2f7"   # Faint pink background
BG_MOUSE        = "#f0f8f4"   # Faint green background

# ─── Dataset Definitions ────────────────────────────────────────────────────
# Each entry: (label, accession, year, n_disease, n_control, note)
human_datasets = [
    ("Chen",     "GSE213621",   2023,  336,  32,  "Duke University"),
    # PRJNA512027 (Gerhard 2018) excluded: L0/S0 library-prep batch
    # confounded with diagnosis. Still used in fibrosis-vs-healthy.
    ("Govaere",  "GSE135251",   2020,  200,  16,  "EU Registry"),
    ("Kozumi",   "GSE167523",   2021,   98,   0,  "NAFL vs NASH"),
    ("Hoang",    "GSE130970",   2019,   52,  26,  "NAS graded"),
    ("Suppli",   "GSE126848",   2019,   28,  29,  "BMI spectrum"),
]

mouse_datasets = [
    ("LIDPAD",       "GSE159911",   2020, 80,  71, "16 timepoints"),
    ("FPC/CDAHFD",   "GSE162876",   2020, 140, 73, "2 diets"),
    ("HFD (Long)",   "GSE274914",   2024, 11,  11, "7w/52w"),
    ("HFD (Short)",  "GSE224069",   2023, 28,  15, "19w"),
    ("GAN",          "GSE225616",   2023, 19,   0, "No ctrl, vehicle"),
    ("Paquette MCD", "GSE156918",   2021, 11,  11, "Cre/FLCN KO"),
    ("Yue MCD",      "GSE205974",   2022,  3,   3, "Validation"),
    ("In-House MCD", "Sanjana Lab", 2024,  6,   6, "3 timepoints"),
]

# Sort by total N (largest on top)
human_datasets.sort(key=lambda x: x[3] + x[4])
mouse_datasets.sort(key=lambda x: x[3] + x[4])


def draw_panel(ax, datasets, disease_color, control_color, species_label, bg_color):
    """Draw a horizontal stacked bar chart for one species panel."""
    
    n = len(datasets)
    y_pos = np.arange(n)
    bar_height = 0.65
    
    # Faint background
    ax.set_facecolor(bg_color)
    
    disease_vals = [d[3] for d in datasets]
    control_vals = [d[4] for d in datasets]
    totals = [d[3] + d[4] for d in datasets]
    max_total = max(totals)
    
    # Draw control bars (left part)
    ctrl_bars = ax.barh(y_pos, control_vals, height=bar_height,
                        color=control_color, edgecolor="white", linewidth=0.5,
                        zorder=3, label="Control")
    
    # Draw disease bars (stacked to the right of control)
    dis_bars = ax.barh(y_pos, disease_vals, left=control_vals, height=bar_height,
                       color=disease_color, edgecolor="white", linewidth=0.5,
                       zorder=3, label="Disease")
    
    # Y-axis labels: Author (Year)
    y_labels = []
    for d in datasets:
        label = f"{d[0]} ({d[2]})"
        y_labels.append(label)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels, fontweight="bold", fontsize=10)
    
    # Annotate: total N + accession on each bar
    for i, d in enumerate(datasets):
        total = d[3] + d[4]
        acc = d[1]
        
        # Total N inside/outside bar
        x_text = total + max_total * 0.02
        ax.text(x_text, y_pos[i] + 0.01, f"n={total}",
                va="center", ha="left", fontsize=8.5, fontweight="bold",
                color=ACCENT_BORDER, zorder=5)
        
        # GEO accession below label
        ax.text(x_text, y_pos[i] - 0.23, acc,
                va="center", ha="left", fontsize=7,
                color=LABEL_GRAY, style="italic", zorder=5)
    
    # Panel title
    ax.set_title(species_label, fontsize=14, fontweight="bold", pad=12,
                 color=disease_color)
    
    # X-axis
    ax.set_xlabel("Number of Samples", fontsize=10, color=LABEL_GRAY)
    ax.set_xlim(0, max_total * 1.45)
    
    # Liang style: no dashed gridlines on stacked bars
    ax.set_axisbelow(True)
    
    # Spine styling
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.spines["left"].set_color("#999999")
    ax.spines["bottom"].set_color("#999999")
    ax.tick_params(axis="y", length=0)
    
    # Invert so largest is on top
    ax.invert_yaxis()


def main():
    # ── Create Figure ────────────────────────────────────────────────────
    fig, (ax_human, ax_mouse) = plt.subplots(
        1, 2, figsize=(14, 5.5),
        gridspec_kw={"wspace": 0.45}
    )
    
    # Draw panels
    draw_panel(ax_human, human_datasets,
               disease_color=HUMAN_DISEASE,
               control_color=HUMAN_CONTROL,
               species_label="Human Patient Cohorts",
               bg_color=BG_HUMAN)
    
    draw_panel(ax_mouse, mouse_datasets,
               disease_color=MOUSE_DISEASE,
               control_color=MOUSE_CONTROL,
               species_label="Mouse Diet Models",
               bg_color=BG_MOUSE)
    
    # ── Legends ──────────────────────────────────────────────────────────
    # Human legend
    h_patches = [
        mpatches.Patch(facecolor=HUMAN_DISEASE, edgecolor="white", label="Disease"),
        mpatches.Patch(facecolor=HUMAN_CONTROL, edgecolor="white", label="Control"),
    ]
    ax_human.legend(handles=h_patches, loc="center right", framealpha=0.9,
                    fontsize=8, edgecolor="#dddddd", title="Condition",
                    title_fontsize=9)
    
    # Mouse legend
    m_patches = [
        mpatches.Patch(facecolor=MOUSE_DISEASE, edgecolor="white", label="Disease/Diet"),
        mpatches.Patch(facecolor=MOUSE_CONTROL, edgecolor="white", label="Control"),
    ]
    ax_mouse.legend(handles=m_patches, loc="center right", framealpha=0.9,
                    fontsize=8, edgecolor="#dddddd", title="Condition",
                    title_fontsize=9)
    
    # ── Suptitle ─────────────────────────────────────────────────────────
    fig.suptitle("RNA-seq Dataset Overview — MASLD Transcriptomic Atlas",
                 fontsize=15, fontweight="bold", y=0.98, color="#333333")
    
    # ── Summary annotations ──────────────────────────────────────────────
    human_total = sum(d[3] + d[4] for d in human_datasets)
    mouse_total = sum(d[3] + d[4] for d in mouse_datasets)
    
    ax_human.text(0.97, 0.97,
                  f"Total: {human_total} samples\n6 cohorts",
                  transform=ax_human.transAxes, ha="right", va="top",
                  fontsize=9, color=HUMAN_DISEASE, fontweight="bold",
                  bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                            edgecolor=HUMAN_DISEASE, alpha=0.85))
    
    ax_mouse.text(0.97, 0.97,
                  f"Total: {mouse_total} samples\n8 datasets, 5 diets",
                  transform=ax_mouse.transAxes, ha="right", va="top",
                  fontsize=9, color=MOUSE_DISEASE, fontweight="bold",
                  bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                            edgecolor=MOUSE_DISEASE, alpha=0.85))
    
    # ── Save ─────────────────────────────────────────────────────────────
    import os
    out_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), "figures")
    os.makedirs(out_dir, exist_ok=True)
    
    pdf_path = os.path.join(out_dir, "dataset_overview.pdf")

    fig.savefig(pdf_path, bbox_inches="tight", dpi=600)
    plt.close(fig)

    print(f"✅ Saved: {pdf_path}")


if __name__ == "__main__":
    main()

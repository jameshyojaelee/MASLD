#!/usr/bin/env python3
"""
Figure 1b (alt) — Multi-Modal Evidence Convergence Wheel, Modular Sunburst Variation
Matches the exact style of the screenshot with proper layout, right-hand legend,
radiating outer bars, and layered structure.
"""

import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import numpy as np
import matplotlib.colors as mc
import os

# ─── Sanjana Lab Publication Standards ───────────────────────────────────────
plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":        10,
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
    "figure.dpi":       150,
})

TEXT_DARK  = "#171A1C" # slightly deeper dark
LABEL_GRAY = "#636E72"
BG_COLOR   = "white"

# Enhanced Palette for deeper contrast
MODALITY_COLORS = {
    "bulk":    "#218C7E",  # deeper teal
    "scrna":   "#D95B3D",  # deeper terracotta
    "spatial": "#5D4FCF",  # deep indigo 
    "atac":    "#00A381",  # deep mint
    "gwas":    "#E7863A",  # deep warm orange 
    "pharma":  "#C79313",  # rich gold
}

# (key, title, l1, l2, l3, weight, datasets)
SECTORS = [
    ("bulk",    "Bulk\nRNA-seq",           "17\ndatasets",   "1,912\nsamples", "Human\n+ Mouse",  17.0, 17),
    ("scrna",   "Single-Cell\nRNA-seq",    "8\ndatasets",   "2.29M+\ncells",  "Human\n+ Mouse",  8.0, 8),
    ("spatial", "Spatial\nOmics",          "2\ndatasets",   "40\nsamples",    "Human\n+ Mouse",  2.0, 2),
    ("atac",    "ATAC-seq",                "4\ndatasets",   "85\nsamples",    "Human\n+ Mouse",  4.0, 4),
    ("pharma",  "Pharmaco-\ngenomics",     "4\ndatasets",  "1,107\ncmpds",   "Cell\nlines",     4.0, 4),
    ("gwas",    "GWAS &\neQTLs",           "12\nGWAS",      "778K+\nsubj.",   "Liver &\nsc-eQTLs", 12.0, 12),
]

def get_text_color(bg_color):
    rgb = mc.to_rgb(bg_color)
    lum = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
    return "white" if lum < 0.65 else TEXT_DARK

def blend_white(color, factor):
    rgb = np.array(mc.to_rgb(color))
    white = np.array([1.0, 1.0, 1.0])
    blended = rgb * (1 - factor) + white * factor
    return mc.to_hex(blended)

def main(no_title=False):
    fig = plt.figure(figsize=(10.5, 7.5))
    fig.patch.set_facecolor(BG_COLOR)
    
    # Left axes: Main sunburst
    ax = fig.add_axes([0.02, 0.05, 0.65, 0.90])
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.55, 1.55)
    ax.set_aspect("equal")
    ax.axis("off")
    
    # Right axes: Half-sunburst legend
    ax_leg = fig.add_axes([0.72, 0.40, 0.25, 0.40])
    ax_leg.set_xlim(-0.3, 1.3)
    ax_leg.set_ylim(-0.1, 1.5)
    ax_leg.set_aspect("equal")
    ax_leg.axis("off")

    # ── Title ────────────────────────────────────────────────────────────────
    if not no_title:
        fig.text(0.5, 0.95, "Multi-omic MASLD Atlas - Sunburst Layout",
                 fontsize=15, fontweight="bold", color=TEXT_DARK,
                 ha="center", va="top")

    # =========================================================================
    # MAIN SUNBURST 
    # =========================================================================

    # Central hub
    hub_r = 0.28
    hub = plt.Circle((0, 0), hub_r, facecolor="white",
                      edgecolor=TEXT_DARK, linewidth=1.5, zorder=10)
    ax.add_patch(hub)
    
    # Decorator ring for hub
    hub_outer = plt.Circle((0, 0), hub_r + 0.015, facecolor="none",
                           edgecolor=LABEL_GRAY, linewidth=0.8, zorder=9, alpha=0.5)
    ax.add_patch(hub_outer)

    ax.text(0, 0, "MASLD\nAtlas", fontsize=13, fontweight="bold",
            color=TEXT_DARK, ha="center", va="center", zorder=11)

    # Rings (Hierarchical thicknesses)
    r_starts = [0.31, 0.58, 0.78]
    r_ends   = [0.58, 0.78, 1.05]
    
    # Reverse shading so inner rings are pale and outer rings are fully vibrant
    layer_shading = [0.45, 0.20, 0.0]
    
    total_weight = sum(s[5] for s in SECTORS)
    gap_deg = 3.5  # Distinct sector gaps
    total_arc = 360 - len(SECTORS) * gap_deg
    
    current_angle = 164.0
    rng = np.random.RandomState(42)  # For outer bars
    
    for i, (key, title, l1, l2, l3, weight, n_data) in enumerate(SECTORS):
        span = (weight / total_weight) * total_arc
        end_angle = current_angle - span
        
        base_color = MODALITY_COLORS[key]
        layer_texts = [title, l1, l2]  # Dropped l3
        
        # 1. 3 Concentric Layers (Inverted Text Order)
        for layer_idx in range(3):
            r_in = r_starts[layer_idx]
            r_out = r_ends[layer_idx]
            
            # Read texts in reverse order: inner ring gets index 2, outer gets index 0
            text_idx = 2 - layer_idx
            text_str = layer_texts[text_idx]
            
            color = blend_white(base_color, layer_shading[layer_idx])
            
            # The Wedge
            wedge = Wedge((0, 0), r_out, end_angle, current_angle, width=r_out - r_in,
                          facecolor=color, edgecolor=BG_COLOR, linewidth=2.0, zorder=5)
            ax.add_patch(wedge)
            
            # Text placement
            r_mid = (r_in + r_out) / 2.0
            mid_angle = (current_angle + end_angle) / 2.0
            
            # Text rotation logic
            rot = mid_angle - 90
            if rot < -90:
                rot += 180
            elif rot > 90:
                rot -= 180
                
            x = r_mid * np.cos(np.radians(mid_angle))
            y = r_mid * np.sin(np.radians(mid_angle))
            
            # Modality is always vibrant base_color -> force white text
            if text_idx == 0:
                t_color = "white"
            else:
                t_color = get_text_color(color)
            
            # Base sizes based on ring
            base_size = 9.0 if text_idx == 0 else (8.0 if text_idx == 1 else 7.0)
            
            # Dynamic scaling for tight slices (span < 15 deg)
            if span < 15:
                scale_factor = max(0.5, span / 15.0)
                f_size = base_size * scale_factor
            else:
                f_size = base_size
                
            f_weight = "bold" if text_idx == 0 else "medium"
                
            ax.text(x, y, text_str, rotation=rot, color=t_color,
                    fontsize=f_size, fontweight=f_weight,
                    ha="center", va="center", zorder=6)
                    
        # 2. Outer bars (Layer 4) representing # datasets
        margin = 1.0
        usable_span = span - 2 * margin
        if n_data > 0:
            bar_span = usable_span / n_data
            for j in range(n_data):
                b_start = current_angle - margin - j * bar_span
                b_end = b_start - bar_span
                
                # Dynamic width and length for spikes
                b_gap = max(0.6, bar_span * 0.18)
                length = 0.15 + rng.uniform(0.0, 0.40)
                
                # Small spacer between ring and spikes
                r_base = r_ends[-1] + 0.02 
                spike_color = base_color
                
                w_bar = Wedge((0, 0), r_base + length, b_end + b_gap, b_start - b_gap, 
                              width=length, facecolor=spike_color, edgecolor=BG_COLOR, 
                              linewidth=1.2, zorder=4)
                ax.add_patch(w_bar)
            
        current_angle = end_angle - gap_deg

    # =========================================================================
    # HALF-SUNBURST LEGEND
    # =========================================================================
    
    cx, cy = 0.0, 0.0
    
    # 'rainbow' legend look
    leg_colors = [MODALITY_COLORS["scrna"], MODALITY_COLORS["bulk"]]
    leg_angles = [(45, 90), (0, 45)]
    labels = ["Counts", "Datasets", "Modality"]  # Updated Order
    
    # Draw Quarter Rings
    for idx in range(3):
        r_in = r_starts[idx]
        r_out = r_ends[idx]
        
        # Two sectors for the quarter ring
        for ang_idx, (th1, th2) in enumerate(leg_angles):
            c = blend_white(leg_colors[ang_idx], layer_shading[idx])
            w = Wedge((cx, cy), r_out, th1, th2, width=r_out - r_in,
                      facecolor=c, edgecolor=BG_COLOR, lw=2.0)
            ax_leg.add_patch(w)
            
        # Add Label to the left
        r_mid = (r_in + r_out) / 2.0
        ax_leg.text(-0.06, r_mid, labels[idx], ha="right", va="center", 
                    fontsize=11, color=TEXT_DARK, fontweight="medium")
        # Add connecting tick
        ax_leg.plot([-0.04, -0.01], [r_mid, r_mid], color="black", lw=1.0)
        
    # Outer bars legend
    ax_leg.text(-0.06, r_ends[-1] + 0.2, "n of datasets", ha="right", va="center", 
                fontsize=11, color=TEXT_DARK, fontweight="medium")
    ax_leg.plot([-0.04, -0.01], [r_ends[-1] + 0.2, r_ends[-1] + 0.2], color="black", lw=1.0)
    
    # Draw a few dummy outer bars at 60 and 30 deg
    for ang, c_idx in zip([30, 60], [1, 0]):
        c = leg_colors[c_idx]
        for offset in [-6, 0, 6]:
            length = 0.18 + rng.uniform(0.0, 0.20)
            w = Wedge((cx, cy), r_ends[-1] + length + 0.02, ang+offset-1.8, ang+offset+1.8, 
                      width=length, facecolor=c, edgecolor=BG_COLOR, lw=1.2)
            ax_leg.add_patch(w)

    # ── Save ─────────────────────────────────────────────────────────────────
    out_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__)))),
        "figures",
    )
    os.makedirs(out_dir, exist_ok=True)

    if no_title:
        pdf_path = os.path.join(out_dir, "fig1b_wheel_sunburst_notitle.pdf")
    else:
        pdf_path = os.path.join(out_dir, "fig1b_wheel_sunburst.pdf")

    fig.savefig(pdf_path, bbox_inches="tight", dpi=600, facecolor=BG_COLOR)
    plt.close(fig)
    print(f"Saved: {pdf_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-title", action="store_true",
                        help="Omit title")
    args = parser.parse_args()
    main(no_title=args.no_title)

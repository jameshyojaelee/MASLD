#!/usr/bin/env python3
"""
Figure 1b (alt) — Multi-Modal Evidence Convergence Wheel, Modular Sunburst Variation
Matches the exact style of the screenshot with proper layout, right-hand legend,
radiating outer bars, and layered structure.
"""

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

TEXT_DARK  = "#2B2B2B" # warm near-black (Liang-style ink, not pure jet)
LABEL_GRAY = "#7A7A7A"
BG_COLOR   = "white"

# Liang et al. 2025 cell-line palette extended to 7 modalities. Each hue
# sits in its own family with ~50% chroma — distinct at panel scale but
# softer than primary jewel tones.
# Canonical Okabe-Ito modality palette (single source of truth)
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from fig1_palette import MODALITY_COLORS

# Species ring: matched to the modality saturation so it sits in the same
# visual register without dominating.
BINARY_COLORS = {
    "Human": "#6FA2C8",  # cornflower
    "Mouse": "#8FBE6F",  # spring sage
}

# (key, title, l1, l2, (n_human, n_mouse), weight, datasets)
# Verified from unified metadata & processed data (2026-03-25 audit)
# scRNA: 283 human donors + 114 mouse donors across 7 source datasets
# GWAS: 13 GWAS (Ghodsian, Chen, UKBB ALT/AST/GGT, PDFF, FinnGen NAFLD/NASH/HCC, BBJ ALT/AST/GGT)
#       + 3 eQTL (GTEx v8, Broadaway, sc-eQTL); ~1.5M+ unique subjects across ancestries
SECTORS = [
    ("bulk",    "RNA-seq",                 "18\ndatasets",    "1,741\nsamples",  (1277, 464),  18.0, 18),
    ("scrna",   "scRNA-seq",               "7\ndatasets",    "2.29M+\ncells",   (283, 114),    8.0, 7),
    ("spatial", "Spatial Omics",           "3\ndatasets",    "60\nsamples",     (60, 0),        4.5, 3),
    ("atac",    "ATAC-seq",                "2\ndatasets",    "30\nsamples",     (18, 12),       4.5, 2),
    ("proteo",  "Proteomics",              "3\ndatasets",    "307\nsamples",    (307, 0),       4.5, 3),
    ("pharma",  "Pharmaco-genomics",       "4\ndatasets",   "1,107\ncmpds",    (1, 0),         4.0, 4),
    ("gwas",    "Genomics",                "13 GWAS\n3 eQTL","1.5M+\nsubj.",  (1, 0),        16.0, 16),
]

# Per-dataset sample counts for proportional outer-bar heights
# Verified from unified metadata & processed data (2026-03-25 audit)
_DATASET_SIZES_RAW = {
    "bulk": [  # 9 Human cohorts + 9 Mouse datasets (sample N from unified_metadata.csv;
               # PRJNA512027 = 185 samples dropped from cohort presentation for
               # L0/S0 library-prep batch confound with diagnosis)
        367, 216, 164, 143, 98, 94, 78, 67, 57,           # Human (GSE213621→GSE126848)
        213, 151, 29, 22, 12, 11, 10, 10, 6,               # Mouse (GSE162876→GSE205974)
    ],
    "scrna": [  # 7 source datasets (donors/samples from sample_manifest.csv)
        158, 117, 67, 21, 20, 10, 4,  # Liver_Atlas, GSE244832, GSE202379, GSE185477, GSE136103, GSE189600, GSE174748
    ],
    "spatial": [35, 15, 10],            # HRA007511, GSE192741, Vu_et_al_2025 (sections/samples)
    "atac": [18, 12],                   # Human_Multiome (18 donors), Mouse_Bulk (12 samples)
    "proteo": [177, 72, 58],            # GSE276114 SomaScan, PXD052937 DIA-MS plasma, PXD051911 DIA-MS liver
    "pharma": [1107, 625, 58, 20],      # LINCS compounds, network drugs, multi-layer, controls
    "gwas": [  # 13 GWAS + 3 eQTL (donors/subjects, in thousands)
        778, 691, 377, 377, 377, 344, 350, 340,   # Ghodsian, Chen, 3×FinnGen, UKBB ALT/AST/GGT
        261, 261, 164, 85, 36,                     # BBJ ALT/AST/GGT, other, PDFF
        1.183, 0.312, 0.208,                       # Broadaway, sc-eQTL, GTEx (thousands)
    ],
}

# Shuffle within each sector (seeded) so bars aren't monotonically sorted
_rng_shuffle = np.random.RandomState(99)
DATASET_SIZES = {}
for _k, _v in _DATASET_SIZES_RAW.items():
    _arr = list(_v)
    _rng_shuffle.shuffle(_arr)
    DATASET_SIZES[_k] = _arr

BAR_LEN_MIN = 0.08   # shortest bar (smallest dataset)
BAR_LEN_MAX = 0.45   # tallest bar  (largest dataset)

# Per-sector multiplier so stored values map to real-world units for labels
# GWAS values are stored in thousands → multiply by 1000 for display
_UNIT_MULT = {"gwas": 1000}

# Per-sector units for axis labels
_SECTOR_UNITS = {
    "bulk":    "samples",
    "scrna":   "donors",
    "spatial": "sections",
    "atac":    "samples",
    "proteo":  "samples",
    "pharma":  "cmpds",
    "gwas":    "subjects",
}

def _fmt_axis_val(key, raw_val):
    """Format a value for a radial axis tick, respecting unit multipliers."""
    real = raw_val * _UNIT_MULT.get(key, 1)
    if real >= 1_000_000:
        return f"{real/1e6:.1f}M"
    elif real >= 10_000:
        return f"{real/1e3:.0f}K"
    elif real >= 1_000:
        return f"{real/1e3:.1f}K"
    elif real < 1:
        return f"{int(real*1000)}"
    else:
        return f"{int(real)}"

def get_text_color(bg_color):
    rgb = mc.to_rgb(bg_color)
    lum = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
    return "white" if lum < 0.65 else TEXT_DARK

def blend_white(color, factor):
    rgb = np.array(mc.to_rgb(color))
    white = np.array([1.0, 1.0, 1.0])
    blended = rgb * (1 - factor) + white * factor
    return mc.to_hex(blended)

def main():
    fig = plt.figure(figsize=(10.5, 7.5))
    fig.patch.set_facecolor(BG_COLOR)

    # Left axes: Main sunburst
    ax = fig.add_axes([0.02, 0.05, 0.65, 0.90])
    ax.set_xlim(-1.65, 1.65)
    ax.set_ylim(-1.65, 1.65)
    ax.set_aspect("equal")
    ax.axis("off")

    # Right axes: Half-sunburst legend (shifted down to avoid modality legend overlap)
    ax_leg = fig.add_axes([0.72, 0.27, 0.25, 0.40])
    ax_leg.set_xlim(-0.3, 1.3)
    ax_leg.set_ylim(-0.8, 1.5)
    ax_leg.set_aspect("equal")
    ax_leg.axis("off")

    # Top-right axes: Modality legend (taller to fit 7 modalities)
    ax_mod = fig.add_axes([0.75, 0.68, 0.2, 0.28])
    ax_mod.set_xlim(-0.05, 1.2)
    ax_mod.set_ylim(-0.15, 1.25)
    ax_mod.axis("off")

    # =========================================================================
    # MODALITY LEGEND
    # =========================================================================
    mod_y = 1.0
    ax_mod.text(0, mod_y + 0.13, "Data Modality", fontsize=11, fontweight="bold", color=TEXT_DARK)
    for key, color in MODALITY_COLORS.items():
        title = next(s[1] for s in SECTORS if s[0] == key).replace("\n", " ")
        patch = plt.Rectangle((0, mod_y - 0.04), 0.09, 0.09, facecolor=color, edgecolor=BG_COLOR)
        ax_mod.add_patch(patch)
        ax_mod.text(0.14, mod_y, title, fontsize=9.5, color=TEXT_DARK, va="center")
        mod_y -= 0.155

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
    # Layer 0: Counts (Inner)
    # Layer 1: Datasets 
    # Layer 2: Binary Species flag (Outermost slim ring)
    r_starts = [0.31, 0.65, 0.95]
    r_ends   = [0.65, 0.95, 1.05]
    
    # Shading: inner ring sits ~65% toward white (lighter wash for the
    # count layer), outer modality ring softened ~15% so it doesn't punch
    # against the muted palette.
    layer_shading = [0.65, 0.15]  # Layer 2 (species) handled separately
    
    total_weight = sum(s[5] for s in SECTORS)
    gap_deg = 3.5  # Distinct sector gaps
    total_arc = 360 - len(SECTORS) * gap_deg
    
    current_angle = 185.0
    rng = np.random.RandomState(42)

    
    for i, (key, title, l1, l2, context, weight, n_data) in enumerate(SECTORS):
        span = (weight / total_weight) * total_arc
        end_angle = current_angle - span
        
        base_color = MODALITY_COLORS[key]
        layer_texts = [l2, l1]  # Inner to Outer texts
        
        # 1. First 2 Concentric Layers (Modality Colors)
        for layer_idx in range(2):
            r_in = r_starts[layer_idx]
            r_out = r_ends[layer_idx]
            
            text_str = layer_texts[layer_idx]
            color = blend_white(base_color, layer_shading[layer_idx])
            
            # The Wedge
            wedge = Wedge((0, 0), r_out, end_angle, current_angle, width=r_out - r_in,
                          facecolor=color, edgecolor=BG_COLOR, linewidth=2.0, zorder=5)
            ax.add_patch(wedge)
            
            # Text placement — adaptive for narrow wedges
            if text_str:
                r_mid = (r_in + r_out) / 2.0
                mid_angle = (current_angle + end_angle) / 2.0

                rot = mid_angle - 90
                if rot < -90:
                    rot += 180
                elif rot > 90:
                    rot -= 180

                x = r_mid * np.cos(np.radians(mid_angle))
                y = r_mid * np.sin(np.radians(mid_angle))

                t_color = get_text_color(color)
                base_size = 9.0 if layer_idx == 0 else 8.0

                if span < 18:
                    f_size = base_size * max(0.55, span / 22.0)
                elif span < 25:
                    f_size = base_size * max(0.7, span / 25.0)
                else:
                    f_size = base_size

                ax.text(x, y, text_str, rotation=rot, color=t_color,
                        fontsize=f_size, fontweight="medium",
                        ha="center", va="center", zorder=6)
        
        # 2. 3rd Layer: Binary Species Flag (Proportionally Scaled)
        total_samples = sum(context)
        if total_samples > 0:
            h_frac = context[0] / total_samples
            m_frac = context[1] / total_samples
            
            # Start from outer bound (current_angle)
            seg_start = current_angle
            
            for f, col_name in zip([h_frac, m_frac], ["Human", "Mouse"]):
                if f > 0:
                    seg_span = f * span
                    seg_end = seg_start - seg_span
                    
                    sp_color = BINARY_COLORS[col_name]
                    w_species = Wedge((0, 0), r_ends[2], seg_end, seg_start, width=r_ends[2] - r_starts[2],
                                      facecolor=sp_color, edgecolor=BG_COLOR, linewidth=1.5, zorder=5)
                    ax.add_patch(w_species)
                    
                    seg_start = seg_end
                          
        # 3. Outer bars (Layer 4) — one bar per dataset, height ∝ dataset size (linear % of max)
        margin = 1.0
        usable_span = span - 2 * margin
        r_base = r_ends[-1] + 0.02
        if n_data > 0:
            sizes = DATASET_SIZES.get(key, [1] * n_data)
            s_max = max(sizes)
            bar_span = usable_span / n_data
            for j in range(n_data):
                b_start = current_angle - margin - j * bar_span
                b_end = b_start - bar_span

                b_gap = max(0.6, bar_span * 0.18)
                # Linear scaling: bar height = fraction of sector max
                frac = sizes[j] / s_max if s_max > 0 else 0.5
                length = BAR_LEN_MIN + frac * (BAR_LEN_MAX - BAR_LEN_MIN)

                spike_color = "#B0B0B0"  # uniform muted gray for all bars

                w_bar = Wedge((0, 0), r_base + length, b_end + b_gap, b_start - b_gap,
                              width=length, facecolor=spike_color, edgecolor=BG_COLOR,
                              linewidth=0.8, alpha=0.55, zorder=4)
                ax.add_patch(w_bar)

            # Per-sector radial axis: line at the leading edge of the sector
            # with a single tick at the per-modality max value, so bar heights
            # are interpretable in the modality's own units (samples / donors /
            # subjects / compounds).
            axis_ang = current_angle - margin
            axis_rad = np.radians(axis_ang)
            r_axis0  = r_base
            r_axis1  = r_base + BAR_LEN_MAX
            # Axis stem
            ax.plot([r_axis0 * np.cos(axis_rad), r_axis1 * np.cos(axis_rad)],
                    [r_axis0 * np.sin(axis_rad), r_axis1 * np.sin(axis_rad)],
                    color="#555555", linewidth=0.7, zorder=6)
            # Outward tick at max, perpendicular to the radial axis
            tick_len = 0.025
            tang_x = -np.sin(axis_rad) * tick_len
            tang_y =  np.cos(axis_rad) * tick_len
            tx0 = r_axis1 * np.cos(axis_rad)
            ty0 = r_axis1 * np.sin(axis_rad)
            ax.plot([tx0 - tang_x, tx0 + tang_x],
                    [ty0 - tang_y, ty0 + tang_y],
                    color="#555555", linewidth=0.7, zorder=6)
            # Max-value tick label, sitting just beyond the axis line
            s_max = max(sizes)
            unit  = _SECTOR_UNITS.get(key, "")
            max_label = f"{_fmt_axis_val(key, s_max)} {unit}".strip()
            r_label = r_axis1 + 0.06
            lab_x = r_label * np.cos(axis_rad) - tang_x * 0.6
            lab_y = r_label * np.sin(axis_rad) - tang_y * 0.6
            # Tangential rotation so label reads along the axis
            rot = (axis_ang - 90) % 360
            if rot > 90 and rot < 270:
                rot -= 180
            ax.text(lab_x, lab_y, max_label, rotation=rot,
                    fontsize=8.0, color="#444444",
                    ha="center", va="center", zorder=6)

        current_angle = end_angle - gap_deg

    # =========================================================================
    # HALF-SUNBURST LEGEND
    # =========================================================================
    
    cx, cy = 0.0, 0.0
    
    leg_colors = [MODALITY_COLORS["scrna"], MODALITY_COLORS["bulk"]]
    leg_angles = [(45, 90), (0, 45)]
    labels = ["Counts", "Datasets"]  
    
    # Draw Quarter Rings for Modality Layers
    for idx in range(2):
        r_in = r_starts[idx]
        r_out = r_ends[idx]
        
        for ang_idx, (th1, th2) in enumerate(leg_angles):
            c = blend_white(leg_colors[ang_idx], layer_shading[idx])
            w = Wedge((cx, cy), r_out, th1, th2, width=r_out - r_in,
                      facecolor=c, edgecolor=BG_COLOR, lw=2.0)
            ax_leg.add_patch(w)
            
        # Label
        r_mid = (r_in + r_out) / 2.0
        ax_leg.text(-0.06, r_mid, labels[idx], ha="right", va="center", 
                    fontsize=11, color=TEXT_DARK, fontweight="medium")
        ax_leg.plot([-0.04, -0.01], [r_mid, r_mid], color="black", lw=1.0)
        
    # Draw Species Ring Legend
    r_in = r_starts[2]
    r_out = r_ends[2]
    sp_legend_colors = [BINARY_COLORS["Human"], BINARY_COLORS["Mouse"]]
    for ang_idx, (th1, th2) in enumerate(leg_angles):
        c = sp_legend_colors[ang_idx]
        w = Wedge((cx, cy), r_out, th1, th2, width=r_out - r_in,
                  facecolor=c, edgecolor=BG_COLOR, lw=2.0)
        ax_leg.add_patch(w)
        
    r_mid = (r_in + r_out) / 2.0
    ax_leg.text(-0.06, r_mid, "Species", ha="right", va="center", 
                fontsize=11, color=TEXT_DARK, fontweight="medium")
    ax_leg.plot([-0.04, -0.01], [r_mid, r_mid], color="black", lw=1.0)
        
    # Outer bars legend
    ax_leg.text(-0.06, r_ends[-1] + 0.2, "n of datasets", ha="right", va="center",
                fontsize=11, color=TEXT_DARK, fontweight="medium")
    ax_leg.plot([-0.04, -0.01], [r_ends[-1] + 0.2, r_ends[-1] + 0.2], color="black", lw=1.0)
    
    # Dummy outer bars (varying heights to illustrate proportional sizing)
    dummy_lengths = [0.42, 0.25, 0.10, 0.38, 0.15, 0.30]
    bar_idx = 0
    for ang, c_idx in zip([30, 60], [1, 0]):
        for offset in [-6, 0, 6]:
            length = dummy_lengths[bar_idx]
            bar_idx += 1
            w = Wedge((cx, cy), r_ends[-1] + length + 0.02, ang+offset-1.8, ang+offset+1.8,
                      width=length, facecolor="#B0B0B0", edgecolor=BG_COLOR, lw=0.8, alpha=0.55)
            ax_leg.add_patch(w)

    # Radial axis line in legend (at right edge of dummy bars)
    leg_axis_ang = 70
    leg_axis_rad = np.radians(leg_axis_ang)
    leg_r0 = r_ends[-1] + 0.02
    leg_r1 = r_ends[-1] + 0.02 + BAR_LEN_MAX
    ax_leg.plot([leg_r0 * np.cos(leg_axis_rad) + cx, leg_r1 * np.cos(leg_axis_rad) + cx],
                [leg_r0 * np.sin(leg_axis_rad) + cy, leg_r1 * np.sin(leg_axis_rad) + cy],
                color="#888888", linewidth=0.5)
    # "max" label at top of legend axis
    ax_leg.text((leg_r1 + 0.04) * np.cos(leg_axis_rad) + cx,
                (leg_r1 + 0.04) * np.sin(leg_axis_rad) + cy,
                "max", fontsize=5, color="#888888", ha="center", va="center")

    # Species Context Legend
    sp_y = -0.5
    ax_leg.text(0, sp_y + 0.15, "Species", fontsize=11, fontweight="bold", color=TEXT_DARK)
    for key, color in BINARY_COLORS.items():
        patch = plt.Rectangle((0, sp_y - 0.05), 0.1, 0.1, facecolor=color, edgecolor=BG_COLOR)
        ax_leg.add_patch(patch)
        ax_leg.text(0.15, sp_y, key, fontsize=10, color=TEXT_DARK, va="center")
        sp_y -= 0.18

    # ── Save ─────────────────────────────────────────────────────────────────
    out_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__)))),
        "figures", "main", "fig1_atlas_overview", "panels",
    )
    os.makedirs(out_dir, exist_ok=True)

    pdf_path = os.path.join(out_dir, "fig1b.pdf")

    fig.savefig(pdf_path, bbox_inches="tight", dpi=600, facecolor=BG_COLOR)
    plt.close(fig)
    print(f"Saved: {pdf_path}")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Publication-quality treemap: relative sample contribution of each dataset
to the MASLD Transcriptomic Atlas.

Inner rectangles sized by total N; color encodes diet/disease model type.
Follows Sanjana Lab Publication Standards.

Output: PDF + PNG (600 DPI) to RNA-seq/figures/
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import squarify
import numpy as np

# ─── Sanjana Lab Publication Standards ───────────────────────────────────────
plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":        6,
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
    "figure.dpi":       150,
})

# ─── Color Definitions ──────────────────────────────────────────────────────
# Human cohorts — shades of pink/purple gradient
HUMAN_COLORS = {
    "MASLD Spectrum":  "#c94d84",   # Medium pink
    "NAFL vs NASH":    "#af3d6e",   # Deep pink
    "NAS Graded":      "#703d74",   # Deep purple
    "BMI Spectrum":    "#962d59",   # Dark pink
    "Duke University": "#dab5d5",   # Light purple
    "Fibrosis":        "#f176ae",   # Light pink
}

# Mouse datasets — shades of cyan/green gradient
MOUSE_COLORS = {
    "MCD":     "#2d652e",   # Deep green
    "HFD":     "#58a843",   # Medium green
    "CDAHFD":  "#095064",   # Deep cyan
    "FPC":     "#237e8e",   # Medium cyan
    "LIDPAD":  "#75c8cd",   # Light cyan
    "GAN":     "#87c168",   # Light green
}

# ─── Dataset Definitions ────────────────────────────────────────────────────
# (label, accession, year, total_n, model_type)
human_datasets = [
    ("Chen",    "GSE213621",   2023, 368, "Duke University"),
    # PRJNA512027 (Gerhard 2018) excluded: L0/S0 library-prep batch
    # confounded with diagnosis. Still used in fibrosis-vs-healthy.
    ("Govaere", "GSE135251",   2020, 216, "Fibrosis"),
    ("Kozumi",  "GSE167523",   2021,  98, "NAFL vs NASH"),
    ("Hoang",   "GSE130970",   2019,  78, "NAS Graded"),
    ("Suppli",  "GSE126848",   2019,  57, "BMI Spectrum"),
]

mouse_datasets = [
    ("LIDPAD",       "GSE159911",   2020, 151, "LIDPAD"),
    ("FPC/CDAHFD",   "GSE162876",   2020, 213, "FPC"),
    ("HFD (Long)",   "GSE274914",   2024,  22, "HFD"),
    ("HFD (Short)",  "GSE224069",   2023,  43, "HFD"),
    ("GAN",          "GSE225616",   2023,  19, "GAN"),
    ("Paquette MCD", "GSE156918",   2021,  22, "MCD"),
    ("Yue MCD",      "GSE205974",   2022,   6, "MCD"),
    ("In-House MCD", "Sanjana Lab", 2024,  12, "MCD"),
]


def draw_treemap(ax, datasets, color_map, title, title_color, species):
    """Draw a treemap for one species panel."""
    # Sort by size descending for layout
    datasets_sorted = sorted(datasets, key=lambda x: x[3], reverse=True)

    labels = []
    sizes = []
    colors = []
    short_labels = []

    for d in datasets_sorted:
        name, acc, year, n, model = d
        # Display by accession (GEO ID) instead of author/descriptive name
        disp = acc
        # Full label
        labels.append(f"{disp}\n({year})\nn={n}")
        # Short label for small boxes
        short_labels.append(f"{disp}\nn={n}")
        sizes.append(n)
        colors.append(color_map.get(model, "#cccccc"))

    # Normalize for squarify
    total = sum(sizes)
    normed = squarify.normalize_sizes(sizes, 100, 60)

    rects = squarify.squarify(normed, 0, 0, 100, 60)

    for rect, label, short_label, color, size in zip(rects, labels, short_labels, colors, sizes):
        x, y, dx, dy = rect["x"], rect["y"], rect["dx"], rect["dy"]

        # Draw rectangle
        patch = mpatches.FancyBboxPatch(
            (x + 0.3, y + 0.3), dx - 0.6, dy - 0.6,
            boxstyle="round,pad=0.3",
            facecolor=color, edgecolor="white", linewidth=2,
            alpha=0.92
        )
        ax.add_patch(patch)

        # ─── Smart Labeling Logic ──────────────────────────────────────────
        # Determine available area and dimensions
        area = dx * dy
        min_dim = min(dx, dy)
        
        # Heuristics for font size and content
        display_label = label
        if area > 400:
            fontsize = 6
        elif area > 150:
            fontsize = 6
        elif area > 60:
            fontsize = 6
        elif area > 25:
             # Very small box - use short label and small font
            fontsize = 6
            display_label = short_label
        else:
            # Extremely small - just show N or maybe nothing if really tiny
            fontsize = 6
            display_label = f"n={size}"

        # If box is too thin, rotate or hide
        if min_dim < 6 and area > 25:
             # If it's long but thin, maybe rotate? (Not doing rotation for now, just shrinking)
             fontsize = 6
        
        if min_dim < 3:
             # Too thin to print anything legible
             display_label = ""

        # Choose text color based on luminance
        text_color = "white" if color in ["#2d652e", "#095064", "#237e8e",
                                           "#af3d6e", "#703d74", "#962d59",
                                           "#c94d84"] else "#333333"

        if display_label:
            ax.text(x + dx / 2, y + dy / 2, display_label,
                    ha="center", va="center", fontsize=fontsize,
                    color=text_color,
                    linespacing=1.2)

    ax.set_xlim(0, 100)
    ax.set_ylim(0, 60)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.axis("off")

    # Total annotation
    ax.text(50, -3, f"Total: {total} samples",
            ha="center", va="center", fontsize=6,
            color=title_color)


def main():
    fig, (ax_human, ax_mouse) = plt.subplots(
        1, 2, figsize=(7.09, 2.66),
        gridspec_kw={"wspace": 0.15}
    )

    draw_treemap(ax_human, human_datasets, HUMAN_COLORS,
                 "Human Patient Cohorts", "#c94d84", "Human")
    draw_treemap(ax_mouse, mouse_datasets, MOUSE_COLORS,
                 "Mouse Diet Models", "#3a7b36", "Mouse")

    # ── Legend for Human ─────────────────────────────────────────────────
    h_patches = [mpatches.Patch(facecolor=c, edgecolor="white", label=k)
                 for k, c in HUMAN_COLORS.items()]
    ax_human.legend(handles=h_patches, loc="lower left",
                    fontsize=6, framealpha=0.9, edgecolor="#dddddd",
                    title="Disease Model", title_fontsize=6, ncol=2)

    # ── Legend for Mouse ─────────────────────────────────────────────────
    m_patches = [mpatches.Patch(facecolor=c, edgecolor="white", label=k)
                 for k, c in MOUSE_COLORS.items()]
    ax_mouse.legend(handles=m_patches, loc="lower left",
                    fontsize=6, framealpha=0.9, edgecolor="#dddddd",
                    title="Diet Type", title_fontsize=6, ncol=2)

    # ── Save ─────────────────────────────────────────────────────────────
    import os
    out_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), "figures")
    os.makedirs(out_dir, exist_ok=True)

    pdf_path = os.path.join(out_dir, "dataset_treemap.pdf")

    # RETIRED 2026-06-23 (user request): stray root-level dataset_treemap.pdf no longer generated.
    # fig.savefig(pdf_path, bbox_inches="tight", dpi=600)
    plt.close(fig)

    print("[caption] Sample Contribution Treemap - MASLD Transcriptomic Atlas. "
          "Left panel: Human Patient Cohorts. Right panel: Mouse Diet Models.")
    print("RETIRED 2026-06-23: dataset_treemap.pdf generation disabled.")


if __name__ == "__main__":
    main()

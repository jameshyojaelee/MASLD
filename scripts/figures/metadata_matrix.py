#!/usr/bin/env python3
"""
Publication-quality metadata overview for the 9 human RNA-seq cohorts
in the MASLD Transcriptomic Atlas. PRJNA512027 (Gerhard 2018) is excluded
from the cohort presentation because its L0/S0 library-prep batch is
perfectly confounded with diagnosis (all 34 controls L0, all 102 NASH S0);
it remains available for fibrosis-vs-healthy contrasts via the pipeline.

Left panel:  OncoPrint-style metadata availability matrix (7 columns)
Right panel: Stacked bar showing sample composition
             (Control / Early-Mild / Advanced-Severe) + QC-pass counts

Sorted by total N (largest on top).  Zebra-striped for readability.
Follows Sanjana Lab Publication Standards.
Output: PDF (600 DPI) to figures/
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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

# ─── Color Definitions ──────────────────────────────────────────────────────
AVAILABLE    = "#C9265E"   # Liang deep magenta — metadata fully available
PARTIAL      = "#F4A674"   # Liang warm peach — partially available / inferred
ABSENT       = "#ECECEC"   # Near-white — not available
TEXT_DARK    = "#2D3436"
LABEL_GRAY   = "#636E72"
ZEBRA        = "#f7f7f7"   # Very faint gray for alternating rows

# Stacked bar colors (disease severity progression)
COL_CONTROL  = "#9E9E9E"   # Liang neutral gray — healthy controls
COL_EARLY    = "#F4A674"   # Liang warm peach — early / mild disease
COL_ADVANCED = "#C9265E"   # Liang deep magenta — advanced / severe disease

# ─── Dataset Definitions (sorted by N, largest first) ────────────────────────
# Tuple fields:
#   0: author           7: n_control       14: fibrosis
#   1: accession        8: n_early         15: nafl_nash_dx
#   2: year             9: n_advanced      16: in_mega
#   3: n_total         10: controls (matrix val)
#   4: n_qc            11: age
#   5: platform        12: sex
#   6: layout          13: nas
#
# Sample composition mapping:
#   Control  = healthy / histologically normal
#   Early    = NAFL, F0, F1, Fibrosis_F0F1, steatosis-only
#   Advanced = NASH, NASH_Fibrosis, F2+
#
# Matrix values: 1 = available, 0.5 = partial / inferred, 0 = absent

datasets = [
    # ── author        accession       yr   N    QC   platform         layout
    #    n_ctrl n_early n_adv
    #    controls age sex nas fibrosis nafl_nash_dx in_mega

    ("Chen",       "GSE213621",   2023, 367, 359, "HiSeq 2500",    "PE",
     68,  97, 202,
     1, 0, 0.5, 0, 0.5, 0, 1),

    ("Govaere",    "GSE135251",   2020, 216, 215, "NextSeq 500",   "PE",
     10,  51, 155,
     1, 0, 0.5, 1, 1, 1, 1),

    # PRJNA512027 (Gerhard 2018) excluded from cohort presentation:
    # L0/S0 library-prep batch is perfectly confounded with diagnosis
    # (all 34 controls L0, all 102 NASH S0). Still loaded by the pipeline
    # for the fibrosis-vs-healthy contrast.

    ("Hoshida",    "GSE193066",   2022, 164, 160, "NextSeq 500",   "SE",
     0,  51, 113,
     0, 1, 1, 1, 1, 1, 0),

    ("Pantano",    "GSE162694",   2021, 143, 142, "HiSeq 3000",    "SE",
     31,  65,  47,
     1, 1, 1, 0.5, 1, 1, 1),

    ("Kozumi",     "GSE167523",   2021,  98,  96, "HiSeq 3000",    "SE",
     0,  51,  47,
     0, 1, 1, 0, 0, 0.5, 0),

    ("Kawamura",   "GSE174478",   2021,  94,  93, "NovaSeq 6000",  "PE",
     0,  36,  58,
     0, 1, 1, 1, 1, 1, 0),

    ("Hoang",      "GSE130970",   2019,  78,  76, "HiSeq 2500",    "PE",
     25,  28,  25,
     1, 1, 1, 1, 1, 1, 1),

    ("Verschuren", "GSE240729",   2023,  67,  63, "NovaSeq 6000",  "PE",
     0,  26,  41,
     0, 0, 0.5, 0, 1, 0, 0),

    ("Suppli",     "GSE126848",   2019,  57,  55, "NextSeq 500",   "SE",
     26,  15,  16,
     1, 0, 1, 0, 0, 0.5, 1),
]

# Column names for the matrix (indices 10-16)
columns = [
    "Healthy\nControls",
    "Age",
    "Sex",
    "NAS\nScore",
    "Fibrosis\nStage",
    "MASL vs\nMASH Dx",
    "Disease vs\nControl DE",
]


def main():
    n_ds  = len(datasets)
    n_col = len(columns)

    # ── Extract fields ────────────────────────────────────────────────────
    matrix     = np.array([[d[i] for i in range(10, 10 + n_col)]
                           for d in datasets])
    n_controls = [d[7]  for d in datasets]
    n_early    = [d[8]  for d in datasets]
    n_advanced = [d[9]  for d in datasets]
    n_totals   = [d[3]  for d in datasets]
    n_qc       = [d[4]  for d in datasets]

    # ── Figure layout ─────────────────────────────────────────────────────
    fig = plt.figure(figsize=(14, 4.5))
    gs  = fig.add_gridspec(1, 2, width_ratios=[3.2, 1.5], wspace=0.05)
    ax_mat = fig.add_subplot(gs[0])
    ax_bar = fig.add_subplot(gs[1])

    # ── Zebra striping (both panels) ──────────────────────────────────────
    for row in range(n_ds):
        if row % 2 == 1:
            ax_mat.axhspan(row - 0.5, row + 0.5, color=ZEBRA, zorder=0)
            ax_bar.axhspan(row - 0.5, row + 0.5, color=ZEBRA, zorder=0)

    # ── Draw matrix cells ─────────────────────────────────────────────────
    cell_w, cell_h = 0.72, 0.88
    for row in range(n_ds):
        for col in range(n_col):
            val = matrix[row, col]
            if val == 1:
                color = AVAILABLE
            elif val == 0.5:
                color = PARTIAL
            else:
                color = ABSENT

            rect = mpatches.FancyBboxPatch(
                (col - cell_w / 2, row - cell_h / 2), cell_w, cell_h,
                boxstyle="round,pad=0.04",
                facecolor=color, edgecolor="white", linewidth=1.5,
            )
            ax_mat.add_patch(rect)

    # ── Matrix axes ───────────────────────────────────────────────────────
    ax_mat.set_xlim(-0.5, n_col - 0.5)
    ax_mat.set_ylim(n_ds - 0.5, -0.5)
    ax_mat.set_xticks(range(n_col))
    ax_mat.set_xticklabels(columns, fontsize=9, fontweight="bold",
                           ha="center", color=TEXT_DARK)
    ax_mat.xaxis.set_ticks_position("top")
    ax_mat.xaxis.set_label_position("top")

    # ── Two-tier row labels (separate text for different font sizes) ──────
    ax_mat.set_yticks(range(n_ds))
    ax_mat.set_yticklabels([""] * n_ds)     # clear default labels
    ax_mat.tick_params(axis="y", length=0)

    label_x = -0.65                         # right edge of label area
    for i, d in enumerate(datasets):
        # Line 1: Accession (Year) — bold, larger
        ax_mat.text(label_x, i - 0.15,
                    f"{d[1]} ({d[2]})",
                    ha="right", va="center", fontsize=9,
                    fontweight="bold", color=TEXT_DARK)
        # Line 2: Platform · Layout — smaller, gray, italic
        ax_mat.text(label_x, i + 0.17,
                    f"{d[5]}  \u00b7  {d[6]}",
                    ha="right", va="center", fontsize=6.5,
                    color=LABEL_GRAY, fontstyle="italic")

    for spine in ax_mat.spines.values():
        spine.set_visible(False)
    ax_mat.tick_params(axis="x", length=0)

    # ── Coverage summary — tight row just below last matrix cell ─────────
    cov_y = n_ds - 0.15
    for ci in range(n_col):
        full    = int(sum(1 for r in range(n_ds) if matrix[r, ci] == 1))
        partial = int(sum(1 for r in range(n_ds) if matrix[r, ci] == 0.5))
        txt = f"{full}+{partial}/9" if partial else f"{full}/9"
        ax_mat.text(ci, cov_y, txt,
                    ha="center", va="top", fontsize=7,
                    color=LABEL_GRAY, fontstyle="italic")

    ax_mat.text(label_x, cov_y, "Coverage:",
                ha="right", va="top", fontsize=7,
                color=LABEL_GRAY, fontstyle="italic", fontweight="bold")

    # ── Legend (metadata status) — well below coverage row ────────────────
    legend_patches = [
        mpatches.Patch(facecolor=AVAILABLE, edgecolor="white",
                       label="Available"),
        mpatches.Patch(facecolor=PARTIAL,   edgecolor="white",
                       label="Partial / Inferred"),
        mpatches.Patch(facecolor=ABSENT,    edgecolor="white",
                       label="Not Available"),
    ]
    ax_mat.legend(handles=legend_patches, loc="lower center",
                  framealpha=0.95, fontsize=8, edgecolor="#dddddd",
                  title="Metadata Status", title_fontsize=9,
                  bbox_to_anchor=(0.5, -0.22), ncol=3)

    # ── Right panel: stacked bar (Control / Early / Advanced) ─────────────
    y_pos      = np.arange(n_ds)
    bar_height = 0.75

    ax_bar.barh(y_pos, n_controls, height=bar_height,
                color=COL_CONTROL, edgecolor="white", linewidth=0.5,
                zorder=3, label="Control")

    ax_bar.barh(y_pos, n_early, left=n_controls, height=bar_height,
                color=COL_EARLY, edgecolor="white", linewidth=0.5,
                zorder=3, label="Early / Mild")

    left_adv = [c + e for c, e in zip(n_controls, n_early)]
    ax_bar.barh(y_pos, n_advanced, left=left_adv, height=bar_height,
                color=COL_ADVANCED, edgecolor="white", linewidth=0.5,
                zorder=3, label="Advanced / Severe")

    # Annotate with n= (QC) on single line
    max_n = max(n_totals)
    for i in range(n_ds):
        ax_bar.text(n_totals[i] + max_n * 0.02, y_pos[i],
                    f"n={n_totals[i]}  ({n_qc[i]} QC)",
                    va="center", ha="left", fontsize=7.5,
                    fontweight="bold", color=TEXT_DARK)

    ax_bar.set_ylim(n_ds - 0.5, -0.5)
    ax_bar.set_yticks([])
    ax_bar.set_xlabel("Samples", fontsize=9, color=LABEL_GRAY)
    ax_bar.set_xlim(0, max_n * 1.50)
    ax_bar.spines["left"].set_visible(False)
    ax_bar.spines["bottom"].set_linewidth(0.8)
    ax_bar.spines["bottom"].set_color("#999999")
    # Liang style: no dashed gridlines on stacked bars
    ax_bar.set_title("Sample Composition", fontsize=11, fontweight="bold",
                     pad=10, color=TEXT_DARK)

    # Bar legend (disease stages)
    bar_patches = [
        mpatches.Patch(facecolor=COL_CONTROL,  edgecolor="white",
                       label="Control"),
        mpatches.Patch(facecolor=COL_EARLY,    edgecolor="white",
                       label="Early / Mild"),
        mpatches.Patch(facecolor=COL_ADVANCED, edgecolor="white",
                       label="Advanced / Severe"),
    ]
    ax_bar.legend(handles=bar_patches, loc="lower right",
                  framealpha=0.95, fontsize=7.5, edgecolor="#dddddd",
                  title="Disease Stage", title_fontsize=8,
                  bbox_to_anchor=(1.0, -0.22))

    # ── Summary box ───────────────────────────────────────────────────────
    total_n  = sum(n_totals)
    total_qc = sum(n_qc)
    ax_bar.text(
        0.97, 0.03,
        f"9 cohorts\n{total_n:,} samples\n{total_qc:,} QC-passing",
        transform=ax_bar.transAxes, ha="right", va="bottom",
        fontsize=9, color=TEXT_DARK, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                  edgecolor="#999999", alpha=0.9),
    )

    # ── Title ─────────────────────────────────────────────────────────────
    fig.suptitle(
        "Human Cohort Metadata Overview \u2014 MASLD Transcriptomic Atlas",
        fontsize=14, fontweight="bold", y=1.02, color=TEXT_DARK,
    )

    # ── Save ──────────────────────────────────────────────────────────────
    out_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__)))),
        "figures",
    )
    os.makedirs(out_dir, exist_ok=True)
    pdf_path = os.path.join(out_dir, "metadata_matrix.pdf")

    # RETIRED 2026-06-23 (user request): this stray root-level metadata_matrix.pdf is no
    # longer generated. The kept metadata table is the fig3 panel from fig_bulkrna_matrix.py
    # -> figures/main/fig3_RNAseq/panels/cohort_metadata_matrix.pdf.
    # fig.savefig(pdf_path, bbox_inches="tight", dpi=600)
    plt.close(fig)
    print("RETIRED 2026-06-23: metadata_matrix.pdf generation disabled (see cohort_metadata_matrix.pdf).")


if __name__ == "__main__":
    main()

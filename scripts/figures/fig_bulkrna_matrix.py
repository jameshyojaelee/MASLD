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

Usage:
  python fig_bulkrna_matrix.py            # standalone (full labels)
  python fig_bulkrna_matrix.py --panel    # fig1 panel version (accession-only y-labels)

Output: PDF (600 DPI) to figures/
"""

import argparse
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
    "font.size":        6,
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
    "axes.spines.top":  False,
    "axes.spines.right": False,
    "figure.dpi":       150,
})

# Equalized publication font sizes (single source of truth for this script).
FS_LABEL    = 6        # tick / row labels, axis annotations
FS_TITLE    = 6        # column headers, panel titles
FS_LEGEND   = 6        # legend body
FS_LEGEND_T = 6        # legend title
FS_SUMMARY  = 6        # summary box

# ─── Color Definitions — Liang et al. (Cas13 lncRNA) palette ────────────────
# Availability matrix
AVAILABLE    = "#2E86AB"   # teal-blue (reads as "yes / present")
PARTIAL      = "#F4A674"   # Liang peach (partial / weak signal)
ABSENT       = "#ECECEC"   # neutral gray (no data)
TEXT_DARK    = "#2D3436"
LABEL_GRAY   = "#636E72"
ZEBRA        = "#f7f7f7"   # Very faint gray for alternating rows

# Stacked bar — disease severity (control = neutral gray per project convention)
# Distinct hue family from left panel's magenta/peach to avoid cross-legend confusion.
COL_CONTROL  = "#9E9E9E"   # canonical control gray (matches fig2 panels)
COL_EARLY    = "#FDB462"   # amber-orange (severity ramp; distinct from left-panel peach)
COL_ADVANCED = "#D62B2B"   # alarm red (distinct from left-panel magenta)

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

    ("Bril",       "GSE162694",   2021, 143, 142, "HiSeq 3000",    "SE",
     31,  65,  47,
     1, 1, 1, 0.5, 1, 1, 1),   # fibrosis=Available: ALL 112 disease patients precisely staged F0-F4 (only the 31 controls are NA). NAS=Partial: 26/112 PATIENTS genuinely lack NAS.

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
]



def main(panel_mode=False):
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
    # Panel mode (fig3a, 2026-07-09 resize to 2.88x1.93in) is much more compact
    # than the standalone supp version, so margins are set EXPLICITLY (absolute
    # inches -> figure-fraction) rather than relying on matplotlib's default
    # rcParams margins, which were sized for the old 6.8x2.6in canvas and would
    # clip the row-accession labels / column headers at this smaller size.
    figsize = (2.88, 1.93) if panel_mode else (7.09, 2.1)
    fig = plt.figure(figsize=figsize)
    if panel_mode:
        gs = fig.add_gridspec(1, 2, width_ratios=[2.5, 1.3], wspace=0.06,
                              left=0.175, right=0.98, top=0.87, bottom=0.40)
    else:
        gs = fig.add_gridspec(1, 2, width_ratios=[2.3, 1.5], wspace=0.05)
    ax_mat = fig.add_subplot(gs[0])
    ax_bar = fig.add_subplot(gs[1])

    # ── Zebra striping (both panels) ──────────────────────────────────────
    for row in range(n_ds):
        if row % 2 == 1:
            ax_mat.axhspan(row - 0.5, row + 0.5, color=ZEBRA, zorder=0)
            ax_bar.axhspan(row - 0.5, row + 0.5, color=ZEBRA, zorder=0)

    # ── Draw matrix cells ─────────────────────────────────────────────────
    cell_w, cell_h = 0.72, 0.78
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
    ax_mat.set_xticklabels(columns, fontsize=FS_TITLE, fontweight="normal",
                           ha="center", color=TEXT_DARK)
    ax_mat.xaxis.set_ticks_position("top")
    ax_mat.xaxis.set_label_position("top")

    # ── Row labels ────────────────────────────────────────────────────────
    ax_mat.set_yticks(range(n_ds))
    ax_mat.set_yticklabels([""] * n_ds)     # clear default labels
    ax_mat.tick_params(axis="y", length=0)

    label_x = -0.65                         # right edge of label area
    for i, d in enumerate(datasets):
        if panel_mode:
            # Panel mode: accession only, single centered line
            ax_mat.text(label_x, i,
                        d[1],
                        ha="right", va="center", fontsize=FS_LABEL,
                        fontweight="normal", color=TEXT_DARK)
        else:
            # Standalone: two-tier labels
            # Line 1: Accession (Year)
            ax_mat.text(label_x, i - 0.15,
                        f"{d[1]} ({d[2]})",
                        ha="right", va="center", fontsize=FS_LABEL,
                        fontweight="normal", color=TEXT_DARK)
            # Line 2: Platform · Layout — smaller, gray
            ax_mat.text(label_x, i + 0.17,
                        f"{d[5]}  \u00b7  {d[6]}",
                        ha="right", va="center", fontsize=FS_LABEL,
                        color=LABEL_GRAY, fontstyle="normal")

    for spine in ax_mat.spines.values():
        spine.set_visible(False)
    ax_mat.tick_params(axis="x", length=0)


    # ── Legend (metadata status) ────────────────────────────────────────────
    # 2026-07-09 fix (round 2): a FIGURE-level anchor (fig.legend, 1st attempt)
    # doesn't reserve layout space and drew on top of other content. A 3-col
    # single-row layout on ax_mat (2nd attempt) was WIDER than ax_mat's own
    # ~1.5in column and got clipped at the panel edge ("Not Available" cut
    # off). Single column (ncol=1): each row only needs room for ONE entry
    # ("Partial / Inferred", the widest, ~0.8in) -- comfortably fits ax_mat's
    # width. 3 stacked rows (~0.3in tall) placed via a margin-derived offset,
    # centered in the bottom=0.20 gridspec margin reserved above.
    legend_patches = [
        mpatches.Patch(facecolor=AVAILABLE, edgecolor="white",
                       label="Available"),
        mpatches.Patch(facecolor=PARTIAL,   edgecolor="white",
                       label="Partial / Inferred"),
        mpatches.Patch(facecolor=ABSENT,    edgecolor="white",
                       label="Not Available"),
    ]
    # 2026-07-09 fix (round 3): rounds 1-2 both underestimated the legend's
    # true rendered height and it clipped off the bottom of the page. Bottom
    # margin raised 0.20->0.28 (generous headroom) and anchored flush just
    # below the axes edge (small fixed gap) rather than "centered" in the
    # margin -- removes the height-estimation step entirely; whatever the
    # true height, it grows down into the now much larger reserved margin.
    ax_mat.legend(handles=legend_patches, loc="upper center",
                  bbox_to_anchor=(0.5, -0.03), bbox_transform=ax_mat.transAxes,
                  framealpha=0.95, fontsize=FS_LEGEND, edgecolor="#dddddd",
                  ncol=1, handlelength=1.1, handleheight=1.0,
                  handletextpad=0.35, borderpad=0.3, labelspacing=0.35)

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

    max_n = max(n_totals)

    ax_bar.set_ylim(n_ds - 0.5, -0.5)
    ax_bar.set_yticks([])
    ax_bar.set_xlabel("Samples", fontsize=FS_TITLE, color=LABEL_GRAY)
    ax_bar.set_xlim(0, max_n * 1.32)
    ax_bar.spines["left"].set_visible(False)
    ax_bar.spines["bottom"].set_linewidth(0.6)
    ax_bar.spines["bottom"].set_color("#999999")

    # Bar legend (disease stages)
    bar_patches = [
        mpatches.Patch(facecolor=COL_CONTROL,  edgecolor="white",
                       label="Control"),
        mpatches.Patch(facecolor=COL_EARLY,    edgecolor="white",
                       label="Early / Mild"),
        mpatches.Patch(facecolor=COL_ADVANCED, edgecolor="white",
                       label="Advanced / Severe"),
    ]
    # 2026-07-09 (round 6): rounds 4-5 kept nudging the offset/margin further
    # but the "Samples" xlabel kept faintly touching the legend's TITLE row --
    # growing the bottom margin further to chase it risks shrinking the 9
    # matrix rows above the Shah-Samworth^Wsingle-line-text floor (7.2pt/row
    # at 6pt) instead. Smarter fix: drop the title entirely (the 3 color
    # swatches read fine in context, next to the bars they describe -- the
    # Metadata Status legend already has no title and is unambiguous the same
    # way), removing a full row + title-gap rather than continuing to grow
    # margins that trade one collision risk for another.
    # 2026-07-09 (round 7): a hairline of the "Samples" xlabel still peeked
    # through the box top (framealpha=0.95 lets a sliver bleed through even
    # at near-zero geometric overlap). Nudged the gap slightly further
    # (-0.20->-0.24) and set framealpha=1.0 (fully opaque) as a safety net --
    # any remaining sub-point overlap is now visually a clean edge, not a
    # ghost.
    ax_bar.legend(handles=bar_patches, loc="upper center",
                  bbox_to_anchor=(0.5, -0.32), bbox_transform=ax_bar.transAxes,
                  framealpha=1.0, fontsize=FS_LEGEND, edgecolor="#dddddd",
                  handlelength=1.1, handleheight=1.0,
                  handletextpad=0.35, borderpad=0.3, labelspacing=0.35)

    # ── Summary box ───────────────────────────────────────────────────────
    # 2026-07-09 (round 5): round 4's single-line "9 cohorts, 1,284 samples"
    # (~1.0in at 6pt) was WIDER than ax_bar's own column (~0.76in), so even
    # right-aligned to ax_bar's edge it spilled left into ax_mat's "Fibrosis
    # Stage" header. Back to 2 lines (longest line "1,284 samples" ~0.58in,
    # comfortably under 0.76in) so it stays within ax_bar's own width.
    total_n  = sum(n_totals)
    total_qc = sum(n_qc)
    ax_bar.text(
        0.97, 1.03,
        f"9 cohorts\n{total_n:,} samples",
        transform=ax_bar.transAxes, ha="right", va="bottom",
        fontsize=FS_SUMMARY, color=TEXT_DARK, fontweight="normal",
    )

    # ── Title ─────────────────────────────────────────────────────────────

    # ── Save ──────────────────────────────────────────────────────────────
    fig_root = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__)))),
        "figures", "main",
    )
    if panel_mode:
        # RESTORED 2026-07-02 (A–J layout): Fig 3A cohort-metadata matrix. Shares
        # its underlying data with Fig 1's bulkrna_metadata_matrix (both from this
        # script) — intentional: the cohort overview opens both the atlas figure
        # (Fig 1) and the RNA-seq figure (Fig 3A, top-left of row 1).
        out_dir = os.path.join(fig_root, "fig3_RNAseq", "panels")
        os.makedirs(out_dir, exist_ok=True)
        pdf_path = os.path.join(out_dir, "fig3a_cohort_metadata_matrix.pdf")
    else:
        # bulkrna_metadata_matrix.pdf is deprecated and must not be recreated.
        print("Note: standalone mode without --panel is deprecated. bulkrna_metadata_matrix.pdf is no longer generated.")
        print("Please run with --panel to generate the canonical fig3a_cohort_metadata_matrix.pdf.")
        plt.close(fig)
        return

    fig.savefig(pdf_path, dpi=600)   # NO bbox_inches="tight": preserve exact figsize for place-at-100%
    plt.close(fig)
    print(f"Saved: {pdf_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", action="store_true",
                        help="Panel mode: accession-only y-labels, saves to bulkrna_metadata_panel.pdf")
    args = parser.parse_args()
    main(panel_mode=args.panel)

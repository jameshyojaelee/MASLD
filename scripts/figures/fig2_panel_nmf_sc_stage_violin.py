#!/usr/bin/env python3
"""fig2_panel_nmf_sc_stage_violin.py

Violin plots: bulk NMF k=6 program scores (from sc.tl.score_genes) in
individual scRNA-seq cells, stratified by MASLD disease stage coarse
(Healthy / Steatosis / Steatohepatitis / Cirrhosis).

Programs shown: P1 Stromal, P2 Inflammatory, P3 Fibrogenic, P6 Kupffer-cell
Excluded: P4 Quiescent-1, P5 Quiescent-2

Input:
  figures/misc/nmf_umap_exploration/subsampled_cells.csv.gz
    (bulk_P1..P6 scores per scRNA cell; 200K subsample)
  Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv
    (disease_stage_coarse per sample)

Output:
  figures/main/fig2_progression_sex/panels/fig2_panel_nmf_sc_stage_violin.pdf
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))

SUBSAMPLE_CSV = BASE / "figures/misc/nmf_umap_exploration/subsampled_cells.csv.gz"
DONOR_META    = BASE / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv"
OUT_PDF       = BASE / "figures/main/fig2_progression_sex/panels/fig2_panel_nmf_sc_stage_violin.pdf"

# ---------------------------------------------------------------------------
# Program metadata (canonical labels from fig2_panel_nmf_programs_lines.R)
# Only 4 non-quiescent programs.
# ---------------------------------------------------------------------------
PROGRAMS = [
    {"col": "bulk_P1", "label": "P1\nStromal",       "color": "#F4A674"},
    {"col": "bulk_P2", "label": "P2\nInflammatory",   "color": "#C9265E"},
    {"col": "bulk_P3", "label": "P3\nFibrogenic",     "color": "#1565C0"},
    {"col": "bulk_P6", "label": "P6\nKupffer-cell",   "color": "#00695C"},
]

# Stage order + display labels
# Cirrhosis is dropped because the only Cirrhosis donors in the atlas come from
# GSE136103 (NPC-enriched, excluded below); after exclusion no Cirrhosis cells remain.
STAGE_ORDER  = ["Healthy", "Steatosis", "Steatohepatitis"]
STAGE_LABELS = ["Healthy", "ST", "SH"]
STAGE_COLORS = {
    "Healthy":        "#9E9E9E",   # neutral gray (control)
    "Steatosis":      "#F4A674",   # warm peach (MASL)
    "Steatohepatitis":"#C9265E",   # Liang magenta (MASH)
}

# Protocol-contaminated datasets (exclude from all cross-stage comparisons):
#  - GSE136103: NPC-enrichment (Ramachandran 2019 Nature), 0% hepatocytes
#  - Liver_Atlas: FACS-sorted macrophages (Guilliams 2022 Cell)
EXCLUDE_DATASETS = {"GSE136103", "Liver_Atlas"}

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print("Loading subsampled cells …", flush=True)
df = pd.read_csv(SUBSAMPLE_CSV)
if "index" in df.columns and "barcode" not in df.columns:
    df = df.rename(columns={"index": "barcode"})

print(f"  cells: {len(df):,}  columns: {list(df.columns)}")

print("Loading donor metadata …", flush=True)
meta = pd.read_csv(DONOR_META, sep="\t")[["sample", "disease_stage_coarse"]]
print(f"  donors: {len(meta):,}")
print(f"  stage counts:\n{meta['disease_stage_coarse'].value_counts().to_string()}")

# Merge disease stage onto cells via sample column
df = df.merge(meta, on="sample", how="left")
missing = df["disease_stage_coarse"].isna().sum()
print(f"  cells after merge: {len(df):,}  missing stage: {missing:,}")
df = df.dropna(subset=["disease_stage_coarse"])
df = df[df["disease_stage_coarse"].isin(STAGE_ORDER)]
# Drop protocol-contaminated datasets (NPC-enriched / FACS-sorted)
before = len(df)
df = df[~df["dataset"].isin(EXCLUDE_DATASETS)]
print(f"  excluded {EXCLUDE_DATASETS}: removed {before - len(df):,} cells")
print(f"  cells with valid stage: {len(df):,}")
print(df["disease_stage_coarse"].value_counts().to_string())

# ---------------------------------------------------------------------------
# Winsorise scores at 1st/99th per program to trim outliers
# ---------------------------------------------------------------------------
for p in PROGRAMS:
    col = p["col"]
    lo, hi = df[col].quantile([0.01, 0.99])
    df[col] = df[col].clip(lo, hi)

# ---------------------------------------------------------------------------
# Figure layout: 1 row × 4 programs
# ---------------------------------------------------------------------------
matplotlib.rcParams.update({
    "font.family":    "Arial",
    "font.size":      7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size":  2.5,
    "ytick.major.size":  2.5,
    "pdf.fonttype": 42,
    "ps.fonttype":  42,
})

n_prog = len(PROGRAMS)
fig_w_mm = 170   # total figure width
fig_h_mm =  65   # total figure height
fig, axes = plt.subplots(
    1, n_prog,
    figsize=(fig_w_mm / 25.4, fig_h_mm / 25.4),
    sharey=False,
)

for ax, prog in zip(axes, PROGRAMS):
    col   = prog["col"]
    label = prog["label"]
    color = prog["color"]

    # Per-stage score vectors in stage order
    data_by_stage = [
        df.loc[df["disease_stage_coarse"] == s, col].values
        for s in STAGE_ORDER
    ]
    n_by_stage = [len(d) for d in data_by_stage]

    # Violin
    parts = ax.violinplot(
        data_by_stage,
        positions=range(len(STAGE_ORDER)),
        widths=0.65,
        showmedians=False,
        showextrema=False,
    )
    for body in parts["bodies"]:
        body.set_facecolor(color)
        body.set_alpha(0.28)
        body.set_edgecolor(color)
        body.set_linewidth(0.4)

    # Median + IQR box overlay
    for i, (vals, stage) in enumerate(zip(data_by_stage, STAGE_ORDER)):
        q25, med, q75 = np.percentile(vals, [25, 50, 75])
        # IQR box
        ax.bar(i, q75 - q25, bottom=q25, width=0.18,
               color=STAGE_COLORS[stage], alpha=0.85, linewidth=0,
               zorder=3)
        # Median line
        ax.plot([i - 0.13, i + 0.13], [med, med],
                color="white", linewidth=1.2, solid_capstyle="round", zorder=4)

    # Axes decoration
    ax.set_xticks(range(len(STAGE_ORDER)))
    ax.set_xticklabels(STAGE_LABELS, fontsize=6.5, fontweight="bold")
    ax.set_xlim(-0.65, len(STAGE_ORDER) - 0.35)

    ax.set_title(label, fontsize=7.5, fontweight="bold", pad=4,
                 color=color)
    ax.set_ylabel("Program score" if ax is axes[0] else "", fontsize=6.5)

    # Light horizontal gridlines
    ax.yaxis.grid(True, linewidth=0.3, color="#E0E0E0", zorder=0)
    ax.set_axisbelow(True)

    # Remove top + right spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # N labels beneath x-tick labels (axes fraction coords; clip_on=False so they
    # appear even if below the axes box)
    for i, n in enumerate(n_by_stage):
        # Map data x-pos to axes fraction
        x_frac = (i - ax.get_xlim()[0]) / (ax.get_xlim()[1] - ax.get_xlim()[0])
        ax.text(x_frac, -0.22, f"n={n:,}",
                ha="center", va="top", fontsize=4.5, color="#757575",
                transform=ax.transAxes, clip_on=False)

fig.suptitle(
    "Bulk NMF program scores in individual scRNA-seq cells by disease stage",
    fontsize=7, y=1.01, fontweight="bold")

plt.tight_layout(w_pad=1.2)

OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_PDF, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"[saved] {OUT_PDF}", flush=True)

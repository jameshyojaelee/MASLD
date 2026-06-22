#!/usr/bin/env python3
"""fig2_panel_sc_stage_v2.py

Two complementary panels for bulk NMF × scRNA-seq cross-modal validation.

Panel A — Cell-type composition across MASLD stages (Option A)
  Stacked bar: fraction of each cell type per disease stage coarse.
  This is the single-cell analog of the bulk NMF stacked bar — the bulk
  programs largely capture composition shifts, not per-cell expression.

Panel C — Hepatocyte-intrinsic bulk NMF program scores, 3 stages (Option C)
  Violin + per-donor pseudobulk dots for programs P1/P2/P3/P6 in
  Hepatocytes only (Cirrhosis excluded: n=5 hepatocytes).
  Shows what each bulk program contributes at the hepatocyte level.

Outputs (figures/main/fig3_RNAseq/panels/):
  nmf_sc_celltype_composition.pdf   (panel A)
  nmf_sc_hep_violin.pdf             (panel C)
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

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))

SUBSAMPLE_CSV = BASE / "figures/misc/nmf_umap_exploration/subsampled_cells.csv.gz"
DONOR_META    = BASE / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv"
OUT_DIR       = BASE / "figures/main/fig3_RNAseq/panels"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Shared constants
# ---------------------------------------------------------------------------
# GSE136103 is an NPC-enrichment experiment (Ramachandran 2019 Nature): 0%
# hepatocytes in every sample by design (FACS-sorted for cholangiocytes /
# endothelial / macrophages). It is the ONLY cirrhosis dataset in the atlas.
# Including it in cross-stage comparisons creates a spurious hepatocyte-depletion
# artifact in Cirrhosis AND inflates Healthy macrophage fraction (28.8% vs <3%
# in all other datasets). Both panels exclude GSE136103 and drop Cirrhosis.
EXCLUDE_DATASETS = {"GSE136103", "Liver_Atlas"}

STAGE_ORDER  = ["Healthy", "Steatosis", "Steatohepatitis"]
STAGE_SHORT  = ["Healthy", "ST", "SH"]
STAGE_COLORS = {
    "Healthy":        "#9E9E9E",
    "Steatosis":      "#F4A674",
    "Steatohepatitis":"#C9265E",
}

# Cell-type palette (consistent with atlas annotations)
CT_ORDER = ["Hepatocytes", "Macrophages", "Fibroblasts",
            "Endothelial cells", "Cholangiocytes"]
CT_COLORS = {
    "Hepatocytes":       "#F4A674",   # warm peach (parenchyma)
    "Macrophages":       "#C9265E",   # magenta (inflammatory)
    "Fibroblasts":       "#1565C0",   # deep blue (fibrogenic)
    "Endothelial cells": "#7B9E87",   # muted green
    "Cholangiocytes":    "#9C6B9E",   # muted violet
}

# NMF programs (excluding Quiescent-1/2)
PROGRAMS = [
    {"col": "bulk_P1", "label": "P1\nStromal",     "color": "#F4A674"},
    {"col": "bulk_P2", "label": "P2\nInflammatory", "color": "#C9265E"},
    {"col": "bulk_P3", "label": "P3\nFibrogenic",   "color": "#1565C0"},
    {"col": "bulk_P6", "label": "P6\nKupffer-cell", "color": "#00695C"},
]

matplotlib.rcParams.update({
    "font.family":    "DejaVu Sans",
    "font.size":      7,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size":  2.5,
    "ytick.major.size":  2.5,
    "pdf.fonttype": 42,
    "ps.fonttype":  42,
})


# ---------------------------------------------------------------------------
# Data loading (shared)
# ---------------------------------------------------------------------------
def load_data():
    df = pd.read_csv(SUBSAMPLE_CSV)
    if "index" in df.columns and "barcode" not in df.columns:
        df = df.rename(columns={"index": "barcode"})
    meta = pd.read_csv(DONOR_META, sep="\t")[["sample", "disease_stage_coarse"]]
    df = df.merge(meta, on="sample", how="left").dropna(subset=["disease_stage_coarse"])
    # Drop NPC-enriched dataset (0% hepatocytes by experimental design)
    before = len(df)
    df = df[~df["dataset"].isin(EXCLUDE_DATASETS)]
    print(f"Excluded {EXCLUDE_DATASETS}: removed {before - len(df):,} cells")
    df = df[df["disease_stage_coarse"].isin(STAGE_ORDER)]
    print(f"Loaded {len(df):,} cells across {STAGE_ORDER}")
    print(df["disease_stage_coarse"].value_counts().reindex(STAGE_ORDER).to_string())
    # Dataset breakdown
    print("\nDatasets per stage:")
    print(df.groupby(["disease_stage_coarse","dataset"]).size()
            .reset_index(name="n_cells").to_string())
    return df


# ---------------------------------------------------------------------------
# Panel A: cell-type composition stacked bar
# ---------------------------------------------------------------------------
def plot_panel_A(df: pd.DataFrame):
    # Composition per stage
    comp = (df.groupby(["disease_stage_coarse", "cell_type"])
              .size()
              .unstack(fill_value=0))
    # Keep only the 5 main cell types; lump any "other" if present
    for ct in CT_ORDER:
        if ct not in comp.columns:
            comp[ct] = 0
    other_cols = [c for c in comp.columns if c not in CT_ORDER]
    if other_cols:
        comp["Other"] = comp[other_cols].sum(axis=1)
        comp = comp.drop(columns=other_cols)
        if "Other" not in CT_ORDER:
            comp["Other"] = comp.get("Other", 0)
    comp = comp[CT_ORDER]
    comp_frac = comp.div(comp.sum(axis=1), axis=0).reindex(STAGE_ORDER)

    # N cells per stage (for annotation)
    n_cells = df.groupby("disease_stage_coarse").size().reindex(STAGE_ORDER)
    n_donors = (df.groupby(["disease_stage_coarse", "sample"])
                  .size()
                  .groupby(level=0)
                  .size()
                  .reindex(STAGE_ORDER))

    fig, ax = plt.subplots(figsize=(75 / 25.4, 72 / 25.4))

    bottom = np.zeros(len(STAGE_ORDER))
    x = np.arange(len(STAGE_ORDER))
    for ct in CT_ORDER:
        vals = comp_frac[ct].values
        color = CT_COLORS.get(ct, "#BDBDBD")
        ax.bar(x, vals, bottom=bottom, width=0.68,
               color=color, linewidth=0, label=ct)
        # Label segments >7% with fraction text
        for xi, (v, b) in enumerate(zip(vals, bottom)):
            if v > 0.07:
                ax.text(xi, b + v / 2, f"{v:.0%}",
                        ha="center", va="center",
                        fontsize=4.5, color="white", fontweight="bold")
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{s}\n(n={n_donors[STAGE_ORDER[i]]})" for i, s in enumerate(STAGE_SHORT)],
        fontsize=6, fontweight="bold")
    ax.set_ylabel("Fraction of cells", fontsize=6.5)
    ax.set_ylim(0, 1)
    ax.set_xlim(-0.55, len(STAGE_ORDER) - 0.45)
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1, decimals=0))
    ax.yaxis.grid(True, linewidth=0.3, color="#E0E0E0", zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.set_title("Cell-type composition\nacross MASLD stages (5 datasets)",
                 fontsize=7.5, fontweight="bold", pad=4)

    handles = [mpatches.Patch(facecolor=CT_COLORS[ct], label=ct) for ct in CT_ORDER]
    ax.legend(handles=handles, fontsize=5, frameon=False,
              loc="lower left", bbox_to_anchor=(1.02, 0.0),
              handlelength=1.0, handletextpad=0.4)

    fig.tight_layout()
    outpath = OUT_DIR / "nmf_sc_celltype_composition.pdf"
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {outpath}")

    # Print composition table for records
    print("\nComposition table (%):")
    print((comp_frac * 100).round(1).to_string())


# ---------------------------------------------------------------------------
# Panel C: hepatocyte-only violin, 3 stages, per-donor pseudobulk overlay
# ---------------------------------------------------------------------------
def plot_panel_C(df: pd.DataFrame):
    # df is already filtered to STAGE_ORDER (Healthy/ST/SH) and excludes GSE136103
    hep = df[df["cell_type"] == "Hepatocytes"].copy()

    # Winsorise at 1/99th per program
    for p in PROGRAMS:
        col = p["col"]
        lo, hi = hep[col].quantile([0.01, 0.99])
        hep[col] = hep[col].clip(lo, hi)

    # Per-donor pseudobulk median (for dot overlay)
    pb = (hep.groupby(["sample", "disease_stage_coarse"])
             [[p["col"] for p in PROGRAMS]]
             .median()
             .reset_index())

    fig, axes = plt.subplots(1, len(PROGRAMS),
                             figsize=(155 / 25.4, 62 / 25.4),
                             sharey=False)

    for ax, prog in zip(axes, PROGRAMS):
        col   = prog["col"]
        label = prog["label"]
        color = prog["color"]

        cell_data = [hep.loc[hep["disease_stage_coarse"] == s, col].values
                     for s in STAGE_ORDER]
        donor_data = [pb.loc[pb["disease_stage_coarse"] == s, col].values
                      for s in STAGE_ORDER]
        n_cells   = [len(d) for d in cell_data]
        n_donors  = [len(d) for d in donor_data]

        # Violin (all cells)
        parts = ax.violinplot(cell_data,
                              positions=range(len(STAGE_ORDER)),
                              widths=0.6,
                              showmedians=False,
                              showextrema=False)
        for body in parts["bodies"]:
            body.set_facecolor(color)
            body.set_alpha(0.20)
            body.set_edgecolor(color)
            body.set_linewidth(0.4)

        # IQR box + median from per-cell data
        for i, (vals, stage) in enumerate(zip(cell_data, STAGE_ORDER)):
            q25, med, q75 = np.percentile(vals, [25, 50, 75])
            ax.bar(i, q75 - q25, bottom=q25, width=0.15,
                   color=STAGE_COLORS[stage], alpha=0.85, linewidth=0, zorder=3)
            ax.plot([i - 0.11, i + 0.11], [med, med],
                    color="white", linewidth=1.3,
                    solid_capstyle="round", zorder=4)

        # Per-donor pseudobulk dots (jittered)
        rng = np.random.default_rng(42)
        for i, (dvals, stage) in enumerate(zip(donor_data, STAGE_ORDER)):
            jitter = rng.uniform(-0.14, 0.14, size=len(dvals))
            ax.scatter(i + jitter, dvals,
                       s=5, color=STAGE_COLORS[stage],
                       alpha=0.65, linewidths=0, zorder=5)
            # Donor-level median crossbar
            dm = np.median(dvals)
            ax.plot([i - 0.16, i + 0.16], [dm, dm],
                    color=STAGE_COLORS[stage], linewidth=1.5,
                    solid_capstyle="butt", zorder=6)

        ax.set_xticks(range(len(STAGE_ORDER)))
        ax.set_xticklabels(
            [f"{s}\n(n={n_donors[j]})" for j, s in enumerate(STAGE_SHORT)],
            fontsize=6, fontweight="bold")
        ax.set_xlim(-0.6, len(STAGE_ORDER) - 0.4)

        ax.set_title(label, fontsize=7.5, fontweight="bold",
                     pad=4, color=color)
        ax.set_ylabel("Score (hepatocytes)" if ax is axes[0] else "",
                      fontsize=6.5)

        ax.yaxis.grid(True, linewidth=0.3, color="#E0E0E0", zorder=0)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Reference zero line
        ax.axhline(0, linewidth=0.5, color="#9E9E9E", linestyle="--",
                   zorder=1, alpha=0.7)

    # Legend for dot layer
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=STAGE_COLORS[s],
               markersize=4, label=s)
        for s in STAGE_ORDER
    ]
    legend_elements.append(
        Line2D([0], [0], color='#9E9E9E', linewidth=1.2, label='Donor median')
    )
    axes[-1].legend(handles=legend_elements, fontsize=5, frameon=False,
                    loc="upper right", bbox_to_anchor=(1.0, 1.0))

    fig.suptitle(
        "Bulk NMF program scores in hepatocytes (Healthy / ST / SH)",
        fontsize=7, y=1.01, fontweight="bold")

    plt.tight_layout(w_pad=1.0)
    outpath = OUT_DIR / "nmf_sc_hep_violin.pdf"
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {outpath}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    df = load_data()

    print("\n=== Panel A: cell-type composition ===")
    plot_panel_A(df)

    print("\n=== Panel C: hepatocyte-only violins ===")
    plot_panel_C(df)

    print("\nDone.")

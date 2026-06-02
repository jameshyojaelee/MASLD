#!/usr/bin/env python
# ============================================================================
# figS_hotspot_heatmap.py
#
# KEY MESSAGE: The canonical Hotspot per-CT gene x gene local-autocorrelation
# heatmap with hierarchical clustering. Modules show as colored ribbons on
# both axes; a second ribbon encodes per-module signed -log10(q_disease) so
# each module's disease relevance is visible in the same panoramic view.
#
# Outputs (PDF only):
#   figures/supplementary/figS_hotspot/figS_hotspot_heatmap.pdf  (main: global)
#   figures/supplementary/figS_hotspot/panels/heatmap_<ct>.pdf   (7 per-CT)
#   figures/supplementary/figS_hotspot/figS_hotspot_heatmap_atlas.pdf
#       (2x4 small-multiples atlas of all 7 cell-type heatmaps)
# ============================================================================
import os
import pickle
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42        # Illustrator-editable
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Helvetica"

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster

ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
HS = ROOT / "Analysis/SingleCell/results_gpu_v2/hotspot_modules"
OUT = ROOT / "figures/supplementary/figS_hotspot"
PANELS = OUT / "panels"
DATA = PANELS / "data"
for d in (OUT, PANELS, DATA): d.mkdir(parents=True, exist_ok=True)

# Liang palette
CONTROL_GRAY = "#9E9E9E"
DISEASE_RED  = "#C9265E"
DISEASE_DEEP = "#A01753"
DOWN_BLUE    = "#1565C0"

# Display names per CT
CT_DISPLAY = {
    "global": "Global atlas",
    "hepatocytes": "Hepatocytes",
    "macrophages": "Macrophages",
    "fibroblasts": "Fibroblasts",
    "endothelial_cells": "Endothelial cells",
    "cholangiocytes": "Cholangiocytes",
    "tcells": "T cells",
}
RUN_ORDER = ["global", "hepatocytes", "macrophages", "fibroblasts",
             "endothelial_cells", "cholangiocytes", "tcells"]


def load_run(ct):
    """Return modules (Series, gene -> module_id), local_corr (DataFrame)."""
    with open(HS / ct / "hotspot_obj.pkl", "rb") as h:
        d = pickle.load(h)
    return d["modules"], d["local_corr"]


def load_disease_beta():
    """phenotype_correlations.tsv → dict[(ct, module_int)] = (beta, q)."""
    df = pd.read_csv(HS / "phenotype_correlations.tsv", sep="\t")
    df = df[df["axis"] == "disease_stage"]
    df["module_int"] = df["module"].str.split("__").str[-1].astype(int)
    return {(r["cell_type"], r["module_int"]): (r["beta"], r["q"])
            for _, r in df.iterrows()}


DISEASE_BETA = load_disease_beta()


def module_palette(n):
    """Categorical palette for module IDs (mod 0 = unassigned = light gray)."""
    base = plt.get_cmap("tab20")(np.linspace(0, 1, 20))
    # Cycle through if more than 20 modules
    cols = [base[i % 20] for i in range(n)]
    return cols


def signed_logq_color(beta, q, vmin=-8, vmax=8):
    if pd.isna(beta) or pd.isna(q): return (0.95, 0.95, 0.95, 1.0)
    s = float(np.sign(beta)) * float(min(-np.log10(max(q, 1e-20)), 15))
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "div_liang", [DOWN_BLUE, "#FFFFFF", DISEASE_RED])
    return cmap(norm(s))


def plot_heatmap(ct, ax_main, ax_mod_top, ax_dis_top, ax_mod_left,
                 ax_dis_left, ax_cbar_main, ax_cbar_dis, title=None,
                 show_strip_labels=True):
    """Render one CT's local-corr heatmap with two side-bars.

    show_strip_labels: when True (default), label the left ribbons "Module"
    and "Disease β". On small panels (atlas grid) labels squeeze into the
    heatmap area; set False there and rely on the figure-level legend."""
    modules, lc = load_run(ct)

    # Reorder genes by hierarchical clustering on local_corr (1 - corr).
    # Use module ID as tiebreaker so module blocks stay contiguous.
    keep = modules > 0           # drop unassigned (module == 0)
    mods_kept = modules[keep]
    lc_kept = lc.loc[keep, keep]

    # Cluster within module then concatenate (preserves module blocks).
    # local_corr is signed Z (NOT Pearson on [-1,1]) — values can exceed |1|.
    # Convert similarity -> non-negative distance by `dmax - sim`.
    ordered_idx = []
    unique_mods = sorted(mods_kept.unique())
    for m in unique_mods:
        gs = mods_kept[mods_kept == m].index.tolist()
        if len(gs) > 1:
            sub = lc_kept.loc[gs, gs].values
            sim = sub.copy()
            np.fill_diagonal(sim, sim.max())   # ignore self in computing dmax
            dmax = sim.max()
            dist = dmax - sub
            np.fill_diagonal(dist, 0)
            # Ensure symmetric + non-negative (numerical noise can leave tiny negatives)
            dist = (dist + dist.T) / 2
            dist[dist < 0] = 0
            condensed = dist[np.triu_indices_from(dist, k=1)]
            Z = linkage(condensed, method="average")
            order = leaves_list(Z)
            ordered_idx.extend([gs[i] for i in order])
        else:
            ordered_idx.extend(gs)

    lc_ord = lc_kept.loc[ordered_idx, ordered_idx].values
    mods_ord = mods_kept.loc[ordered_idx].values

    # Cap colorscale at the 2-98 percentile for visual stability
    vmax = float(np.percentile(np.abs(lc_ord[lc_ord != 0]), 98))
    vmax = max(vmax, 4.0)
    cmap_main = mcolors.LinearSegmentedColormap.from_list(
        "lc_corr", [DOWN_BLUE, "#FFFFFF", DISEASE_RED])
    norm_main = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    im = ax_main.imshow(lc_ord, cmap=cmap_main, norm=norm_main,
                        aspect="auto", interpolation="nearest")
    ax_main.set_xticks([]); ax_main.set_yticks([])
    for sp in ax_main.spines.values(): sp.set_linewidth(0.3)

    # Module ribbon (left + top), with axis labels for the legend
    n = len(unique_mods)
    pal = module_palette(n)
    mod_to_col = {m: pal[i] for i, m in enumerate(unique_mods)}
    mod_strip = np.array([mcolors.to_rgb(mod_to_col[m]) for m in mods_ord])
    ax_mod_top.imshow(mod_strip[np.newaxis, :, :], aspect="auto")
    ax_mod_top.set_xticks([]); ax_mod_top.set_yticks([])
    for sp in ax_mod_top.spines.values(): sp.set_visible(False)
    ax_mod_left.imshow(mod_strip[:, np.newaxis, :], aspect="auto")
    ax_mod_left.set_xticks([]); ax_mod_left.set_yticks([])
    for sp in ax_mod_left.spines.values(): sp.set_visible(False)
    if show_strip_labels:
        ax_mod_left.set_ylabel("Module", fontsize=6, rotation=90, labelpad=2)

    # Disease β ribbon (per-gene = per-module value)
    dis_colors = []
    for m in mods_ord:
        beta, q = DISEASE_BETA.get((ct, int(m)), (np.nan, np.nan))
        dis_colors.append(signed_logq_color(beta, q))
    dis_strip = np.array(dis_colors)
    ax_dis_top.imshow(dis_strip[np.newaxis, :, :3], aspect="auto")
    ax_dis_top.set_xticks([]); ax_dis_top.set_yticks([])
    for sp in ax_dis_top.spines.values(): sp.set_visible(False)
    ax_dis_left.imshow(dis_strip[:, np.newaxis, :3], aspect="auto")
    ax_dis_left.set_xticks([]); ax_dis_left.set_yticks([])
    for sp in ax_dis_left.spines.values(): sp.set_visible(False)
    if show_strip_labels:
        ax_dis_left.set_ylabel("Disease β", fontsize=6, rotation=90, labelpad=2)

    # Module-block separator lines (gentle)
    cum = 0
    sep = []
    for m in unique_mods:
        cnt = (mods_ord == m).sum()
        cum += cnt
        sep.append(cum)
    for s in sep[:-1]:
        ax_main.axhline(s - 0.5, color="white", linewidth=0.15, alpha=0.7)
        ax_main.axvline(s - 0.5, color="white", linewidth=0.15, alpha=0.7)

    # Module-ID labels above the top ribbon with leader segments.
    # Show labels only for modules large enough to read at panel size;
    # smaller modules are listed in the side data CSV instead.
    cum = 0
    n_genes_total = len(mods_ord)
    # min-width threshold scales with figure: ~2% of x-axis or 15 genes, whichever larger
    min_genes_for_label = max(15, int(0.02 * n_genes_total))
    for m in unique_mods:
        cnt = (mods_ord == m).sum()
        center = cum + cnt / 2
        if cnt >= min_genes_for_label:
            # Bold black text above the ribbon (so legible regardless of strip color)
            ax_mod_top.annotate(str(int(m)),
                                xy=(center, 1.0), xycoords=("data", "axes fraction"),
                                xytext=(center, 1.85), textcoords=("data", "axes fraction"),
                                ha="center", va="bottom",
                                fontsize=6, fontweight="bold", color="black",
                                arrowprops=dict(arrowstyle="-", color="0.6",
                                                lw=0.25, shrinkA=0, shrinkB=0))
        cum += cnt

    # Main colorbar
    if ax_cbar_main is not None:
        cbm = plt.colorbar(im, cax=ax_cbar_main, orientation="horizontal")
        cbm.set_label("Local correlation Z", fontsize=5)
        cbm.ax.tick_params(labelsize=4, length=2, width=0.3)
        cbm.outline.set_linewidth(0.3)

    # Disease β colorbar
    if ax_cbar_dis is not None:
        # Build a separate scalar mappable for the disease scale
        norm_d = mcolors.TwoSlopeNorm(vmin=-8, vcenter=0, vmax=8)
        cmap_d = mcolors.LinearSegmentedColormap.from_list(
            "div_liang", [DOWN_BLUE, "#FFFFFF", DISEASE_RED])
        sm = plt.cm.ScalarMappable(norm=norm_d, cmap=cmap_d)
        cbd = plt.colorbar(sm, cax=ax_cbar_dis, orientation="horizontal")
        cbd.set_label("sign(β) × −log10(q_disease)", fontsize=5)
        cbd.ax.tick_params(labelsize=4, length=2, width=0.3)
        cbd.outline.set_linewidth(0.3)

    if title:
        ax_main.set_title(title, fontsize=7, fontweight="bold", pad=12)

    return n, len(ordered_idx), unique_mods, mods_ord


def render_single(ct, outpath, with_colorbars=True):
    """Single CT figure: gene x gene heatmap + 2 sidebars + colorbars.
    Extra top whitespace is reserved for offset module-ID labels."""
    fig = plt.figure(figsize=(5.4, 5.6))
    # GridSpec: [side ribbons] [heatmap] | colorbar row at bottom
    gs = gridspec.GridSpec(
        4, 4,
        width_ratios=[0.4, 0.4, 5, 0.05],
        height_ratios=[0.4, 0.4, 5, 0.6],
        wspace=0.05, hspace=0.05,
        left=0.08, right=0.98, top=0.82, bottom=0.07)

    ax_mod_top  = fig.add_subplot(gs[0, 2])
    ax_dis_top  = fig.add_subplot(gs[1, 2])
    ax_mod_left = fig.add_subplot(gs[2, 0])
    ax_dis_left = fig.add_subplot(gs[2, 1])
    ax_main     = fig.add_subplot(gs[2, 2])
    ax_cbar     = fig.add_subplot(gs[3, 1:3]) if with_colorbars else None
    ax_cbar_dis = fig.add_subplot(gs[3, 0:1]) if with_colorbars else None
    # Move disease cbar to a more sensible position (under the heatmap, right half)
    if with_colorbars:
        ax_cbar.set_position([0.4, 0.04, 0.30, 0.018])
        ax_cbar_dis.set_position([0.72, 0.04, 0.22, 0.018])

    # Reserve clean space above the top ribbons for offset labels (no axis frame)
    for sp in ax_mod_top.spines.values(): sp.set_visible(False)
    ax_mod_top.set_clip_on(False)

    n_mods, n_genes, uniq, mods_ord = plot_heatmap(
        ct, ax_main, ax_mod_top, ax_dis_top, ax_mod_left, ax_dis_left,
        ax_cbar, ax_cbar_dis, title=None)
    # Use suptitle so it sits above the offset-label band, not on the strips
    fig.suptitle(
        f"{CT_DISPLAY[ct]} — {n_mods} modules, {n_genes} genes",
        fontsize=8, fontweight="bold", y=0.96)
    fig.savefig(outpath, format="pdf", bbox_inches="tight")
    plt.close(fig)
    # Side data: ordered modules + disease beta
    rec = pd.DataFrame({
        "module": mods_ord,
        "disease_beta": [DISEASE_BETA.get((ct, int(m)), (np.nan, np.nan))[0]
                          for m in mods_ord],
        "disease_q": [DISEASE_BETA.get((ct, int(m)), (np.nan, np.nan))[1]
                       for m in mods_ord],
    })
    return rec


def render_atlas_grid(outpath):
    """Small-multiples 2x4 atlas of all 7 cell-type heatmaps.
    Module-ID labels are dropped at this size — see per-CT panels in panels/
    for legible module labels. Legend at top of the figure."""
    fig = plt.figure(figsize=(9.0, 5.6))
    # Generous hspace so the per-panel titles have headroom above the strips
    gs = gridspec.GridSpec(
        2, 4, wspace=0.20, hspace=0.85,
        left=0.04, right=0.99, top=0.88, bottom=0.08)
    for i, ct in enumerate(RUN_ORDER):
        r, c = i // 4, i % 4
        sub = gridspec.GridSpecFromSubplotSpec(
            3, 3,
            subplot_spec=gs[r, c],
            width_ratios=[0.12, 0.12, 3],
            height_ratios=[0.12, 0.12, 3],
            wspace=0.03, hspace=0.03)
        ax_mod_top  = fig.add_subplot(sub[0, 2])
        ax_dis_top  = fig.add_subplot(sub[1, 2])
        ax_mod_left = fig.add_subplot(sub[2, 0])
        ax_dis_left = fig.add_subplot(sub[2, 1])
        ax_main     = fig.add_subplot(sub[2, 2])
        n_mods, n_genes, _, _ = plot_heatmap(
            ct, ax_main, ax_mod_top, ax_dis_top,
            ax_mod_left, ax_dis_left, None, None, title=None,
            show_strip_labels=False)
        # Hide auto-generated module-ID label annotations (too small to read here)
        for child in list(ax_mod_top.texts):
            try: child.remove()
            except Exception: pass
        # Title anchored to the TOPMOST strip (ax_mod_top) so it sits above
        # the whole stack of ribbons + heatmap, never overlapping them.
        ax_mod_top.set_title(
            f"{CT_DISPLAY[ct]}\n{n_mods} modules · {n_genes} genes",
            fontsize=7, fontweight="bold", pad=4)
    # Top shared legend
    fig.text(0.04, 0.96,
             "Each panel: gene × gene local-autocorrelation Z. Top + left "
             "ribbons = module ID (cycling palette). Inner ribbons = signed "
             "−log10(q_disease) per module. See panels/heatmap_<ct>.pdf for "
             "labelled per-CT views.",
             fontsize=6, color="0.25", wrap=True)
    fig.savefig(outpath, format="pdf", bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Render everything
# ----------------------------------------------------------------------------
records = []
for ct in RUN_ORDER:
    out = PANELS / f"heatmap_{ct}.pdf"
    print(f"[{ct}] rendering -> {out}")
    rec = render_single(ct, out)
    rec["cell_type"] = ct
    records.append(rec)

pd.concat(records).to_csv(DATA / "heatmap_ordered_modules.csv", index=False)

# Main composite: global as the comprehensive view
print("rendering main composite (global)")
render_single("global", OUT / "figS_hotspot_heatmap.pdf")

# Atlas grid
print("rendering atlas grid")
render_atlas_grid(OUT / "figS_hotspot_heatmap_atlas.pdf")

print("done")

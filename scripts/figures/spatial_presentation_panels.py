#!/usr/bin/env python3
"""
spatial_presentation_panels.py — Presentation-quality spatial transcriptomics panels.

Generates individual high-resolution panels AND a combined composite figure:
  1. Cell-type deconvolution spatial maps (cell2location)
  2. Liver zonation gradient (periportal ↔ pericentral)
  3. Disease trajectory pseudotime on tissue
  4. Spatial domains with cell-type composition
  5. Spatially variable gene heatmap (Moran's I, disease-emergent/lost/stable)
  6. COMMOT ligand-receptor signaling flow
  7. Gene expression spatial maps (MASLD-relevant genes)
  8. SVG volcano: Moran's I healthy vs MASLD
  9. Zonation-stratified cell-type composition
 10. Disease communication shift (waterfall)

Outputs to: Analysis/Spatial/results/presentation_panels/
  - Individual panels: panel_01_*.png/pdf through panel_10_*.png/pdf
  - Combined composite: spatial_composite.png/pdf

SLURM: --partition=cpu --cpus=8 --mem=64G --time=4:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import Normalize, LinearSegmentedColormap, TwoSlopeNorm
from matplotlib import cm
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyArrowPatch
from scipy.stats import spearmanr, mannwhitneyu
import warnings
warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SPATIAL_DIR = PROJECT_ROOT / "Analysis" / "Spatial"
RESULTS_DIR = SPATIAL_DIR / "results"
OUT_DIR = RESULTS_DIR / "presentation_panels"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Presentation-quality rcParams ──────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 12,
    "axes.titlesize": 16,
    "axes.labelsize": 14,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.15,
    "axes.linewidth": 1.0,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 5,
    "ytick.major.size": 5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# ── Color palettes ─────────────────────────────────────────────────────────────
CONDITION_COLORS = {"Healthy": "#9E9E9E", "Steatotic": "#C9265E"}
CELL_TYPE_COLORS = {
    "Hepatocytes": "#1B4F72", "Fibroblasts": "#E67E22",
    "Endothelial cells": "#27AE60", "Macrophages": "#C0392B",
    "Cholangiocytes": "#8E44AD", "T cells": "#2980B9",
    "B cells": "#F39C12", "Neutrophils": "#7F8C8D",
    "Mono+mono derived cells": "#16A085", "Plasma cells": "#E74C3C",
    "Circulating NK/NKT": "#D4AC0D", "Resident NK": "#AF7AC5",
    "cDC1s": "#48C9B0", "cDC2s": "#EC7063", "pDCs": "#85929E",
    "Basophils": "#D5DBDB",
}
ZONATION_CMAP = LinearSegmentedColormap.from_list(
    "zonation", ["#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"]
)
SVG_CATEGORY_COLORS = {
    "stable": "#7F8C8D",
    "disease_emergent_SVG": "#C0392B",
    "disease_lost_SVG": "#2E86AB",
}
DOMAIN_COLORS = ["#E74C3C", "#3498DB", "#2ECC71", "#F39C12", "#9B59B6",
                 "#1ABC9C", "#E67E22", "#34495E"]

C2L_PREFIX = "c2l_q05cell_abundance_w_sf_means_per_cluster_mu_fg_"


def save_panel(fig, name, formats=("png", "pdf")):
    """Save figure panel in multiple formats."""
    for fmt in formats:
        fig.savefig(OUT_DIR / f"{name}.{fmt}", format=fmt)
    plt.close(fig)
    print(f"  ✓ Saved: {name}")


def pick_representative_sample(adata, condition):
    """Pick the sample with the most spots for a condition."""
    mask = adata.obs["condition"] == condition
    return adata.obs.loc[mask, "sample_id"].value_counts().index[0]


def add_scalebar(ax, coords, length_um=500, pixel_per_um=1.0, fontsize=9):
    """Add a scale bar to spatial plot."""
    xmin, xmax = coords[:, 0].min(), coords[:, 0].max()
    ymin, ymax = coords[:, 1].min(), coords[:, 1].max()
    bar_len = (xmax - xmin) * 0.15
    x0 = xmin + (xmax - xmin) * 0.05
    y0 = ymax - (ymax - ymin) * 0.05
    ax.plot([x0, x0 + bar_len], [y0, y0], "k-", lw=2.5)
    ax.text(x0 + bar_len / 2, y0 - (ymax - ymin) * 0.03,
            f"{length_um} µm", ha="center", fontsize=fontsize, fontweight="bold")


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 1: Cell-type deconvolution spatial maps
# ═══════════════════════════════════════════════════════════════════════════════
def panel_01_deconvolution(adata):
    """cell2location abundance maps — key cell types, Healthy vs Steatotic.

    Uses per-spot proportions (fraction of total cell abundance) instead of
    raw abundance to ensure fair visual comparison across samples with
    different spot counts. Picks size-matched samples and shares the
    colorbar range across conditions for each cell type.
    """
    print("[Panel 1] Cell-type deconvolution spatial maps")
    cell_types = ["Hepatocytes", "Macrophages", "Fibroblasts",
                  "Endothelial cells", "Cholangiocytes", "T cells"]
    c2l_cols = {ct: f"{C2L_PREFIX}{ct}" for ct in cell_types
                if f"{C2L_PREFIX}{ct}" in adata.obs.columns}
    if not c2l_cols:
        print("  SKIP: No cell2location columns found")
        return

    conditions = ["Healthy", "Steatotic"]

    # Pick size-matched samples (closest spot counts) instead of largest
    def pick_balanced_sample(adata, cond, target_n=None):
        mask = adata.obs["condition"] == cond
        counts = adata.obs.loc[mask, "sample_id"].value_counts()
        counts = counts[counts > 0]
        if target_n is not None:
            return counts.index[(counts - target_n).abs().argmin()]
        return counts.index[0]

    # Pick steatotic first (fewer spots), then match healthy
    steatotic_rep = pick_balanced_sample(adata, "Steatotic")
    steatotic_n = (adata.obs["condition"] == "Steatotic").sum() and \
                  ((adata.obs["condition"] == "Steatotic") & (adata.obs["sample_id"] == steatotic_rep)).sum()
    healthy_rep = pick_balanced_sample(adata, "Healthy", target_n=steatotic_n)
    reps = {"Healthy": healthy_rep, "Steatotic": steatotic_rep}
    for c, r in reps.items():
        n = ((adata.obs["condition"] == c) & (adata.obs["sample_id"] == r)).sum()
        print(f"  {c}: sample {r} ({n} spots)")

    # Compute per-spot proportions: fraction of total cell abundance at each spot
    all_c2l = list(c2l_cols.values())
    prop_df = adata.obs[all_c2l].copy()
    row_totals = prop_df.sum(axis=1)
    row_totals = row_totals.replace(0, 1)  # avoid division by zero
    prop_df = prop_df.div(row_totals, axis=0)

    # Pre-compute shared colorbar range per cell type (across both conditions)
    shared_vmax = {}
    for ct, col in c2l_cols.items():
        vals_all = []
        for cond in conditions:
            mask = (adata.obs["condition"] == cond) & (adata.obs["sample_id"] == reps[cond])
            vals_all.append(prop_df.loc[mask, col].values)
        combined = np.concatenate(vals_all)
        shared_vmax[ct] = np.percentile(combined, 97) if combined.max() > 0 else 1

    fig, axes = plt.subplots(len(conditions), len(c2l_cols), figsize=(16, 6.5))
    fig.subplots_adjust(wspace=0.05, hspace=0.25)
    fig.suptitle("Cell-Type Deconvolution (cell2location)", fontsize=22, fontweight="bold", y=1.01)

    for i, cond in enumerate(conditions):
        mask = (adata.obs["condition"] == cond) & (adata.obs["sample_id"] == reps[cond])
        coords = adata[mask].obsm["spatial"]
        for j, (ct, col) in enumerate(c2l_cols.items()):
            ax = axes[i, j]
            vals = prop_df.loc[mask, col].values
            vmax = shared_vmax[ct]
            sc_map = ax.scatter(
                coords[:, 0], coords[:, 1], c=vals, cmap="magma",
                s=6, vmin=0, vmax=vmax, edgecolors="none", rasterized=True
            )
            ax.set_aspect("equal")
            ax.axis("off")
            if i == 0:
                ax.set_title(ct, fontsize=16, fontweight="bold",
                             color=CELL_TYPE_COLORS.get(ct, "black"), pad=4)
            if j == 0:
                ax.text(-0.05, 0.5, cond, transform=ax.transAxes, fontsize=16,
                        fontweight="bold", va="center", ha="right", rotation=90,
                        color=CONDITION_COLORS[cond])
            cb = fig.colorbar(sc_map, ax=ax, shrink=0.55, pad=0.01, aspect=10)
            cb.ax.tick_params(labelsize=9)
            if i == 1 and j == len(c2l_cols) - 1:
                cb.set_label("Proportion", fontsize=10)

    save_panel(fig, "panel_01_cell_type_deconvolution", formats=("pdf",))


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 2: Liver zonation gradient
# ═══════════════════════════════════════════════════════════════════════════════
def panel_02_zonation(adata):
    """Periportal ↔ Pericentral zonation gradient on tissue sections."""
    print("[Panel 2] Liver zonation gradient")
    if "zonation_score" not in adata.obs.columns:
        print("  SKIP: No zonation_score")
        return

    conditions = ["Healthy", "Steatotic"]
    reps = {c: pick_representative_sample(adata, c) for c in conditions}

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle("Hepatic Zonation Gradient", fontsize=20, fontweight="bold", y=1.02)

    for i, cond in enumerate(conditions):
        ax = axes[i]
        mask = (adata.obs["condition"] == cond) & (adata.obs["sample_id"] == reps[cond])
        coords = adata[mask].obsm["spatial"]
        vals = adata.obs.loc[mask, "zonation_score"].values
        vabs = max(abs(np.nanpercentile(vals, 2)), abs(np.nanpercentile(vals, 98)))

        sc_map = ax.scatter(
            coords[:, 0], coords[:, 1], c=vals, cmap=ZONATION_CMAP,
            s=12, vmin=-vabs, vmax=vabs, edgecolors="none", rasterized=True
        )
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(cond, fontsize=16, fontweight="bold", color=CONDITION_COLORS[cond])
        cb = fig.colorbar(sc_map, ax=ax, shrink=0.7, pad=0.03, aspect=20)
        cb.set_label("Periportal ← → Pericentral", fontsize=11)
        cb.ax.tick_params(labelsize=10)

    fig.tight_layout()
    save_panel(fig, "panel_02_zonation_gradient")


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 3: Disease trajectory pseudotime
# ═══════════════════════════════════════════════════════════════════════════════
def panel_03_trajectory(adata_traj):
    """Disease pseudotime on steatotic tissue sections."""
    print("[Panel 3] Disease trajectory pseudotime")
    if "dpt_pseudotime" not in adata_traj.obs.columns:
        print("  SKIP: No dpt_pseudotime")
        return

    fig, axes = plt.subplots(1, 3, figsize=(22, 7),
                              gridspec_kw={"width_ratios": [1, 1, 1.3]})
    fig.suptitle("Disease Trajectory — Steatotic Liver", fontsize=20, fontweight="bold", y=1.02)

    # Pick two steatotic samples
    samples = adata_traj.obs["sample_id"].value_counts().index[:2].tolist()

    for i, samp in enumerate(samples):
        ax = axes[i]
        mask = adata_traj.obs["sample_id"] == samp
        coords = adata_traj[mask].obsm["spatial"]
        vals = adata_traj.obs.loc[mask, "dpt_pseudotime"].values

        sc_map = ax.scatter(
            coords[:, 0], coords[:, 1], c=vals, cmap="viridis",
            s=12, vmin=0, vmax=1, edgecolors="none", rasterized=True
        )
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(f"Sample {samp}", fontsize=14, fontweight="bold")
        cb = fig.colorbar(sc_map, ax=ax, shrink=0.7, pad=0.03, aspect=20)
        cb.set_label("Pseudotime", fontsize=11)

    # Gene dynamics along trajectory
    ax3 = axes[2]
    gene_dyn_path = RESULTS_DIR / "trajectory" / "gene_dynamics_along_trajectory.csv"
    if gene_dyn_path.exists():
        dyn = pd.read_csv(gene_dyn_path)
        genes_to_plot = ["FASN", "SCD", "CYP2E1", "COL1A1", "ALB", "HSD17B13",
                         "GLUL", "HAL"]
        genes_avail = [g for g in genes_to_plot if g in dyn["gene"].values][:6]

        cmap_genes = plt.cm.get_cmap("tab10", len(genes_avail))
        for gi, gene in enumerate(genes_avail):
            gdata = dyn[dyn["gene"] == gene].sort_values("bin")
            ax3.plot(gdata["mean_dpt"], gdata["mean_expr"], "-o",
                     color=cmap_genes(gi), lw=2.5, ms=5, label=gene)
            ax3.fill_between(gdata["mean_dpt"],
                             gdata["mean_expr"] - gdata["std_expr"],
                             gdata["mean_expr"] + gdata["std_expr"],
                             color=cmap_genes(gi), alpha=0.12)
        ax3.set_xlabel("Disease Pseudotime", fontsize=14)
        ax3.set_ylabel("Mean Expression", fontsize=14)
        ax3.set_title("Gene Dynamics Along Trajectory", fontsize=14, fontweight="bold")
        ax3.legend(fontsize=10, frameon=True, framealpha=0.9, edgecolor="gray",
                   loc="upper left", ncol=2)
        ax3.grid(True, alpha=0.2, ls="--")

    fig.tight_layout()
    save_panel(fig, "panel_03_disease_trajectory")


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 4: Spatial domains + composition
# ═══════════════════════════════════════════════════════════════════════════════
def panel_04_domains(adata_dom):
    """Spatial domain maps with cell-type composition bars."""
    print("[Panel 4] Spatial domains + composition")
    if "spatial_domain" not in adata_dom.obs.columns:
        print("  SKIP: No spatial_domain")
        return

    conditions = ["Healthy", "Steatotic"]
    reps = {c: pick_representative_sample(adata_dom, c) for c in conditions}
    n_domains = adata_dom.obs["spatial_domain"].nunique()
    domain_cmap = mcolors.ListedColormap(DOMAIN_COLORS[:n_domains])

    fig = plt.figure(figsize=(22, 8))
    gs = GridSpec(1, 3, width_ratios=[1, 1, 1.2], wspace=0.15)
    fig.suptitle("Spatial Domains & Cell-Type Composition", fontsize=20, fontweight="bold", y=1.02)

    for i, cond in enumerate(conditions):
        ax = fig.add_subplot(gs[0, i])
        mask = (adata_dom.obs["condition"] == cond) & (adata_dom.obs["sample_id"] == reps[cond])
        coords = adata_dom[mask].obsm["spatial"]
        domains = adata_dom.obs.loc[mask, "spatial_domain"].astype(int).values

        sc_map = ax.scatter(
            coords[:, 0], coords[:, 1], c=domains, cmap=domain_cmap,
            s=12, edgecolors="none", rasterized=True,
            vmin=-0.5, vmax=n_domains - 0.5
        )
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(cond, fontsize=16, fontweight="bold", color=CONDITION_COLORS[cond])

    # Composition stacked bars
    comp_path = RESULTS_DIR / "domains" / "domain_composition.csv"
    if comp_path.exists():
        ax_comp = fig.add_subplot(gs[0, 2])
        df = pd.read_csv(comp_path)
        if "spatial_domain" in df.columns:
            ct_cols = [c for c in df.columns if c != "spatial_domain" and c != "Unnamed: 0"]
            ct_cols_clean = []
            for c in ct_cols:
                clean = c.replace("q05cell_abundance_w_sf_means_per_cluster_mu_fg_", "")
                clean = clean.replace("c2l_q05cell_abundance_w_sf_means_per_cluster_mu_fg_", "")
                ct_cols_clean.append(clean)
            df_plot = df[["spatial_domain"] + ct_cols].copy()
            df_plot.columns = ["Domain"] + ct_cols_clean
            df_plot = df_plot.set_index("Domain")

            # Normalize rows
            row_sums = df_plot.sum(axis=1)
            df_plot = df_plot.div(row_sums, axis=0)

            bottom = np.zeros(len(df_plot))
            x = np.arange(len(df_plot))
            for ct in ct_cols_clean:
                color = CELL_TYPE_COLORS.get(ct, "#CCCCCC")
                ax_comp.barh(x, df_plot[ct].values, left=bottom, color=color,
                             edgecolor="white", linewidth=0.5, label=ct, height=0.7)
                bottom += df_plot[ct].values

            ax_comp.set_yticks(x)
            ax_comp.set_yticklabels([f"Domain {i}" for i in df_plot.index], fontsize=11)
            ax_comp.set_xlabel("Proportion", fontsize=13)
            ax_comp.set_title("Cell-Type Composition", fontsize=14, fontweight="bold")
            ax_comp.set_xlim(0, 1)
            ax_comp.invert_yaxis()

            # Legend outside
            handles, labels = ax_comp.get_legend_handles_labels()
            # Filter to top 8 cell types by mean proportion
            mean_props = {ct: df_plot[ct].mean() for ct in ct_cols_clean}
            top_cts = sorted(mean_props, key=mean_props.get, reverse=True)[:8]
            h_filt = [h for h, l in zip(handles, labels) if l in top_cts]
            l_filt = [l for l in labels if l in top_cts]
            ax_comp.legend(h_filt, l_filt, fontsize=8, loc="lower right",
                           frameon=True, framealpha=0.9, ncol=2)

    fig.tight_layout()
    save_panel(fig, "panel_04_spatial_domains")


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 5: SVG heatmap — Moran's I comparison
# ═══════════════════════════════════════════════════════════════════════════════
def panel_05_svg_heatmap():
    """Spatially variable gene heatmap with disease categorization."""
    print("[Panel 5] SVG Moran's I heatmap")
    svg_path = RESULTS_DIR / "svg" / "differential_svgs.csv"
    if not svg_path.exists():
        print("  SKIP: No SVG data")
        return

    svg = pd.read_csv(svg_path).rename(columns={"Unnamed: 0": "gene"})

    fig, axes = plt.subplots(1, 3, figsize=(15, 6),
                              gridspec_kw={"width_ratios": [1.5, 1, 1]})
    fig.subplots_adjust(wspace=0.45)
    fig.suptitle("Spatially Variable Genes in MASLD", fontsize=22, fontweight="bold", y=1.03)

    # (a) Scatter: Moran's I healthy vs steatotic
    ax = axes[0]
    for cat, color in SVG_CATEGORY_COLORS.items():
        mask = svg["category"] == cat
        n = mask.sum()
        label = cat.replace("_", " ").replace("SVG", "").strip().title()
        ax.scatter(svg.loc[mask, "morans_I_healthy"], svg.loc[mask, "morans_I_masld"],
                   c=color, s=12 if cat == "stable" else 30, alpha=0.4 if cat == "stable" else 0.8,
                   edgecolors="none", label=f"{label} (n={n})", rasterized=True, zorder=2 if cat == "stable" else 3)

    # Highlight key genes
    highlight_genes = ["SCD", "COL1A1", "CYP7A1", "HP", "SERPINE1", "CCL19", "CYP2E1", "GLUL"]
    for g in highlight_genes:
        row = svg[svg["gene"] == g]
        if not row.empty:
            x, y = row["morans_I_healthy"].values[0], row["morans_I_masld"].values[0]
            cat = row["category"].values[0]
            ax.annotate(g, (x, y), fontsize=12, fontweight="bold",
                        textcoords="offset points", xytext=(8, 5),
                        arrowprops=dict(arrowstyle="->", color="black", lw=1.0),
                        color=SVG_CATEGORY_COLORS.get(cat, "black"), zorder=5)

    lim = max(svg["morans_I_healthy"].max(), svg["morans_I_masld"].max()) * 1.1
    ax.plot([0, lim], [0, lim], "--", color="gray", alpha=0.5, lw=1)
    ax.set_xlabel("Moran's I (Healthy)", fontsize=16)
    ax.set_ylabel("Moran's I (Steatotic)", fontsize=16)
    ax.set_title("Spatial Autocorrelation Shift", fontsize=16, fontweight="bold")
    ax.legend(fontsize=11, frameon=True, framealpha=0.9, loc="upper left")
    ax.tick_params(axis="both", labelsize=12)
    ax.grid(True, alpha=0.15, ls="--")

    # (b) Top disease-emergent genes — horizontal bar
    ax2 = axes[1]
    emergent = svg[svg["category"] == "disease_emergent_SVG"].nlargest(15, "morans_I_masld")
    y_pos = np.arange(len(emergent))
    bars = ax2.barh(y_pos, emergent["morans_I_masld"].values,
                    color="#C0392B", alpha=0.85, edgecolor="white", height=0.7)
    ax2.barh(y_pos, emergent["morans_I_healthy"].values,
             color="#2E86AB", alpha=0.5, edgecolor="white", height=0.7)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(emergent["gene"].values, fontsize=12, fontweight="bold")
    ax2.set_xlabel("Moran's I", fontsize=15)
    ax2.set_title("Disease-Emergent SVGs", fontsize=16, fontweight="bold", color="#C0392B")
    ax2.invert_yaxis()
    ax2.legend(["Steatotic", "Healthy"], fontsize=11, loc="lower right")
    ax2.tick_params(axis="x", labelsize=12)

    # (c) Top disease-lost genes
    ax3 = axes[2]
    lost = svg[svg["category"] == "disease_lost_SVG"].nlargest(15, "morans_I_healthy")
    y_pos2 = np.arange(len(lost))
    ax3.barh(y_pos2, lost["morans_I_healthy"].values,
             color="#2E86AB", alpha=0.85, edgecolor="white", height=0.7)
    ax3.barh(y_pos2, lost["morans_I_masld"].values,
             color="#C0392B", alpha=0.5, edgecolor="white", height=0.7)
    ax3.set_yticks(y_pos2)
    ax3.set_yticklabels(lost["gene"].values, fontsize=12, fontweight="bold")
    ax3.set_xlabel("Moran's I", fontsize=15)
    ax3.set_title("Disease-Lost SVGs", fontsize=16, fontweight="bold", color="#2E86AB")
    ax3.invert_yaxis()
    ax3.legend(["Healthy", "Steatotic"], fontsize=11, loc="lower right")
    ax3.tick_params(axis="x", labelsize=12)

    save_panel(fig, "panel_05_svg_analysis", formats=("pdf",))


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 6: COMMOT signaling flow on tissue
# ═══════════════════════════════════════════════════════════════════════════════
def panel_06_commot_signaling():
    """COMMOT ligand-receptor signaling spatial maps + differential waterfall."""
    print("[Panel 6] COMMOT signaling")

    # Load COMMOT objects
    h_path = RESULTS_DIR / "commot" / "adata_commot_Healthy.h5ad"
    s_path = RESULTS_DIR / "commot" / "adata_commot_Steatotic.h5ad"
    if not h_path.exists() or not s_path.exists():
        print("  SKIP: No COMMOT data")
        return

    adata_h = sc.read_h5ad(h_path)
    adata_s = sc.read_h5ad(s_path)

    # Key MASLD-relevant pathways
    pathways = ["s-TGFB1-TGFBR1_TGFBR2", "s-HGF-MET", "s-VEGFA-FLT1",
                "s-C3-C3AR1", "s-IGF1-IGF1R", "s-CSF1-CSF1R"]
    available = [p for p in pathways if p in adata_h.obsm["commot-CellChat-sum-sender"].columns]

    if not available:
        print("  SKIP: No matching COMMOT pathways")
        return

    n_paths = min(len(available), 4)
    fig, axes = plt.subplots(2, n_paths, figsize=(6 * n_paths, 12))
    fig.suptitle("Spatial Signaling Patterns (COMMOT)", fontsize=20, fontweight="bold", y=1.02)

    for j, pw in enumerate(available[:n_paths]):
        pw_label = pw.replace("s-", "").replace("_", "/")

        # Healthy
        ax_h = axes[0, j] if n_paths > 1 else axes[0]
        coords_h = adata_h.obsm["spatial"]
        vals_h = adata_h.obsm["commot-CellChat-sum-sender"][pw].values
        vmax = np.percentile(np.concatenate([vals_h,
                adata_s.obsm["commot-CellChat-sum-sender"][pw].values]), 97)
        vmax = max(vmax, 1e-6)

        sc1 = ax_h.scatter(coords_h[:, 0], coords_h[:, 1], c=vals_h,
                           cmap="YlOrRd", s=8, vmin=0, vmax=vmax,
                           edgecolors="none", rasterized=True)
        ax_h.set_aspect("equal")
        ax_h.axis("off")
        ax_h.set_title(f"{pw_label}\nHealthy", fontsize=12, fontweight="bold",
                       color=CONDITION_COLORS["Healthy"])

        # Steatotic
        ax_s = axes[1, j] if n_paths > 1 else axes[1]
        coords_s = adata_s.obsm["spatial"]
        vals_s = adata_s.obsm["commot-CellChat-sum-sender"][pw].values

        sc2 = ax_s.scatter(coords_s[:, 0], coords_s[:, 1], c=vals_s,
                           cmap="YlOrRd", s=8, vmin=0, vmax=vmax,
                           edgecolors="none", rasterized=True)
        ax_s.set_aspect("equal")
        ax_s.axis("off")
        ax_s.set_title(f"Steatotic", fontsize=12, fontweight="bold",
                       color=CONDITION_COLORS["Steatotic"])

        fig.colorbar(sc2, ax=[ax_h, ax_s], shrink=0.4, pad=0.02, aspect=15,
                     label="Sender strength")

    fig.tight_layout()
    save_panel(fig, "panel_06_commot_signaling")


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 7: Gene expression spatial maps
# ═══════════════════════════════════════════════════════════════════════════════
def panel_07_gene_expression(adata):
    """Spatial gene expression maps for MASLD-relevant genes."""
    print("[Panel 7] Gene expression spatial maps")

    genes = ["CYP2E1", "GLUL", "ALB", "SCD", "COL1A1", "HAL"]
    available = [g for g in genes if g in adata.var_names]
    if not available:
        print("  SKIP: No MASLD genes found in var_names")
        return

    conditions = ["Healthy", "Steatotic"]
    reps = {c: pick_representative_sample(adata, c) for c in conditions}
    n_genes = len(available)

    fig, axes = plt.subplots(len(conditions), n_genes, figsize=(5 * n_genes, 10))
    fig.suptitle("Spatial Gene Expression — MASLD-Relevant Genes", fontsize=20,
                 fontweight="bold", y=1.02)

    for i, cond in enumerate(conditions):
        mask = (adata.obs["condition"] == cond) & (adata.obs["sample_id"] == reps[cond])
        coords = adata[mask].obsm["spatial"]
        for j, gene in enumerate(available):
            ax = axes[i, j] if len(conditions) > 1 else axes[j]
            # Get expression from X matrix
            gene_idx = list(adata.var_names).index(gene)
            expr = np.asarray(adata[mask].X[:, gene_idx].todense()).flatten() \
                if hasattr(adata[mask].X, "todense") else adata[mask].X[:, gene_idx].flatten()

            vmax = np.percentile(expr, 97) if expr.max() > 0 else 1
            sc_map = ax.scatter(
                coords[:, 0], coords[:, 1], c=expr, cmap="Reds",
                s=10, vmin=0, vmax=vmax, edgecolors="none", rasterized=True
            )
            ax.set_aspect("equal")
            ax.axis("off")
            if i == 0:
                ax.set_title(gene, fontsize=16, fontweight="bold", fontstyle="italic")
            if j == 0:
                ax.text(-0.08, 0.5, cond, transform=ax.transAxes, fontsize=14,
                        fontweight="bold", va="center", ha="right", rotation=90,
                        color=CONDITION_COLORS[cond])
            fig.colorbar(sc_map, ax=ax, shrink=0.6, pad=0.02, aspect=12)

    fig.tight_layout()
    save_panel(fig, "panel_07_gene_expression_maps")


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 8: Hepatocyte-intrinsic validation scatter
# ═══════════════════════════════════════════════════════════════════════════════
def panel_08_hep_validation():
    """Bulk deconvolution attribution vs spatial hepatocyte FC."""
    print("[Panel 8] Hepatocyte-intrinsic validation")
    hep_path = RESULTS_DIR / "cell2location" / "hep_intrinsic_validation.csv"
    if not hep_path.exists():
        print("  SKIP: No hep validation data")
        return

    df = pd.read_csv(hep_path)
    x = df["attribution_raw"]
    y = df["spatial_fc"]
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid].values, y[valid].values
    validated = df.loc[valid, "validated"].values if "validated" in df.columns else np.zeros(len(x), dtype=bool)

    fig, ax = plt.subplots(figsize=(9, 8))

    # Background non-validated
    ax.scatter(x[~validated], y[~validated], s=40, alpha=0.35, c="#BDC3C7",
               edgecolors="white", linewidths=0.3, rasterized=True, zorder=2,
               label=f"Not validated (n={int((~validated).sum())})")
    # Validated
    ax.scatter(x[validated], y[validated], s=55, alpha=0.7, c="#2E86AB",
               edgecolors="white", linewidths=0.5, rasterized=True, zorder=3,
               label=f"Validated (n={int(validated.sum())})")

    # Regression line
    z = np.polyfit(x, y, 1)
    xline = np.linspace(x.min(), x.max(), 100)
    ax.plot(xline, np.polyval(z, xline), "--", c="#D64933", lw=2, alpha=0.8, zorder=4)

    rho, p = spearmanr(x, y)
    ax.axhline(1.0, ls=":", lw=1, c="gray", alpha=0.5)
    ax.set_xlabel("Bulk Deconvolution Attribution Score", fontsize=14)
    ax.set_ylabel("Spatial Hepatocyte Fold Change", fontsize=14)
    ax.set_title("Hepatocyte-Intrinsic Gene Validation\n(Bulk vs Spatial)", fontsize=16, fontweight="bold")

    # Stats box
    ax.annotate(f"Spearman ρ = {rho:.3f}\np = {p:.1e}\nn = {len(x)} genes",
                xy=(0.04, 0.96), xycoords="axes fraction", fontsize=13, va="top",
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="#2E86AB",
                          alpha=0.9, lw=1.5))
    ax.legend(fontsize=11, loc="lower right", frameon=True, framealpha=0.9,
              edgecolor="gray", markerscale=1.2)
    ax.grid(True, alpha=0.15, ls="--")

    # Label top validated genes
    top_genes = df.loc[valid & df["validated"]].nlargest(5, "spatial_fc")
    for _, row in top_genes.iterrows():
        ax.annotate(row["gene"], (row["attribution_raw"], row["spatial_fc"]),
                    fontsize=9, fontweight="bold", fontstyle="italic",
                    textcoords="offset points", xytext=(8, 5),
                    arrowprops=dict(arrowstyle="->", color="black", lw=0.8))

    fig.tight_layout()
    save_panel(fig, "panel_08_hep_validation")


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 9: Zonation-stratified cell-type proportions
# ═══════════════════════════════════════════════════════════════════════════════
def panel_09_zonation_celltypes(adata):
    """Cell-type proportions across zonation bins."""
    print("[Panel 9] Zonation-stratified cell types")
    if "zonation_bin" not in adata.obs.columns:
        print("  SKIP: No zonation_bin")
        return

    c2l_cols = [c for c in adata.obs.columns if c.startswith(C2L_PREFIX)]
    if not c2l_cols:
        print("  SKIP: No c2l columns")
        return

    # Key cell types for display
    display_cts = ["Hepatocytes", "Macrophages", "Fibroblasts",
                   "Endothelial cells", "Cholangiocytes", "T cells"]

    fig, axes = plt.subplots(1, 2, figsize=(18, 8))
    fig.suptitle("Cell-Type Abundance Across Zonation", fontsize=20, fontweight="bold", y=1.02)

    bin_order = ["PP1", "PP2", "Mid", "PC2", "PC1"]

    for ci, cond in enumerate(["Healthy", "Steatotic"]):
        ax = axes[ci]
        mask = adata.obs["condition"] == cond
        data_cond = adata.obs.loc[mask]

        for ct_name in display_cts:
            col = f"{C2L_PREFIX}{ct_name}"
            if col not in data_cond.columns:
                continue
            means = []
            sems = []
            for b in bin_order:
                vals = data_cond.loc[data_cond["zonation_bin"] == b, col].values
                means.append(np.mean(vals))
                sems.append(np.std(vals) / np.sqrt(len(vals)) if len(vals) > 1 else 0)

            color = CELL_TYPE_COLORS.get(ct_name, "#CCCCCC")
            ax.plot(range(len(bin_order)), means, "-o", color=color, lw=2.5, ms=7,
                    label=ct_name, zorder=3)
            ax.fill_between(range(len(bin_order)),
                            [m - s for m, s in zip(means, sems)],
                            [m + s for m, s in zip(means, sems)],
                            color=color, alpha=0.12)

        ax.set_xticks(range(len(bin_order)))
        ax.set_xticklabels(bin_order, fontsize=12)
        ax.set_xlabel("Zonation (Periportal → Pericentral)", fontsize=13)
        ax.set_ylabel("Mean Abundance (cell2location)", fontsize=13)
        ax.set_title(cond, fontsize=16, fontweight="bold", color=CONDITION_COLORS[cond])
        ax.legend(fontsize=9, frameon=True, framealpha=0.9, loc="upper right", ncol=2)
        ax.grid(True, alpha=0.15, ls="--")

    fig.tight_layout()
    save_panel(fig, "panel_09_zonation_celltypes")


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 10: Communication shift waterfall
# ═══════════════════════════════════════════════════════════════════════════════
def panel_10_communication_waterfall():
    """Waterfall plot of differential communication (Steatotic vs Healthy)."""
    print("[Panel 10] Communication shift waterfall")
    comm_path = RESULTS_DIR / "commot" / "differential_communication.csv"
    if not comm_path.exists():
        print("  SKIP: No differential communication data")
        return

    df = pd.read_csv(comm_path)
    if "total_communication_fc" not in df.columns:
        print("  SKIP: No total_communication_fc")
        return

    # Filter to pathways with meaningful signal
    df = df.dropna(subset=["total_communication_fc"])
    df = df[df["total_communication_fc"] > 0]
    df["log2FC"] = np.log2(df["total_communication_fc"])
    df = df.sort_values("log2FC")

    # Take top and bottom pathways
    top_up = df.nlargest(15, "log2FC")
    top_down = df.nsmallest(15, "log2FC")
    plot_df = pd.concat([top_down, top_up]).drop_duplicates()
    plot_df = plot_df.sort_values("log2FC")

    fig, ax = plt.subplots(figsize=(8, 9))
    y_pos = np.arange(len(plot_df))
    colors = ["#C0392B" if v > 0 else "#2E86AB" for v in plot_df["log2FC"]]

    ax.barh(y_pos, plot_df["log2FC"], color=colors, edgecolor="white",
            linewidth=0.5, height=0.7, alpha=0.85)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df["pathway"].values, fontsize=13, fontweight="bold")
    ax.axvline(0, color="black", lw=1)
    ax.set_xlabel("log₂ Fold Change (Steatotic / Healthy)", fontsize=16)
    ax.set_title("Differential Spatial Communication\n(COMMOT Ligand-Receptor Pathways)",
                 fontsize=18, fontweight="bold")
    ax.tick_params(axis="x", labelsize=13)

    # Annotations
    ax.text(0.98, 0.02, "↑ Upregulated\nin disease", transform=ax.transAxes,
            fontsize=14, ha="right", va="bottom", color="#C0392B", fontweight="bold")
    ax.text(0.02, 0.02, "↓ Downregulated\nin disease", transform=ax.transAxes,
            fontsize=14, ha="left", va="bottom", color="#2E86AB", fontweight="bold")
    ax.grid(True, axis="x", alpha=0.15, ls="--")

    save_panel(fig, "panel_10_communication_waterfall", formats=("pdf",))


# ═══════════════════════════════════════════════════════════════════════════════
# COMPOSITE FIGURE
# ═══════════════════════════════════════════════════════════════════════════════
def make_composite():
    """Assemble key panels into a composite overview figure."""
    print("\n[Composite] Assembling overview figure from saved panels...")

    # Read saved panel PNGs
    panel_files = sorted(OUT_DIR.glob("panel_*.png"))
    if len(panel_files) < 3:
        print("  SKIP: Not enough panels saved")
        return

    from PIL import Image

    imgs = []
    for pf in panel_files:
        try:
            imgs.append((pf.stem, Image.open(pf)))
        except Exception:
            pass

    if not imgs:
        return

    # Simple 2-column grid layout
    n_cols = 2
    n_rows = (len(imgs) + 1) // 2
    max_w = max(img.width for _, img in imgs)
    max_h = max(img.height for _, img in imgs)

    # Scale to uniform cell size
    cell_w = 2400
    cell_h = 1600
    margin = 60

    canvas_w = n_cols * cell_w + (n_cols + 1) * margin
    canvas_h = n_rows * cell_h + (n_rows + 1) * margin + 120  # title space
    canvas = Image.new("RGB", (canvas_w, canvas_h), "white")

    for idx, (name, img) in enumerate(imgs):
        row = idx // n_cols
        col = idx % n_cols
        # Scale image to fit cell
        scale = min(cell_w / img.width, cell_h / img.height)
        new_w = int(img.width * scale)
        new_h = int(img.height * scale)
        img_resized = img.resize((new_w, new_h), Image.LANCZOS)

        x = col * cell_w + (col + 1) * margin + (cell_w - new_w) // 2
        y = row * cell_h + (row + 1) * margin + 120 + (cell_h - new_h) // 2
        canvas.paste(img_resized, (x, y))

    canvas.save(OUT_DIR / "spatial_composite.png", dpi=(300, 300))
    print(f"  ✓ Saved: spatial_composite.png ({canvas_w}x{canvas_h})")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════
def main():
    print("=" * 70)
    print("SPATIAL TRANSCRIPTOMICS — PRESENTATION PANELS")
    print("=" * 70)
    print(f"Output directory: {OUT_DIR}")
    print()

    # Load data objects
    print("Loading data...")
    adata_zon = sc.read_h5ad(RESULTS_DIR / "preprocessed" / "merged_spatial.h5ad")

    # Try zonation-enriched object
    zon_path = RESULTS_DIR / "zonation" / "spatial_with_zonation.h5ad"
    adata_z = sc.read_h5ad(zon_path) if zon_path.exists() else None

    # Domains
    dom_path = RESULTS_DIR / "domains" / "spatial_with_domains.h5ad"
    adata_d = sc.read_h5ad(dom_path) if dom_path.exists() else None

    # Trajectory
    traj_path = RESULTS_DIR / "trajectory" / "spatial_with_trajectory.h5ad"
    adata_t = sc.read_h5ad(traj_path) if traj_path.exists() else None

    print(f"  Merged spatial: {adata_zon.shape}")
    if adata_z is not None:
        print(f"  Zonation: {adata_z.shape}")
    if adata_d is not None:
        print(f"  Domains: {adata_d.shape}")
    if adata_t is not None:
        print(f"  Trajectory: {adata_t.shape}")
    print()

    # Use zonation-enriched object where available (has c2l columns)
    adata_main = adata_z if adata_z is not None else adata_zon

    # Generate panels
    panel_01_deconvolution(adata_main)
    panel_02_zonation(adata_main)
    panel_03_trajectory(adata_t)
    panel_04_domains(adata_d)
    panel_05_svg_heatmap()
    panel_06_commot_signaling()
    panel_07_gene_expression(adata_zon)
    panel_08_hep_validation()
    panel_09_zonation_celltypes(adata_main)
    panel_10_communication_waterfall()

    # Composite
    try:
        make_composite()
    except ImportError:
        print("  PIL not available — skipping composite assembly")

    print()
    print("=" * 70)
    print(f"ALL PANELS SAVED TO: {OUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()

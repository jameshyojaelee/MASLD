#!/usr/bin/env python3
"""
07_spatial_figures.py — Publication-quality spatial transcriptomics figures.

Main figure: 6 panels (deconvolution, zonation, hep validation, domains,
domain composition, cell type proportions).
Supplementary figure: 8 panels (QC, nhood enrichment, SVG volcano, trajectory,
zonation heatmap, co-expression, differential L-R, enrichment).
Standalone: enrichment bar, disruption boxplot, SVG comparison scatter,
top hep-intrinsic genes, communication matrix.

SLURM: --partition=cpu --cpus=8 --mem=64G --time=2:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.colors as mcolors
from matplotlib.colors import Normalize
from matplotlib import cm
from scipy.stats import pearsonr, spearmanr
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    PROJECT_ROOT, RESULTS_DIR, C2L_PREFIX, load_config, load_spatial_adata,
    load_deconvolved_adata, print_header,
)

# ── Nature-compatible plotting defaults ──────────────────────────────────────
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# F248: route figures into the canonical figures/ tree (figS05_epigenomic_spatial),
# NOT Analysis/Spatial/results/figures (which is outside the manuscript figure
# tree and not under the FIG_* constants in scripts/figures/load_figure_data.R).
FIG_DIR = PROJECT_ROOT / "figures" / "supplementary" / "figS05_epigenomic_spatial"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Color palette consistent with project publication_theme.R
CONDITION_COLORS = {
    "Healthy": "#4DAF4A", "Steatotic": "#E41A1C",
    "MASLD": "#E41A1C", "MASLD_spectrum": "#E41A1C",
}
ZONATION_CMAP = "RdBu_r"
CELL_TYPE_COLORS = {
    "Hepatocytes": "#1f77b4", "Fibroblasts": "#ff7f0e",
    "Endothelial cells": "#2ca02c", "Macrophages": "#d62728",
    "Cholangiocytes": "#9467bd", "T cells": "#8c564b",
    "B cells": "#e377c2", "Neutrophils": "#7f7f7f",
    "Circulating NK/NKT": "#bcbd22", "Plasma cells": "#17becf",
    "Mono+mono derived cells": "#aec7e8", "Resident NK": "#ffbb78",
    "cDC1s": "#98df8a", "cDC2s": "#ff9896", "pDCs": "#c5b0d5",
    "Basophils": "#c49c94",
}
def strip_c2l(name):
    """Strip cell2location column prefix."""
    return name[len(C2L_PREFIX):] if isinstance(name, str) and name.startswith(C2L_PREFIX) else name


def save_fig(fig, name, formats=("pdf",)):
    """Save figure as PDF only.

    F248: PNG output is prohibited by the project figure invariants (PDF only,
    no PNG). Default dropped from ("pdf","png") to ("pdf",).
    """
    formats = tuple(f for f in formats if f != "png")  # never emit PNG
    for fmt in formats:
        fig.savefig(FIG_DIR / f"{name}.{fmt}", format=fmt)
    plt.close(fig)
    print(f"  Saved: {name} ({', '.join(formats)})")


def _add_colorbar(fig, ax, mappable, label="", shrink=0.6, pad=0.02):
    """Add a compact colorbar."""
    cb = fig.colorbar(mappable, ax=ax, shrink=shrink, pad=pad, aspect=15)
    cb.ax.tick_params(labelsize=4, width=0.3)
    if label:
        cb.set_label(label, fontsize=5)
    return cb


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN FIGURE PANELS
# ═══════════════════════════════════════════════════════════════════════════════

def panel_deconvolution(adata, fig, gs):
    """Panel A: cell2location abundance maps — 6 cell types × 2 conditions.

    Uses one representative sample per condition. Includes colorbars.
    """
    cell_types = ["Hepatocytes", "Fibroblasts", "Endothelial cells",
                  "Macrophages", "Cholangiocytes", "T cells"]
    abundances = adata.obsm.get("q05_cell_abundance_w_sf")
    if abundances is None:
        return

    # Map cell type names to columns
    col_map = {}
    for ct in cell_types:
        matches = [c for c in abundances.columns if ct.lower() in c.lower()]
        if matches:
            col_map[ct] = matches[0]
    available = list(col_map.keys())
    if not available:
        return

    conditions = sorted(adata.obs["condition"].unique())

    # One representative sample per condition (most spots)
    rep = {}
    for cond in conditions:
        rep[cond] = adata.obs[adata.obs["condition"] == cond]["sample_id"].value_counts().index[0]

    inner = gs.subgridspec(len(conditions), len(available), wspace=0.08, hspace=0.2)

    for i, cond in enumerate(conditions):
        mask = (adata.obs["condition"] == cond) & (adata.obs["sample_id"] == rep[cond])
        coords = adata[mask].obsm["spatial"]
        for j, ct in enumerate(available):
            ax = fig.add_subplot(inner[i, j])
            vals = abundances.loc[mask, col_map[ct]].values
            vmax = np.percentile(vals, 97) if vals.max() > 0 else 1
            sc_map = ax.scatter(
                coords[:, 0], coords[:, 1], c=vals,
                cmap="magma", s=1.2, vmin=0, vmax=vmax,
                edgecolors="none", rasterized=True,
            )
            ax.set_aspect("equal")
            ax.axis("off")
            if i == 0:
                ax.set_title(ct, fontsize=5.5, fontweight="bold")
            if j == 0:
                ax.set_ylabel(cond, fontsize=6, rotation=90, labelpad=8)
                ax.yaxis.set_visible(True)
            # Colorbar on last column
            if j == len(available) - 1:
                _add_colorbar(fig, ax, sc_map, shrink=0.5, pad=0.05)


def panel_zonation(adata, fig, gs):
    """Panel B: zonation gradient on spatial coords with colorbars."""
    if "zonation_score" not in adata.obs.columns:
        return

    conditions = sorted(adata.obs["condition"].unique())
    inner = gs.subgridspec(1, len(conditions), wspace=0.15)

    for i, cond in enumerate(conditions):
        ax = fig.add_subplot(inner[0, i])
        mask = adata.obs["condition"] == cond
        coords = adata[mask].obsm["spatial"]
        vals = adata.obs.loc[mask, "zonation_score"].values

        vabs = max(abs(np.nanpercentile(vals, 5)), abs(np.nanpercentile(vals, 95)))
        sc_map = ax.scatter(
            coords[:, 0], coords[:, 1], c=vals,
            cmap=ZONATION_CMAP, s=1.0, vmin=-vabs, vmax=vabs,
            edgecolors="none", rasterized=True,
        )
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(f"{cond}", fontsize=6, fontweight="bold")
        cb = _add_colorbar(fig, ax, sc_map, label="PP ← → PC", shrink=0.5)


def panel_hep_validation(fig, gs):
    """Panel C: Bulk attribution_raw vs spatial_fc with regression line."""
    hep_path = RESULTS_DIR / "cell2location" / "hep_intrinsic_validation.csv"
    if not hep_path.exists():
        return

    ax = fig.add_subplot(gs)
    df = pd.read_csv(hep_path)
    x = df.get("attribution_raw")
    y = df.get("spatial_fc")
    if x is None or y is None:
        return

    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid].values, y[valid].values
    if len(x) < 5:
        return

    rho, p = spearmanr(x, y)
    validated = df.loc[valid, "validated"].values if "validated" in df.columns else np.zeros(len(x), dtype=bool)

    # Scatter with color by validation status
    ax.scatter(x[~validated], y[~validated], s=6, alpha=0.4, c="#AAAAAA",
               edgecolors="none", rasterized=True, label=f"Not validated (n={int((~validated).sum())})", zorder=2)
    ax.scatter(x[validated], y[validated], s=8, alpha=0.6, c="#2166AC",
               edgecolors="none", rasterized=True, label=f"Validated (n={int(validated.sum())})", zorder=3)

    # Regression line
    z = np.polyfit(x, y, 1)
    xline = np.linspace(x.min(), x.max(), 100)
    ax.plot(xline, np.polyval(z, xline), "--", c="#B2182B", lw=0.8, alpha=0.8, zorder=4)

    ax.axhline(1.0, ls=":", lw=0.4, c="gray", alpha=0.6)
    ax.set_xlabel("Bulk attribution score", fontsize=6)
    ax.set_ylabel("Spatial hepatocyte FC", fontsize=6)
    ax.set_title("Hepatocyte-intrinsic validation", fontsize=6, fontweight="bold")
    ax.annotate(f"Spearman ρ = {rho:.2f}\np = {p:.1e}\nn = {len(x)} genes",
                xy=(0.03, 0.97), xycoords="axes fraction", fontsize=5, va="top",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8, lw=0.3))
    ax.legend(fontsize=4.5, loc="lower right", frameon=True, framealpha=0.8,
              edgecolor="gray", markerscale=1.5, handletextpad=0.3)
    ax.tick_params(axis="both", labelsize=5)


def panel_domains(adata, fig, gs):
    """Panel D: Spatial domain maps with domain legend."""
    if "spatial_domain" not in adata.obs.columns:
        return

    conditions = sorted(adata.obs["condition"].unique())
    n_domains = adata.obs["spatial_domain"].nunique()
    cmap = plt.cm.get_cmap("Set2", n_domains)
    inner = gs.subgridspec(1, len(conditions) + 1, wspace=0.1,
                           width_ratios=[1] * len(conditions) + [0.3])

    for i, cond in enumerate(conditions):
        ax = fig.add_subplot(inner[0, i])
        mask = adata.obs["condition"] == cond
        coords = adata[mask].obsm["spatial"]
        domains = adata.obs.loc[mask, "spatial_domain"].astype(int).values
        ax.scatter(
            coords[:, 0], coords[:, 1], c=domains,
            cmap=cmap, s=1.0, edgecolors="none", rasterized=True,
            vmin=-0.5, vmax=n_domains - 0.5,
        )
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(f"{cond}", fontsize=6, fontweight="bold")

    # Legend panel
    ax_leg = fig.add_subplot(inner[0, -1])
    ax_leg.axis("off")
    for d in range(n_domains):
        ax_leg.scatter([], [], c=[cmap(d)], s=15, label=f"D{d}")
    ax_leg.legend(fontsize=4.5, loc="center", frameon=False, title="Domain",
                  title_fontsize=5, handletextpad=0.2, labelspacing=0.4)


def panel_domain_composition(fig, gs):
    """Panel E: Domain cell type composition — stacked bar with non-hep zoom."""
    comp_path = RESULTS_DIR / "domains" / "domain_composition.csv"
    if not comp_path.exists():
        return

    df = pd.read_csv(comp_path)
    if "spatial_domain" not in df.columns:
        return

    rename_map = {c: strip_c2l(c) for c in df.columns if c != "spatial_domain"}
    df = df.rename(columns=rename_map)
    ct_cols = [c for c in df.columns if c != "spatial_domain"]
    df = df.set_index("spatial_domain")

    inner = gs.subgridspec(1, 2, wspace=0.35, width_ratios=[1, 1])

    # Full composition
    ax1 = fig.add_subplot(inner[0, 0])
    bottom = np.zeros(len(df))
    for ct in ct_cols:
        color = CELL_TYPE_COLORS.get(ct, "#CCCCCC")
        ax1.bar(range(len(df)), df[ct].values, bottom=bottom, color=color,
                edgecolor="none", width=0.7, label=ct)
        bottom += df[ct].values
    ax1.set_xlabel("Domain", fontsize=5)
    ax1.set_ylabel("Proportion", fontsize=5)
    ax1.set_title("All cell types", fontsize=5.5, fontweight="bold")
    ax1.set_xticks(range(len(df)))
    ax1.set_xticklabels([f"D{i}" for i in df.index], fontsize=4)
    ax1.tick_params(axis="y", labelsize=4)
    ax1.set_ylim(0, 1.05)

    # Non-hepatocyte zoom
    ax2 = fig.add_subplot(inner[0, 1])
    non_hep_cols = [c for c in ct_cols if c != "Hepatocytes"]
    bottom2 = np.zeros(len(df))
    for ct in non_hep_cols:
        color = CELL_TYPE_COLORS.get(ct, "#CCCCCC")
        ax2.bar(range(len(df)), df[ct].values, bottom=bottom2, color=color,
                edgecolor="none", width=0.7, label=ct)
        bottom2 += df[ct].values
    ax2.set_xlabel("Domain", fontsize=5)
    ax2.set_ylabel("Proportion (excl. Hep)", fontsize=5)
    ax2.set_title("Non-hepatocyte types", fontsize=5.5, fontweight="bold")
    ax2.set_xticks(range(len(df)))
    ax2.set_xticklabels([f"D{i}" for i in df.index], fontsize=4)
    ax2.tick_params(axis="y", labelsize=4)
    ax2.legend(fontsize=3.5, loc="upper right", frameon=True, framealpha=0.9,
               ncol=1, handletextpad=0.2, labelspacing=0.2,
               edgecolor="gray", borderpad=0.3)


def panel_cell_type_proportions(fig, gs):
    """Panel F: Per-sample cell type proportions grouped by condition."""
    props_path = RESULTS_DIR / "cell2location" / "spatial_cell_type_proportions.csv"
    if not props_path.exists():
        return

    ax = fig.add_subplot(gs)
    df = pd.read_csv(props_path, index_col=0)
    df.columns = [strip_c2l(c) for c in df.columns]

    # Get sample metadata for condition coloring
    qc_path = RESULTS_DIR / "qc" / "spot_qc_summary.csv"
    if qc_path.exists():
        qc = pd.read_csv(qc_path)
        sample_cond = dict(zip(qc["sample_id"], qc["condition"]))
    else:
        sample_cond = {}

    # Sort by condition then sample
    df["condition"] = df.index.map(lambda x: sample_cond.get(x, "Unknown"))
    df = df.sort_values("condition")
    cond_labels = df["condition"].values
    df = df.drop(columns=["condition"])

    # Sort cell types by mean abundance
    ct_order = df.mean().sort_values(ascending=True).index.tolist()

    bottom = np.zeros(len(df))
    for ct in ct_order:
        color = CELL_TYPE_COLORS.get(ct, "#CCCCCC")
        ax.barh(range(len(df)), df[ct].values, left=bottom, color=color,
                edgecolor="none", height=0.7, label=ct)
        bottom += df[ct].values

    ax.set_yticks(range(len(df)))
    ylabels = [f"{s} ({c[0]})" for s, c in zip(df.index, cond_labels)]
    ax.set_yticklabels(ylabels, fontsize=4.5)
    ax.set_xlabel("Proportion", fontsize=5)
    ax.set_title("Cell type proportions per sample", fontsize=5.5, fontweight="bold")
    ax.tick_params(axis="x", labelsize=4)
    ax.set_xlim(0, 1.05)

    # Compact legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], fontsize=3, loc="lower right",
              frameon=True, framealpha=0.9, ncol=2, handletextpad=0.2,
              labelspacing=0.15, edgecolor="gray", borderpad=0.3)


# ═══════════════════════════════════════════════════════════════════════════════
# SUPPLEMENTARY PANELS
# ═══════════════════════════════════════════════════════════════════════════════

def supp_spot_qc(fig, gs):
    """S1: Per-sample QC metrics (genes, counts, pct_mt)."""
    qc_path = RESULTS_DIR / "qc" / "spot_qc_summary.csv"
    if not qc_path.exists():
        return

    df = pd.read_csv(qc_path)
    metrics = [("median_genes", "Median genes/spot"), ("median_counts", "Median UMIs/spot"),
               ("median_pct_mt", "% Mitochondrial")]
    available = [(col, label) for col, label in metrics if col in df.columns]
    if not available:
        return

    inner = gs.subgridspec(1, len(available), wspace=0.4)
    for i, (col, label) in enumerate(available):
        ax = fig.add_subplot(inner[0, i])
        colors = [CONDITION_COLORS.get(c, "#999999") for c in df["condition"]]
        ax.bar(range(len(df)), df[col], color=colors, edgecolor="none", width=0.6)
        ax.set_xticks(range(len(df)))
        ax.set_xticklabels(df["sample_id"], fontsize=4, rotation=45, ha="right")
        ax.set_ylabel(label, fontsize=5)
        ax.set_title(label, fontsize=5.5)
        ax.tick_params(axis="y", labelsize=4)

        # Add condition legend on first panel
        if i == 0:
            for cond, color in CONDITION_COLORS.items():
                if cond in df["condition"].values:
                    ax.scatter([], [], c=color, s=15, label=cond)
            ax.legend(fontsize=4, loc="upper right", frameon=False)


def supp_nhood_enrichment(fig, gs):
    """S2: Neighborhood enrichment heatmaps with clean labels."""
    import seaborn as sns

    conditions = []
    for pattern in sorted(RESULTS_DIR.glob("communication/nhood_enrichment_*_zscore.csv")):
        cond = pattern.stem.replace("nhood_enrichment_", "").replace("_zscore", "")
        conditions.append((cond, pd.read_csv(pattern, index_col=0)))

    if not conditions:
        return

    inner = gs.subgridspec(1, len(conditions), wspace=0.4)
    for i, (cond, df) in enumerate(conditions):
        ax = fig.add_subplot(inner[0, i])
        df.index = [strip_c2l(x) for x in df.index]
        df.columns = [strip_c2l(x) for x in df.columns]
        vmax = max(abs(df.values.min()), abs(df.values.max()))
        sns.heatmap(df, cmap="RdBu_r", center=0, ax=ax, vmin=-vmax, vmax=vmax,
                    xticklabels=True, yticklabels=True,
                    cbar_kws={"shrink": 0.5, "label": "z-score"},
                    linewidths=0.2, linecolor="white",
                    square=True)
        ax.set_title(f"{cond}", fontsize=5.5, fontweight="bold")
        ax.tick_params(axis="both", labelsize=3.5)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")


def supp_svg_volcano(fig, gs):
    """S3: SVG volcano — Moran's I vs -log10(FDR) with gene labels."""
    svg_files = sorted(RESULTS_DIR.glob("svg/svgs_*.csv"))
    if not svg_files:
        return

    inner = gs.subgridspec(1, len(svg_files), wspace=0.3)
    for i, svg_path in enumerate(svg_files):
        ax = fig.add_subplot(inner[0, i])
        df = pd.read_csv(svg_path, index_col=0)
        cond = svg_path.stem.replace("svgs_", "")

        if "I" not in df.columns or "pval_norm" not in df.columns:
            continue

        x = df["I"]
        y = -np.log10(df["pval_norm"].clip(lower=1e-300))
        sig = df.get("svg", pd.Series(False, index=df.index))

        ax.scatter(x[~sig], y[~sig], s=1, c="#CCCCCC", alpha=0.4,
                   edgecolors="none", rasterized=True)
        ax.scatter(x[sig], y[sig], s=2, c="#D62728", alpha=0.6,
                   edgecolors="none", rasterized=True)

        # Label top SVGs
        top = df[sig].nlargest(5, "I")
        for gene, row in top.iterrows():
            ax.annotate(gene, (row["I"], -np.log10(max(row["pval_norm"], 1e-300))),
                        fontsize=3.5, alpha=0.8, ha="left",
                        xytext=(3, 2), textcoords="offset points")

        n_sig = sig.sum()
        ax.set_xlabel("Moran's I", fontsize=5)
        ax.set_ylabel("-log₁₀(FDR)", fontsize=5)
        ax.set_title(f"{cond} (n={n_sig} SVGs)", fontsize=5.5, fontweight="bold")
        ax.tick_params(axis="both", labelsize=4)


def supp_trajectory(fig, gs):
    """S4: DPT pseudotime distribution + gene dynamics."""
    traj_path = RESULTS_DIR / "trajectory" / "spatial_pseudotime.csv"
    dynamics_path = RESULTS_DIR / "trajectory" / "gene_dynamics_along_trajectory.csv"
    if not traj_path.exists():
        return

    inner = gs.subgridspec(1, 2, wspace=0.35)

    ax = fig.add_subplot(inner[0, 0])
    traj = pd.read_csv(traj_path, index_col=0)
    if "dpt_pseudotime" in traj.columns:
        ax.hist(traj["dpt_pseudotime"].dropna(), bins=50, color="#2166AC",
                alpha=0.7, edgecolor="none")
        ax.set_xlabel("Pseudotime", fontsize=5)
        ax.set_ylabel("# Spots", fontsize=5)
        ax.set_title("DPT distribution", fontsize=5.5, fontweight="bold")
        ax.tick_params(axis="both", labelsize=4)

    if dynamics_path.exists():
        ax2 = fig.add_subplot(inner[0, 1])
        dyn = pd.read_csv(dynamics_path)
        colors = plt.cm.tab10(np.linspace(0, 1, 6))
        for idx, gene in enumerate(dyn["gene"].unique()[:6]):
            gd = dyn[dyn["gene"] == gene]
            ax2.plot(gd["mean_dpt"], gd["mean_expr"], label=gene, lw=1.0, color=colors[idx])
        ax2.legend(fontsize=4, loc="upper right", frameon=True, framealpha=0.8,
                   edgecolor="gray")
        ax2.set_xlabel("Pseudotime", fontsize=5)
        ax2.set_ylabel("Mean expression", fontsize=5)
        ax2.set_title("Gene dynamics", fontsize=5.5, fontweight="bold")
        ax2.tick_params(axis="both", labelsize=4)


def supp_zonation_heatmap(fig, gs):
    """S5: Zonation gene expression heatmap — top PP/PC genes across bins."""
    zon_path = RESULTS_DIR / "zonation" / "deg_zonation_classification.csv"
    if not zon_path.exists():
        return

    ax = fig.add_subplot(gs)
    df = pd.read_csv(zon_path, index_col=0)

    # Get top periportal and pericentral genes by absolute spearman_rho
    pp = df[df["zonation_class"] == "Periportal-enriched"].nlargest(15, "kruskal_stat")
    pc = df[df["zonation_class"] == "Pericentral-enriched"].nlargest(15, "kruskal_stat")
    selected = pd.concat([pp, pc])

    if len(selected) < 3:
        ax.text(0.5, 0.5, "Insufficient zonation\ngenes for heatmap",
                ha="center", va="center", fontsize=6, transform=ax.transAxes)
        return

    bin_cols = ["mean_PP1", "mean_PP2", "mean_Mid", "mean_PC2", "mean_PC1"]
    available_bins = [c for c in bin_cols if c in selected.columns]
    if not available_bins:
        return

    heatmap_data = selected[available_bins].copy()
    heatmap_data.index = selected["gene"].values

    # Z-score normalize each gene (row)
    row_means = heatmap_data.mean(axis=1)
    row_stds = heatmap_data.std(axis=1).clip(lower=1e-6)
    heatmap_data = heatmap_data.sub(row_means, axis=0).div(row_stds, axis=0)

    import seaborn as sns
    vmax = max(abs(heatmap_data.values.min()), abs(heatmap_data.values.max()))
    sns.heatmap(heatmap_data, cmap="RdBu_r", center=0, ax=ax, vmin=-vmax, vmax=vmax,
                xticklabels=[c.replace("mean_", "") for c in available_bins],
                yticklabels=True,
                cbar_kws={"shrink": 0.5, "label": "Z-score"},
                linewidths=0.3, linecolor="white")
    ax.set_title(f"Zonation gene expression\n({len(pp)} PP + {len(pc)} PC genes)", fontsize=5.5, fontweight="bold")
    ax.tick_params(axis="both", labelsize=4)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

    # Color-code y-axis labels
    for i, label in enumerate(ax.get_yticklabels()):
        gene = label.get_text()
        if gene in pp["gene"].values:
            label.set_color("#2166AC")  # Blue = periportal
        else:
            label.set_color("#B2182B")  # Red = pericentral


def supp_coexpression_modules(fig, gs):
    """S6: Co-expression module sizes with condition comparison."""
    module_files = sorted(RESULTS_DIR.glob("coexpression/modules_*.csv"))
    if not module_files:
        return

    inner = gs.subgridspec(1, len(module_files[:2]), wspace=0.4)
    for i, mod_path in enumerate(module_files[:2]):
        ax = fig.add_subplot(inner[0, i])
        cond = mod_path.stem.replace("modules_", "")
        df = pd.read_csv(mod_path, index_col=0)
        if "Module" not in df.columns:
            continue

        mod_sizes = df["Module"].value_counts().sort_index()
        colors = plt.cm.Set3(np.linspace(0, 1, len(mod_sizes)))
        ax.bar(range(len(mod_sizes)), mod_sizes.values, color=colors,
               edgecolor="none", width=0.7)
        ax.set_xlabel("Module", fontsize=5)
        ax.set_ylabel("# Genes", fontsize=5)
        ax.set_title(f"{cond} ({len(mod_sizes)} modules, {len(df)} genes)",
                     fontsize=5, fontweight="bold")
        ax.tick_params(axis="both", labelsize=4)
        ax.set_xticks(range(len(mod_sizes)))
        ax.set_xticklabels([str(m) for m in mod_sizes.index], fontsize=3.5)


def supp_differential_lr(fig, gs):
    """S7: Differential L-R — top emergent and lost pairs."""
    lr_path = RESULTS_DIR / "communication" / "differential_lr_pairs.csv"
    if not lr_path.exists():
        return

    ax = fig.add_subplot(gs)
    df = pd.read_csv(lr_path, index_col=0)
    if len(df) == 0 or "delta_mean_expr" not in df.columns:
        return

    df = df.dropna(subset=["delta_mean_expr"])
    # Top 10 emergent and top 10 lost by absolute delta
    emergent = df[df["category"] == "disease_emergent"].nlargest(10, "delta_mean_expr")
    lost = df[df["category"] == "disease_lost"].nsmallest(10, "delta_mean_expr")
    top = pd.concat([lost, emergent])

    if len(top) == 0:
        return

    # Clean labels
    labels = []
    for _, row in top.iterrows():
        lr = str(row["lr_pair"]).replace("(", "").replace(")", "").replace("'", "")
        src = str(row.get("source", "?"))[:8]
        tgt = str(row.get("target", "?"))[:8]
        labels.append(f"{lr}\n({src}→{tgt})")

    colors = ["#B2182B" if c == "disease_emergent" else "#2166AC" for c in top["category"]]

    y_pos = range(len(top))
    ax.barh(y_pos, top["delta_mean_expr"].values, color=colors, edgecolor="none", height=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=3.5)
    ax.axvline(0, ls="-", lw=0.5, c="black")
    ax.set_xlabel("Δ mean expression (disease − healthy)", fontsize=5)
    ax.set_title(f"Top differential L-R pairs\n({len(emergent)} emergent, {len(lost)} lost)",
                 fontsize=5.5, fontweight="bold")
    ax.tick_params(axis="x", labelsize=4)

    # Legend
    from matplotlib.patches import Patch
    ax.legend([Patch(fc="#B2182B"), Patch(fc="#2166AC")],
              ["Disease-emergent", "Disease-lost"],
              fontsize=4, loc="lower right", frameon=True, framealpha=0.9)


def supp_enrichment(fig, gs):
    """S8: Enrichment bar chart — SVG enrichment across gene sets, BH-corrected."""
    enrich_path = RESULTS_DIR / "integration" / "spatial_enrichment_tests.csv"
    if not enrich_path.exists():
        return

    ax = fig.add_subplot(gs)
    df = pd.read_csv(enrich_path)
    if len(df) == 0:
        return

    # Only show SVG-level enrichments (not zonation cross-tests which clutter)
    svg_tests = df[~df["gene_set"].str.contains("enriched_in_")]
    zon_tests = df[df["gene_set"].str.contains("enriched_in_")]

    # Combine — put SVG tests first, then zonation
    plot_df = pd.concat([svg_tests.sort_values("odds_ratio"),
                         zon_tests.sort_values("odds_ratio")])

    log2_or = np.log2(plot_df["odds_ratio"].clip(lower=0.01))
    p_col = "fisher_padj_bh" if "fisher_padj_bh" in plot_df.columns else "fisher_pval"

    # Color by significance
    colors = []
    for _, row in plot_df.iterrows():
        p = row[p_col]
        if p < 0.001:
            colors.append("#B2182B")
        elif p < 0.01:
            colors.append("#EF8A62")
        elif p < 0.05:
            colors.append("#67A9CF")
        else:
            colors.append("#D1D1D1")

    # Clean labels
    clean_labels = []
    for name in plot_df["gene_set"]:
        name = name.replace("_", " ").replace("enriched in ", "∩ ")
        if "Periportal" in name:
            name = name.replace("Periportal ", "PP ∩ ")
        if "Pericentral" in name:
            name = name.replace("Pericentral ", "PC ∩ ")
        clean_labels.append(name)

    ax.barh(range(len(plot_df)), log2_or.values, color=colors, edgecolor="none", height=0.7)
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(clean_labels, fontsize=4)
    ax.axvline(0, ls="-", lw=0.5, c="black")
    ax.set_xlabel("log₂(Odds Ratio)", fontsize=5)
    ax.set_title("Spatial enrichment tests", fontsize=5.5, fontweight="bold")
    ax.tick_params(axis="x", labelsize=4)

    # Annotate significance
    for i, (_, row) in enumerate(plot_df.iterrows()):
        p = row[p_col]
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        if sig:
            x_pos = log2_or.values[i]
            offset = 0.08
            ax.text(x_pos + offset, i, sig, fontsize=4, va="center", color="#333333")

    # Significance legend
    from matplotlib.patches import Patch
    ax.legend(
        [Patch(fc="#B2182B"), Patch(fc="#EF8A62"), Patch(fc="#67A9CF"), Patch(fc="#D1D1D1")],
        ["p < 0.001", "p < 0.01", "p < 0.05", "ns"],
        fontsize=3.5, loc="lower right", frameon=True, framealpha=0.9,
        title="BH-adjusted", title_fontsize=4, edgecolor="gray")


# ═══════════════════════════════════════════════════════════════════════════════
# STANDALONE FIGURES
# ═══════════════════════════════════════════════════════════════════════════════

def plot_enrichment_bar():
    """Standalone enrichment bar chart (SVG-level only, cleaner)."""
    enrich_path = RESULTS_DIR / "integration" / "spatial_enrichment_tests.csv"
    if not enrich_path.exists():
        return

    df = pd.read_csv(enrich_path)
    # SVG-level only
    df = df[~df["gene_set"].str.contains("enriched_in_")]
    if len(df) == 0:
        return

    df = df.sort_values("odds_ratio")
    p_col = "fisher_padj_bh" if "fisher_padj_bh" in df.columns else "fisher_pval"
    log2_or = np.log2(df["odds_ratio"].clip(lower=0.01))
    colors = ["#2166AC" if p < 0.05 else "#D1D1D1" for p in df[p_col]]

    fig, ax = plt.subplots(figsize=(90 / 25.4, 55 / 25.4))
    ax.barh(range(len(df)), log2_or.values, color=colors, edgecolor="none", height=0.6)
    labels = [n.replace("_", " ") for n in df["gene_set"]]
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(labels, fontsize=5.5)
    ax.axvline(0, ls="-", lw=0.5, c="black")
    ax.set_xlabel("log₂(Odds Ratio)", fontsize=6)
    ax.set_title("SVG enrichment in bulk gene sets", fontsize=7, fontweight="bold")
    ax.tick_params(axis="x", labelsize=5)

    for i, (_, row) in enumerate(df.iterrows()):
        sig = "***" if row[p_col] < 0.001 else "**" if row[p_col] < 0.01 else "*" if row[p_col] < 0.05 else "ns"
        ax.text(log2_or.values[i] + 0.08, i, sig, fontsize=5, va="center")

    save_fig(fig, "enrichment_bar")


def plot_zonation_disruption():
    """Standalone: per-gene zonation disruption for known markers."""
    dis_path = RESULTS_DIR / "zonation" / "zonation_disruption_scores.csv"
    if not dis_path.exists():
        return

    df = pd.read_csv(dis_path)
    if "disruption_score" not in df.columns or "gene" not in df.columns:
        return

    # Per-gene mean disruption
    gene_disrupt = df.groupby("gene")["disruption_score"].agg(["mean", "std", "count"])
    gene_disrupt = gene_disrupt.sort_values("mean", ascending=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(120 / 25.4, 60 / 25.4),
                                    gridspec_kw={"width_ratios": [1.5, 1]})

    # Left: Per-gene bar chart
    colors = ["#B2182B" if m > 1 else "#2166AC" for m in gene_disrupt["mean"]]
    ax1.barh(range(len(gene_disrupt)), gene_disrupt["mean"], xerr=gene_disrupt["std"],
             color=colors, edgecolor="none", height=0.6, error_kw={"lw": 0.5, "capsize": 1.5})
    ax1.set_yticks(range(len(gene_disrupt)))
    ax1.set_yticklabels(gene_disrupt.index, fontsize=5)
    ax1.axvline(1, ls="--", lw=0.5, c="gray")
    ax1.set_xlabel("Disruption score", fontsize=6)
    ax1.set_title("Per-gene zonation disruption\n(Steatotic vs Healthy)", fontsize=6, fontweight="bold")
    ax1.tick_params(axis="x", labelsize=5)

    # Right: Overall boxplot
    ax2.boxplot(df["disruption_score"].values, widths=0.5, patch_artist=True,
                boxprops=dict(facecolor="#E41A1C", alpha=0.5),
                medianprops=dict(color="black", lw=1))
    ax2.scatter(np.ones(len(df)) + np.random.normal(0, 0.04, len(df)),
                df["disruption_score"].values, s=3, c="#333333", alpha=0.4, zorder=3)
    ax2.axhline(1, ls="--", lw=0.5, c="gray")
    ax2.set_ylabel("Disruption score", fontsize=6)
    ax2.set_title("All markers", fontsize=6, fontweight="bold")
    ax2.set_xticklabels(["Steatotic"], fontsize=5)
    ax2.tick_params(axis="y", labelsize=5)

    plt.tight_layout()
    save_fig(fig, "zonation_disruption_boxplot")


def plot_svg_condition_comparison():
    """Standalone: Scatter plot comparing Moran's I between conditions."""
    svg_h_path = RESULTS_DIR / "svg" / "svgs_Healthy.csv"
    svg_m_path = RESULTS_DIR / "svg" / "svgs_Steatotic.csv"
    if not svg_h_path.exists() or not svg_m_path.exists():
        return

    df_h = pd.read_csv(svg_h_path, index_col=0)
    df_m = pd.read_csv(svg_m_path, index_col=0)
    shared = df_h.index.intersection(df_m.index)
    if len(shared) < 10:
        return

    x = df_h.loc[shared, "I"]
    y = df_m.loc[shared, "I"]
    svg_h_only = df_h.loc[shared, "svg"] & ~df_m.loc[shared, "svg"]
    svg_m_only = ~df_h.loc[shared, "svg"] & df_m.loc[shared, "svg"]
    svg_both = df_h.loc[shared, "svg"] & df_m.loc[shared, "svg"]
    svg_none = ~df_h.loc[shared, "svg"] & ~df_m.loc[shared, "svg"]

    rho, p = spearmanr(x, y)

    fig, ax = plt.subplots(figsize=(80 / 25.4, 80 / 25.4))
    ax.scatter(x[svg_none], y[svg_none], s=1, c="#D1D1D1", alpha=0.3,
               edgecolors="none", rasterized=True, label="Neither")
    ax.scatter(x[svg_both], y[svg_both], s=4, c="#2166AC", alpha=0.5,
               edgecolors="none", rasterized=True, label=f"Both (n={svg_both.sum()})")
    ax.scatter(x[svg_h_only], y[svg_h_only], s=4, c="#4DAF4A", alpha=0.6,
               edgecolors="none", rasterized=True, label=f"Healthy only (n={svg_h_only.sum()})")
    ax.scatter(x[svg_m_only], y[svg_m_only], s=4, c="#E41A1C", alpha=0.6,
               edgecolors="none", rasterized=True, label=f"Steatotic only (n={svg_m_only.sum()})")

    # Diagonal
    lim = max(x.max(), y.max()) * 1.05
    ax.plot([0, lim], [0, lim], "--", c="gray", lw=0.5, alpha=0.5)

    # Label top differential SVGs
    diff_svg_path = RESULTS_DIR / "svg" / "differential_svgs.csv"
    if diff_svg_path.exists():
        diff = pd.read_csv(diff_svg_path, index_col=0)
        top_diff = diff.nlargest(5, "delta_I").index.tolist() + diff.nsmallest(5, "delta_I").index.tolist()
        for gene in top_diff:
            if gene in shared:
                ax.annotate(gene, (x[gene], y[gene]), fontsize=3.5, alpha=0.8,
                            xytext=(3, 3), textcoords="offset points")

    ax.set_xlabel("Moran's I (Healthy)", fontsize=6)
    ax.set_ylabel("Moran's I (Steatotic)", fontsize=6)
    ax.set_title("Spatial variability comparison", fontsize=7, fontweight="bold")
    ax.annotate(f"ρ = {rho:.2f}, p = {p:.1e}\nn = {len(shared)}",
                xy=(0.03, 0.97), xycoords="axes fraction", fontsize=5, va="top",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8, lw=0.3))
    ax.legend(fontsize=4.5, loc="lower right", frameon=True, framealpha=0.9,
              edgecolor="gray", markerscale=2)
    ax.tick_params(axis="both", labelsize=5)
    ax.set_xlim(-0.02, lim)
    ax.set_ylim(-0.02, lim)

    save_fig(fig, "svg_condition_comparison")


def plot_top_hep_intrinsic():
    """Standalone: Top validated hepatocyte-intrinsic genes by spatial FC."""
    hep_path = RESULTS_DIR / "cell2location" / "hep_intrinsic_validation.csv"
    if not hep_path.exists():
        return

    df = pd.read_csv(hep_path)
    if "spatial_fc" not in df.columns or "gene" not in df.columns:
        return

    # Top 30 by spatial_fc
    validated = df[df.get("validated", False) == True]
    if len(validated) < 5:
        return
    top = validated.nlargest(30, "spatial_fc")

    fig, ax = plt.subplots(figsize=(90 / 25.4, 90 / 25.4))
    top_sorted = top.sort_values("spatial_fc")
    colors = ["#2166AC" if p < 0.001 else "#67A9CF" if p < 0.01 else "#D1E5F0"
              for p in top_sorted["wilcoxon_padj_bh"]]
    ax.barh(range(len(top_sorted)), top_sorted["spatial_fc"].values,
            color=colors, edgecolor="none", height=0.7)
    ax.set_yticks(range(len(top_sorted)))
    ax.set_yticklabels(top_sorted["gene"].values, fontsize=4.5)
    ax.axvline(1.0, ls="--", lw=0.5, c="gray")
    ax.set_xlabel("Spatial hepatocyte FC", fontsize=6)
    ax.set_title(f"Top validated hep-intrinsic genes\n(n={len(validated)}/{len(df)} validated)",
                 fontsize=6.5, fontweight="bold")
    ax.tick_params(axis="x", labelsize=5)

    from matplotlib.patches import Patch
    ax.legend([Patch(fc="#2166AC"), Patch(fc="#67A9CF"), Patch(fc="#D1E5F0")],
              ["p < 0.001", "p < 0.01", "p < 0.05"],
              fontsize=4, loc="lower right", frameon=True, framealpha=0.9,
              title="BH p-value", title_fontsize=4.5)

    save_fig(fig, "top_hep_intrinsic_genes")


def plot_communication_matrix():
    """Standalone: Cell-cell communication heatmap comparing conditions."""
    net_files = sorted(RESULTS_DIR.glob("communication/network_*.csv"))
    if len(net_files) < 2:
        return

    fig, axes = plt.subplots(1, len(net_files), figsize=(140 / 25.4, 65 / 25.4))
    if len(net_files) == 1:
        axes = [axes]

    for i, net_path in enumerate(net_files):
        ax = axes[i]
        cond = net_path.stem.replace("network_", "")
        df = pd.read_csv(net_path, index_col=0)
        if len(df) == 0:
            continue

        # Build adjacency matrix
        cell_types = sorted(set(df["source"].tolist() + df["target"].tolist()))
        n = len(cell_types)
        mat = pd.DataFrame(0.0, index=cell_types, columns=cell_types)
        for _, row in df.iterrows():
            mat.loc[row["source"], row["target"]] = row["weight"]

        import seaborn as sns
        sns.heatmap(mat, cmap="YlOrRd", ax=ax, annot=True, fmt=".0f",
                    annot_kws={"size": 4}, linewidths=0.3, linecolor="white",
                    cbar_kws={"shrink": 0.5, "label": "# L-R pairs"},
                    square=True)
        ax.set_title(f"{cond}", fontsize=6, fontweight="bold")
        ax.tick_params(axis="both", labelsize=4)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

    plt.tight_layout()
    save_fig(fig, "communication_matrix")


# ═══════════════════════════════════════════════════════════════════════════════
# ASSEMBLY
# ═══════════════════════════════════════════════════════════════════════════════

def assemble_main_figure(adata):
    """Main spatial figure: 6 panels (A-F)."""
    print("\n  Assembling main figure...")
    fig = plt.figure(figsize=(180 / 25.4, 300 / 25.4))
    gs = GridSpec(6, 2, figure=fig, hspace=0.4, wspace=0.3,
                  height_ratios=[1.5, 0.8, 1, 1, 1, 1])

    # A: Deconvolution (full width)
    panel_deconvolution(adata, fig, gs[0, :])

    # B: Zonation (full width)
    panel_zonation(adata, fig, gs[1, :])

    # C: Hep validation
    panel_hep_validation(fig, gs[2, 0])

    # D: Spatial domains
    panel_domains(adata, fig, gs[2, 1])

    # E: Domain composition (full width)
    panel_domain_composition(fig, gs[3, :])

    # F: Cell type proportions (full width)
    panel_cell_type_proportions(fig, gs[4, :])

    # Panel labels
    label_positions = [
        ("A", 0.02, 0.97), ("B", 0.02, 0.73), ("C", 0.02, 0.58),
        ("D", 0.52, 0.58), ("E", 0.02, 0.43), ("F", 0.02, 0.28),
    ]
    for label, x, y in label_positions:
        fig.text(x, y, label, fontsize=10, fontweight="bold", va="top")

    save_fig(fig, "Fig_spatial_main")


def assemble_supplementary(adata):
    """Supplementary spatial figure: 8 panels (S1-S8)."""
    print("\n  Assembling supplementary figure...")
    fig = plt.figure(figsize=(180 / 25.4, 480 / 25.4))
    gs = GridSpec(8, 1, figure=fig, hspace=0.5,
                  height_ratios=[0.7, 1, 0.8, 0.8, 1, 0.8, 1, 1.2])

    supp_spot_qc(fig, gs[0])            # S1
    supp_nhood_enrichment(fig, gs[1])    # S2
    supp_svg_volcano(fig, gs[2])         # S3
    supp_trajectory(fig, gs[3])          # S4
    supp_zonation_heatmap(fig, gs[4])    # S5
    supp_coexpression_modules(fig, gs[5])  # S6
    supp_differential_lr(fig, gs[6])     # S7
    supp_enrichment(fig, gs[7])          # S8

    # Panel labels
    for i, label in enumerate([f"S{j+1}" for j in range(8)]):
        fig.text(0.02, 0.99 - i * 0.125, label, fontsize=9, fontweight="bold", va="top")

    save_fig(fig, "Fig_spatial_supp")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    print_header("07: Spatial Publication Figures")

    config = load_config()

    # Load data
    try:
        adata = load_deconvolved_adata()
        print(f"  Loaded deconvolved: {adata.n_obs} spots")

        dom_path = RESULTS_DIR / "domains" / "spatial_domains.csv"
        if dom_path.exists():
            domains = pd.read_csv(dom_path, index_col=0)
            if "spatial_domain" in domains.columns:
                shared = adata.obs.index.intersection(domains.index)
                adata.obs.loc[shared, "spatial_domain"] = domains.loc[shared, "spatial_domain"]
    except FileNotFoundError:
        try:
            adata = load_spatial_adata()
            print(f"  Loaded preprocessed: {adata.n_obs} spots")
        except FileNotFoundError:
            print("  ERROR: No spatial data found.")
            return

    # Load zonation
    zon_path = RESULTS_DIR / "zonation" / "zonation_scores.csv"
    if zon_path.exists():
        zon = pd.read_csv(zon_path, index_col=0)
        if "zonation_score" in zon.columns:
            shared = adata.obs.index.intersection(zon.index)
            adata.obs.loc[shared, "zonation_score"] = zon.loc[shared, "zonation_score"]

    # Composite figures
    assemble_main_figure(adata)
    assemble_supplementary(adata)

    # Standalone figures
    plot_enrichment_bar()
    plot_zonation_disruption()
    plot_svg_condition_comparison()
    plot_top_hep_intrinsic()
    plot_communication_matrix()

    print_header("07: Complete")


if __name__ == "__main__":
    main()

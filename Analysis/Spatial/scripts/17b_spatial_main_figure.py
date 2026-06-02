#!/usr/bin/env python3
"""
17b_spatial_main_figure.py — Generate the main spatial figure (6 panels, 2x3).

Panels:
  (a) gsMap risk heatmap: tissue colored by -log10(p), Healthy vs Steatotic
  (b) Multi-GWAS zonation enrichment heatmap: 7 GWAS traits x zonation bins
  (c) Niche trajectory spatial maps: NT score, Healthy vs Disease
  (d) Disease trajectory shift: KDE + cell-type composition along trajectory
  (e) TGFB signaling flow: direction arrows on tissue
  (f) Multi-method convergence UpSet plot

Outputs:
  - scripts/figures/output/fig_spatial_main.pdf (combined figure)
  - scripts/figures/output/fig_spatial_panel_{a-f}.pdf (individual panels)

SLURM: --partition=cpu --cpus=8 --mem=64G --time=4:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize
from matplotlib import cm
import matplotlib.colors as mcolors
from matplotlib.patches import FancyArrowPatch

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    PROJECT_ROOT, RESULTS_DIR, load_config, load_spatial_adata,
    load_deconvolved_adata, print_header, print_step,
)

# ── Nature-compatible plotting defaults ──────────────────────────────────────
plt.rcParams.update({
    "font.size": 8,
    "axes.titlesize": 10,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica"],
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

FIG_OUTPUT_DIR = PROJECT_ROOT / "scripts/figures/output"
FIG_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CONDITION_COLORS = {
    "Healthy": "#4DAF4A",
    "Steatotic": "#E41A1C",
    "MASLD": "#E41A1C",
}


def save_fig(fig, name, directory=None, formats=("pdf",)):
    """Save figure in specified formats."""
    out = directory or FIG_OUTPUT_DIR
    out.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        fig.savefig(out / f"{name}.{fmt}", format=fmt)
    plt.close(fig)
    print(f"  Saved: {name} ({', '.join(formats)})")


def _add_colorbar(fig, ax, mappable, label="", shrink=0.6, pad=0.02):
    """Add a compact colorbar."""
    cb = fig.colorbar(mappable, ax=ax, shrink=shrink, pad=pad, aspect=15)
    cb.ax.tick_params(labelsize=5, width=0.3)
    if label:
        cb.set_label(label, fontsize=6)
    return cb


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def panel_a_gsmap_risk(fig, gs):
    """Panel (a): gsMap risk heatmap — tissue colored by -log10(p), side-by-side.

    Loads gsMap per-spot results and spatial coordinates, renders scatter plots
    for Healthy vs Steatotic with shared colorbar.
    """
    ax_container = fig.add_subplot(gs)
    ax_container.set_axis_off()

    # Try loading gsMap per-spot results
    gsmap_dir = RESULTS_DIR / "gsmap"
    spot_results = None
    for candidate in [
        gsmap_dir / "gsmap_spot_pvalues.csv",
        gsmap_dir / "per_spot_gsmap_results.csv",
        gsmap_dir / "gsmap_results_per_spot.csv",
    ]:
        if candidate.exists():
            spot_results = pd.read_csv(candidate, index_col=0)
            break

    if spot_results is None:
        ax_container.text(
            0.5, 0.5, "gsMap per-spot results\nnot found",
            ha="center", va="center", fontsize=8,
            transform=ax_container.transAxes, color="gray",
        )
        return

    # Identify NAFLD/MASLD GWAS p-value column
    pval_col = None
    for candidate in ["NAFLD_pvalue", "nafld_pvalue", "pvalue", "p_value"]:
        if candidate in spot_results.columns:
            pval_col = candidate
            break
    if pval_col is None:
        # Use first numeric column
        num_cols = spot_results.select_dtypes(include=[np.number]).columns
        if len(num_cols) > 0:
            pval_col = num_cols[0]
        else:
            ax_container.text(
                0.5, 0.5, "No p-value column\nfound in gsMap results",
                ha="center", va="center", fontsize=8,
                transform=ax_container.transAxes, color="gray",
            )
            return

    # Need spatial coordinates and condition labels
    coord_cols = [c for c in spot_results.columns if c in ["x", "y", "array_row", "array_col"]]
    cond_col = None
    for candidate in ["condition", "group", "sample_condition"]:
        if candidate in spot_results.columns:
            cond_col = candidate
            break

    if len(coord_cols) < 2 or cond_col is None:
        # Load from spatial adata
        try:
            adata = load_deconvolved_adata()
        except FileNotFoundError:
            try:
                adata = load_spatial_adata()
            except FileNotFoundError:
                ax_container.text(
                    0.5, 0.5, "Cannot load\nspatial coordinates",
                    ha="center", va="center", fontsize=8,
                    transform=ax_container.transAxes, color="gray",
                )
                return

        shared_idx = spot_results.index.intersection(adata.obs_names)
        if len(shared_idx) == 0:
            return

        spot_results = spot_results.loc[shared_idx]
        spot_results["x"] = adata[shared_idx].obsm["spatial"][:, 0]
        spot_results["y"] = adata[shared_idx].obsm["spatial"][:, 1]
        if cond_col is None:
            spot_results["condition"] = adata.obs.loc[shared_idx, "condition"].values
            cond_col = "condition"

    # Compute -log10(p)
    neg_log10_p = -np.log10(spot_results[pval_col].clip(lower=1e-300))

    conditions = sorted(spot_results[cond_col].unique())
    inner = gs.subgridspec(1, len(conditions), wspace=0.15)

    vmax = np.percentile(neg_log10_p.values, 99)
    norm = Normalize(vmin=0, vmax=vmax)

    for i, cond in enumerate(conditions):
        ax = fig.add_subplot(inner[0, i])
        mask = spot_results[cond_col] == cond
        x = spot_results.loc[mask, "x"].values
        y = spot_results.loc[mask, "y"].values
        vals = neg_log10_p.loc[mask].values

        sc_map = ax.scatter(
            x, y, c=vals, cmap="YlOrRd", s=1.2, norm=norm,
            edgecolors="none", rasterized=True,
        )
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(cond, fontsize=8, fontweight="bold")

        if i == len(conditions) - 1:
            _add_colorbar(fig, ax, sc_map, label="$-\\log_{10}(p)$", shrink=0.5, pad=0.05)


def panel_b_gwas_zonation(fig, gs):
    """Panel (b): Multi-GWAS zonation enrichment heatmap.

    Rows = GWAS traits, Columns = zonation bins (PP1 -> PC1).
    """
    import seaborn as sns

    ax = fig.add_subplot(gs)

    heatmap_path = RESULTS_DIR / "gsmap" / "gwas_zonation_heatmap.csv"
    if not heatmap_path.exists():
        # Fallback: try alternative locations
        alt_path = RESULTS_DIR / "gsmap" / "gwas_zonation_enrichment.csv"
        if alt_path.exists():
            heatmap_path = alt_path
        else:
            ax.text(
                0.5, 0.5, "gwas_zonation_heatmap.csv\nnot found",
                ha="center", va="center", fontsize=8,
                transform=ax.transAxes, color="gray",
            )
            return

    df = pd.read_csv(heatmap_path, index_col=0)

    # Ensure column order is PP -> PC
    bin_order = ["PP1", "PP2", "Mid", "PC2", "PC1"]
    available_bins = [b for b in bin_order if b in df.columns]
    if not available_bins:
        # Try without PP/PC prefix
        available_bins = [c for c in df.columns if c not in ["gwas", "trait"]]
    if not available_bins:
        ax.text(0.5, 0.5, "No zonation bins found", ha="center", va="center",
                fontsize=8, transform=ax.transAxes, color="gray")
        return

    plot_df = df[available_bins]

    vmax = max(abs(plot_df.values.min()), abs(plot_df.values.max()))
    sns.heatmap(
        plot_df, cmap="RdBu_r", center=0, ax=ax,
        vmin=-vmax, vmax=vmax,
        xticklabels=True, yticklabels=True,
        cbar_kws={"shrink": 0.6, "label": "Enrichment\n($-\\log_{10}(p)$ * sign)"},
        linewidths=0.3, linecolor="white",
        annot=True, fmt=".1f", annot_kws={"size": 5},
    )
    ax.set_title("GWAS risk enrichment\nacross zonation", fontsize=8, fontweight="bold")
    ax.set_xlabel("Zonation bin (PP -> PC)", fontsize=7)
    ax.set_ylabel("")
    ax.tick_params(axis="both", labelsize=6)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)


def panel_c_niche_trajectory(fig, gs):
    """Panel (c): Niche trajectory spatial maps — NT score, Healthy vs Disease."""
    ax_container = fig.add_subplot(gs)
    ax_container.set_axis_off()

    nt_path = RESULTS_DIR / "ontrac" / "niche_trajectory_scores.csv"
    if not nt_path.exists():
        ax_container.text(
            0.5, 0.5, "niche_trajectory_scores.csv\nnot found",
            ha="center", va="center", fontsize=8,
            transform=ax_container.transAxes, color="gray",
        )
        return

    nt = pd.read_csv(nt_path, index_col=0)

    # Need spatial coords and condition
    coord_cols = [c for c in nt.columns if c in ["x", "y"]]
    nt_score_col = None
    for candidate in ["NT_score", "nt_score", "niche_trajectory_score", "trajectory_score"]:
        if candidate in nt.columns:
            nt_score_col = candidate
            break
    if nt_score_col is None:
        num_cols = nt.select_dtypes(include=[np.number]).columns
        non_coord = [c for c in num_cols if c not in ["x", "y", "array_row", "array_col"]]
        if non_coord:
            nt_score_col = non_coord[0]
        else:
            return

    cond_col = None
    for candidate in ["condition", "group", "sample_condition"]:
        if candidate in nt.columns:
            cond_col = candidate
            break

    if len(coord_cols) < 2 or cond_col is None:
        # Load from spatial adata
        try:
            adata = load_deconvolved_adata()
        except FileNotFoundError:
            try:
                adata = load_spatial_adata()
            except FileNotFoundError:
                return

        shared_idx = nt.index.intersection(adata.obs_names)
        if len(shared_idx) == 0:
            return
        nt = nt.loc[shared_idx]
        nt["x"] = adata[shared_idx].obsm["spatial"][:, 0]
        nt["y"] = adata[shared_idx].obsm["spatial"][:, 1]
        if cond_col is None:
            nt["condition"] = adata.obs.loc[shared_idx, "condition"].values
            cond_col = "condition"

    conditions = sorted(nt[cond_col].unique())
    inner = gs.subgridspec(1, len(conditions), wspace=0.15)

    vals = nt[nt_score_col].values
    vmin, vmax = np.nanpercentile(vals, [2, 98])
    norm = Normalize(vmin=vmin, vmax=vmax)

    for i, cond in enumerate(conditions):
        ax = fig.add_subplot(inner[0, i])
        mask = nt[cond_col] == cond
        x = nt.loc[mask, "x"].values
        y = nt.loc[mask, "y"].values
        v = nt.loc[mask, nt_score_col].values

        sc_map = ax.scatter(
            x, y, c=v, cmap="viridis", s=1.2, norm=norm,
            edgecolors="none", rasterized=True,
        )
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(cond, fontsize=8, fontweight="bold")

        if i == len(conditions) - 1:
            _add_colorbar(fig, ax, sc_map, label="NT score", shrink=0.5, pad=0.05)


def panel_d_trajectory_shift(fig, gs):
    """Panel (d): Disease trajectory shift — KDE + cell-type composition.

    Left: NT score KDE overlay (Healthy vs Steatotic)
    Right: Stacked area of cell-type composition along trajectory
    """
    inner = gs.subgridspec(1, 2, wspace=0.35)

    # ── Left: KDE of NT scores ───────────────────────────────────────────
    ax_kde = fig.add_subplot(inner[0, 0])

    nt_path = RESULTS_DIR / "ontrac" / "niche_trajectory_scores.csv"
    disease_comp_path = RESULTS_DIR / "ontrac" / "disease_trajectory_comparison.csv"

    plotted_kde = False
    if disease_comp_path.exists():
        dc = pd.read_csv(disease_comp_path)
        cond_col = None
        nt_col = None
        for candidate in ["condition", "group"]:
            if candidate in dc.columns:
                cond_col = candidate
                break
        for candidate in ["NT_score", "nt_score", "niche_trajectory_score"]:
            if candidate in dc.columns:
                nt_col = candidate
                break
        if cond_col and nt_col:
            for cond in sorted(dc[cond_col].unique()):
                color = CONDITION_COLORS.get(cond, "#999999")
                vals = dc.loc[dc[cond_col] == cond, nt_col].dropna()
                if len(vals) > 10:
                    from scipy.stats import gaussian_kde
                    kde = gaussian_kde(vals, bw_method=0.3)
                    x_grid = np.linspace(vals.min(), vals.max(), 200)
                    ax_kde.plot(x_grid, kde(x_grid), lw=1.5, color=color, label=cond)
                    ax_kde.fill_between(x_grid, kde(x_grid), alpha=0.2, color=color)
                    plotted_kde = True
    elif nt_path.exists():
        nt = pd.read_csv(nt_path, index_col=0)
        cond_col = None
        nt_col = None
        for candidate in ["condition", "group"]:
            if candidate in nt.columns:
                cond_col = candidate
                break
        for candidate in ["NT_score", "nt_score", "niche_trajectory_score", "trajectory_score"]:
            if candidate in nt.columns:
                nt_col = candidate
                break
        if cond_col and nt_col:
            for cond in sorted(nt[cond_col].unique()):
                color = CONDITION_COLORS.get(cond, "#999999")
                vals = nt.loc[nt[cond_col] == cond, nt_col].dropna()
                if len(vals) > 10:
                    from scipy.stats import gaussian_kde
                    kde = gaussian_kde(vals, bw_method=0.3)
                    x_grid = np.linspace(vals.min(), vals.max(), 200)
                    ax_kde.plot(x_grid, kde(x_grid), lw=1.5, color=color, label=cond)
                    ax_kde.fill_between(x_grid, kde(x_grid), alpha=0.2, color=color)
                    plotted_kde = True

    if plotted_kde:
        ax_kde.set_xlabel("Niche trajectory score", fontsize=7)
        ax_kde.set_ylabel("Density", fontsize=7)
        ax_kde.set_title("NT score distribution", fontsize=8, fontweight="bold")
        ax_kde.legend(fontsize=6, frameon=True, framealpha=0.9, edgecolor="gray")
    else:
        ax_kde.text(0.5, 0.5, "NT score data\nnot available",
                    ha="center", va="center", fontsize=8,
                    transform=ax_kde.transAxes, color="gray")

    # ── Right: Stacked area of cell-type composition along trajectory ────
    ax_comp = fig.add_subplot(inner[0, 1])

    comp_path = RESULTS_DIR / "ontrac" / "trajectory_composition.csv"
    if not comp_path.exists():
        comp_path = RESULTS_DIR / "ontrac" / "niche_celltype_composition.csv"

    if comp_path.exists():
        comp = pd.read_csv(comp_path, index_col=0)
        # Expect columns = cell types, index = trajectory bins
        ct_cols = [c for c in comp.columns if c not in ["bin", "trajectory_bin", "NT_bin"]]
        if ct_cols:
            x_vals = np.arange(len(comp))
            # Sort cell types by mean abundance (largest at bottom)
            ct_order = comp[ct_cols].mean().sort_values(ascending=False).index.tolist()
            colors = plt.cm.Set3(np.linspace(0, 1, len(ct_order)))
            ax_comp.stackplot(
                x_vals,
                *[comp[ct].values for ct in ct_order],
                labels=ct_order,
                colors=colors,
                alpha=0.85,
            )
            ax_comp.set_xlabel("Trajectory position", fontsize=7)
            ax_comp.set_ylabel("Proportion", fontsize=7)
            ax_comp.set_title("Cell-type composition\nalong trajectory", fontsize=8, fontweight="bold")
            ax_comp.set_xlim(0, len(comp) - 1)
            ax_comp.set_ylim(0, 1)
            ax_comp.legend(
                fontsize=4.5, loc="center left", bbox_to_anchor=(1.02, 0.5),
                frameon=True, framealpha=0.9, edgecolor="gray",
                handletextpad=0.3, labelspacing=0.2,
            )
        else:
            ax_comp.text(0.5, 0.5, "No cell-type columns found",
                         ha="center", va="center", fontsize=8,
                         transform=ax_comp.transAxes, color="gray")
    else:
        ax_comp.text(0.5, 0.5, "trajectory_composition.csv\nnot found",
                     ha="center", va="center", fontsize=8,
                     transform=ax_comp.transAxes, color="gray")


def panel_e_tgfb_flow(fig, gs):
    """Panel (e): TGFB signaling flow — direction arrows on tissue.

    Shows quiver plot of TGFB signaling direction overlaid on spatial scatter,
    side-by-side Healthy vs Steatotic.
    """
    ax_container = fig.add_subplot(gs)
    ax_container.set_axis_off()

    commot_dir = RESULTS_DIR / "commot"

    # Try loading per-condition direction files
    direction_files = {}
    for cond in ["Healthy", "Steatotic", "MASLD"]:
        for pattern in [
            f"direction_TGFB_{cond}.csv",
            f"direction_TGFB1_{cond}.csv",
            f"commot_direction_TGFB_{cond}.csv",
        ]:
            path = commot_dir / pattern
            if path.exists():
                direction_files[cond] = pd.read_csv(path, index_col=0)
                break

    if not direction_files:
        # Try loading combined signaling results
        combined_path = commot_dir / "signaling_direction_vectors.csv"
        if combined_path.exists():
            combined = pd.read_csv(combined_path, index_col=0)
            # Filter to TGFB pathway
            pathway_col = None
            for candidate in ["pathway", "ligand", "lr_pair"]:
                if candidate in combined.columns:
                    pathway_col = candidate
                    break
            if pathway_col:
                tgfb_mask = combined[pathway_col].str.contains("TGFB", case=False, na=False)
                combined = combined[tgfb_mask]

            cond_col = None
            for candidate in ["condition", "group"]:
                if candidate in combined.columns:
                    cond_col = candidate
                    break
            if cond_col:
                for cond in combined[cond_col].unique():
                    direction_files[cond] = combined[combined[cond_col] == cond]

    if not direction_files:
        ax_container.text(
            0.5, 0.5, "TGFB direction data\nnot found",
            ha="center", va="center", fontsize=8,
            transform=ax_container.transAxes, color="gray",
        )
        return

    conditions = sorted(direction_files.keys())
    inner = gs.subgridspec(1, len(conditions), wspace=0.15)

    for i, cond in enumerate(conditions):
        ax = fig.add_subplot(inner[0, i])
        df = direction_files[cond]

        # Identify coordinate and direction columns
        x_col = "x" if "x" in df.columns else ("array_col" if "array_col" in df.columns else None)
        y_col = "y" if "y" in df.columns else ("array_row" if "array_row" in df.columns else None)
        dx_col = None
        dy_col = None
        for candidate in ["sender_vx", "dx", "direction_x", "u", "vx"]:
            if candidate in df.columns:
                dx_col = candidate
                break
        for candidate in ["sender_vy", "dy", "direction_y", "v", "vy"]:
            if candidate in df.columns:
                dy_col = candidate
                break

        if x_col is None or y_col is None or dx_col is None or dy_col is None:
            ax.text(0.5, 0.5, f"Missing columns\nfor {cond}",
                    ha="center", va="center", fontsize=7,
                    transform=ax.transAxes, color="gray")
            continue

        x = df[x_col].values
        y = df[y_col].values
        dx = df[dx_col].values
        dy = df[dy_col].values

        # Magnitude for coloring
        mag = np.sqrt(dx ** 2 + dy ** 2)
        mag_norm = Normalize(vmin=0, vmax=np.percentile(mag, 95))

        # Background scatter (all spots, gray)
        ax.scatter(x, y, s=0.5, c="#EEEEEE", edgecolors="none", rasterized=True, zorder=1)

        # Subsample for readability (max 500 arrows)
        n_arrows = min(500, len(df))
        if len(df) > n_arrows:
            idx = np.random.RandomState(42).choice(len(df), n_arrows, replace=False)
        else:
            idx = np.arange(len(df))

        quiver = ax.quiver(
            x[idx], y[idx], dx[idx], dy[idx],
            mag[idx], cmap="Reds", norm=mag_norm,
            scale_units="xy", angles="xy",
            width=0.003, headwidth=3, headlength=4,
            alpha=0.8, zorder=2,
        )
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(f"TGFB flow\n{cond}", fontsize=8, fontweight="bold")

        if i == len(conditions) - 1:
            _add_colorbar(fig, ax, quiver, label="Flow magnitude", shrink=0.5, pad=0.05)


def panel_f_upset(fig, gs):
    """Panel (f): Multi-method convergence UpSet plot.

    Shows intersection sizes for 4 methods: ONTraC TAG, gsMap risk,
    COMMOT regulated, MultiSP domain.
    """
    ax = fig.add_subplot(gs)

    conv_path = RESULTS_DIR / "spatial_integration" / "cross_method_convergence.csv"
    if not conv_path.exists():
        ax.text(
            0.5, 0.5, "cross_method_convergence.csv\nnot found\n(run 17a first)",
            ha="center", va="center", fontsize=8,
            transform=ax.transAxes, color="gray",
        )
        return

    conv = pd.read_csv(conv_path, index_col=0)

    method_cols = {
        "ONTraC": "ontrac_evidence",
        "gsMap": "gsmap_evidence",
        "COMMOT": "commot_evidence",
        "MultiSP": "multisp_evidence",
    }

    # Check which methods are present
    available = {name: col for name, col in method_cols.items() if col in conv.columns}
    if not available:
        ax.text(0.5, 0.5, "No method evidence\ncolumns found",
                ha="center", va="center", fontsize=8,
                transform=ax.transAxes, color="gray")
        return

    method_names = list(available.keys())
    method_col_list = list(available.values())
    n_methods = len(method_names)

    # Compute all intersection sizes
    # Each intersection is a binary tuple indicating membership
    from itertools import product
    intersections = {}
    for combo in product([False, True], repeat=n_methods):
        if not any(combo):
            continue  # skip the empty set
        mask = pd.Series(True, index=conv.index)
        for j, include in enumerate(combo):
            if include:
                mask &= conv[method_col_list[j]].fillna(False).astype(bool)
            else:
                mask &= ~conv[method_col_list[j]].fillna(False).astype(bool)
        count = mask.sum()
        if count > 0:
            intersections[combo] = count

    if not intersections:
        ax.text(0.5, 0.5, "No intersections\nwith evidence",
                ha="center", va="center", fontsize=8,
                transform=ax.transAxes, color="gray")
        return

    # Sort by cardinality (number of methods in intersection), then by size
    sorted_combos = sorted(
        intersections.keys(),
        key=lambda c: (-sum(c), -intersections[c]),
    )

    # Limit to top 15 intersections for readability
    sorted_combos = sorted_combos[:15]
    sizes = [intersections[c] for c in sorted_combos]
    n_bars = len(sorted_combos)

    # Create UpSet-style layout manually
    # Top: bar chart of intersection sizes
    # Bottom: dot matrix showing which methods are included
    ax.clear()

    bar_height = 0.7
    matrix_height = 0.3
    bar_bottom = matrix_height

    # Bar chart (top portion)
    bar_positions = np.arange(n_bars)
    n_sets_per_bar = [sum(c) for c in sorted_combos]
    bar_colors = []
    for n_s in n_sets_per_bar:
        if n_s >= 4:
            bar_colors.append("#B2182B")
        elif n_s == 3:
            bar_colors.append("#EF8A62")
        elif n_s == 2:
            bar_colors.append("#67A9CF")
        else:
            bar_colors.append("#D1D1D1")

    # Normalize positions to axes coordinates
    ax.set_xlim(-0.5, n_bars - 0.5)
    max_size = max(sizes) * 1.15
    ax.set_ylim(-n_methods - 0.5, max_size)

    # Draw bars
    ax.bar(bar_positions, sizes, color=bar_colors, edgecolor="none", width=0.6, zorder=3)

    # Add count labels on bars
    for x_pos, size in zip(bar_positions, sizes):
        ax.text(x_pos, size + max_size * 0.02, str(size), ha="center", va="bottom",
                fontsize=5, fontweight="bold")

    # Draw horizontal separator
    ax.axhline(0, color="gray", lw=0.5, zorder=2)

    # Draw dot matrix below
    for i, combo in enumerate(sorted_combos):
        active = [j for j, inc in enumerate(combo) if inc]
        inactive = [j for j, inc in enumerate(combo) if not inc]

        # Draw inactive dots (gray)
        for j in inactive:
            ax.plot(i, -(j + 0.5), "o", color="#D1D1D1", markersize=5, zorder=3)

        # Draw active dots (black) and connecting line
        for j in active:
            ax.plot(i, -(j + 0.5), "o", color="#333333", markersize=5, zorder=4)

        if len(active) > 1:
            y_range = [-(a + 0.5) for a in active]
            ax.plot([i, i], [min(y_range), max(y_range)], "-", color="#333333",
                    lw=1.5, zorder=3)

    # Method labels on y-axis (left of dot matrix)
    for j, name in enumerate(method_names):
        ax.text(-0.7, -(j + 0.5), name, ha="right", va="center", fontsize=6)

    ax.set_ylabel("Intersection size", fontsize=7)
    ax.set_title("Multi-method convergence", fontsize=8, fontweight="bold")
    ax.set_xticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)


# ═══════════════════════════════════════════════════════════════════════════════
# ASSEMBLY
# ═══════════════════════════════════════════════════════════════════════════════

def assemble_main_figure():
    """Assemble the 6-panel main spatial figure (2 rows x 3 columns)."""
    print_step("Assembling main spatial figure (2x3 layout)")

    # 183mm wide = Nature double column
    fig = plt.figure(figsize=(183 / 25.4, 140 / 25.4))
    gs = GridSpec(
        2, 3, figure=fig,
        hspace=0.45, wspace=0.35,
        height_ratios=[1, 1],
        width_ratios=[1, 1, 1],
    )

    # Row 1
    print_step("Panel (a): gsMap risk heatmap")
    panel_a_gsmap_risk(fig, gs[0, 0])

    print_step("Panel (b): GWAS zonation enrichment")
    panel_b_gwas_zonation(fig, gs[0, 1])

    print_step("Panel (c): Niche trajectory maps")
    panel_c_niche_trajectory(fig, gs[0, 2])

    # Row 2
    print_step("Panel (d): Trajectory shift")
    panel_d_trajectory_shift(fig, gs[1, 0])

    print_step("Panel (e): TGFB signaling flow")
    panel_e_tgfb_flow(fig, gs[1, 1])

    print_step("Panel (f): Multi-method UpSet")
    panel_f_upset(fig, gs[1, 2])

    # Panel labels
    label_positions = [
        ("a", 0.01, 0.98), ("b", 0.34, 0.98), ("c", 0.67, 0.98),
        ("d", 0.01, 0.48), ("e", 0.34, 0.48), ("f", 0.67, 0.48),
    ]
    for label, x, y in label_positions:
        fig.text(x, y, label, fontsize=12, fontweight="bold", va="top",
                 fontfamily="sans-serif")

    save_fig(fig, "fig_spatial_main")


def save_individual_panels():
    """Save each panel as a separate PDF for flexible layout."""
    panel_funcs = [
        ("fig_spatial_panel_a", panel_a_gsmap_risk, (89, 65)),
        ("fig_spatial_panel_b", panel_b_gwas_zonation, (89, 65)),
        ("fig_spatial_panel_c", panel_c_niche_trajectory, (89, 65)),
        ("fig_spatial_panel_d", panel_d_trajectory_shift, (89, 65)),
        ("fig_spatial_panel_e", panel_e_tgfb_flow, (89, 65)),
        ("fig_spatial_panel_f", panel_f_upset, (89, 65)),
    ]

    for name, func, (w_mm, h_mm) in panel_funcs:
        print_step(f"Individual panel: {name}")
        fig = plt.figure(figsize=(w_mm / 25.4, h_mm / 25.4))
        gs = GridSpec(1, 1, figure=fig)
        func(fig, gs[0, 0])
        save_fig(fig, name)


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    print_header("17b: Spatial Main Figure")

    config = load_config()

    # Assemble combined figure
    assemble_main_figure()

    # Save individual panels
    save_individual_panels()

    print_header("17b: Complete")


if __name__ == "__main__":
    main()

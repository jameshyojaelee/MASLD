#!/usr/bin/env python3
"""
Supplementary ATAC-seq / scATAC-seq panels (Python-based).
Generates panels that require h5ad data: QC metrics, UMAPs, cell type
proportions, marker gene activity, label transfer confusion matrix.

Output: figures/supplementary/figS05_epigenomic_spatial/panel_*.pdf
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
ATAC = os.path.join(BASE, "Analysis", "ATAC", "Human_Multiome")
SNAP_DIR = os.path.join(ATAC, "results", "snapatac2")
LABEL_TRANSFER_DIR = os.path.join(ATAC, "results", "label_transfer")
PER_DONOR = os.path.join(SNAP_DIR, "per_donor")
FRAG_DIR = os.path.join(ATAC, "results", "fragments")

def _get_processed_h5ad():
    """Return path to processed h5ad, preferring label-transferred version."""
    lt = os.path.join(LABEL_TRANSFER_DIR, "snapatac2_label_transferred.h5ad")
    if os.path.exists(lt):
        return lt
    return os.path.join(SNAP_DIR, "snapatac2_processed.h5ad")
OUT = os.path.join(BASE, "figures", "supplementary", "figS05_epigenomic_spatial")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------------------
# Style — match publication_theme.R magenta/blue palette
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "Helvetica",
    "font.size": 6,
    "axes.titlesize": 6,
    "axes.labelsize": 6,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "pdf.fonttype": 42,       # Illustrator-editable
    "ps.fonttype": 42,
    "axes.linewidth": 0.3,
    "xtick.major.width": 0.3,
    "ytick.major.width": 0.3,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

MASLD_UP = "#C9265E"
MASLD_DN = "#1565C0"
MASLD_NS = "#BDBDBD"
MASLD_CTRL = "#9E9E9E"
MASLD_DIS = "#C9265E"

# Cell type palette — label-transfer-corrected (CellTypist names)
# Fallback keys for gene-activity names retained for backwards compatibility
CT_PALETTE = {
    "Hepatocytes":           "#0D47A1",
    "Fibroblasts":           "#F57F17",
    "Macrophages":           "#C2185B",
    "Plasma cells":          "#5D4037",
    "Endothelial cells":     "#2E7D32",
    "Cholangiocytes":        "#1565C0",
    "T cells":               "#7B1FA2",
    "Circulating NK/NKT":    "#9C27B0",
    "Resident NK":           "#AB47BC",
    "Low_confidence":        "#BDBDBD",
    # Gene-activity names (for backwards compatibility with old h5ad)
    "Hepatocyte":            "#0D47A1",
    "Endothelial":           "#2E7D32",
    "Stellate_Cell":         "#F57F17",
    "Plasma_Cell":           "#5D4037",
    "Macrophage":            "#C2185B",
    "LSEC":                  "#00695C",
    "Cholangiocyte":         "#1565C0",
    "Kupffer_Cell":          "#E91E63",
    "NK_T_Cell":             "#7B1FA2",
    "B_Cell":                "#9C27B0",
}

# Width in inches matching publication_theme.R
FIG_FULL = 180 / 25.4
FIG_HALF = 88 / 25.4


def save_panel(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ===================================================================
# Step 1: Load per-donor h5ad files for QC metrics
# ===================================================================
def load_per_donor_qc():
    """Load QC metrics from all per-donor h5ad files."""
    import anndata as ad

    records = []
    for fname in sorted(os.listdir(PER_DONOR)):
        if not fname.endswith(".h5ad"):
            continue
        donor_id = fname.replace(".h5ad", "")
        path = os.path.join(PER_DONOR, fname)
        adata = ad.read_h5ad(path, backed="r")
        obs = adata.obs.copy()
        obs["donor_id"] = donor_id
        records.append(obs)
        adata.file.close()
    df = pd.concat(records, ignore_index=True)
    return df


# ===================================================================
# Step 2: Load main processed h5ad for UMAP + cell type
# ===================================================================
def load_main_adata():
    """Load the main processed scATAC h5ad (backed) for UMAP coords + obs."""
    import anndata as ad

    path = _get_processed_h5ad()
    adata = ad.read_h5ad(path, backed="r")
    umap = adata.obsm["X_umap"].copy()
    # Prefer label-transferred cell types if available
    ct_col = "cell_type_transferred" if "cell_type_transferred" in adata.obs.columns else "cell_type"
    obs = adata.obs[[ct_col, "condition", "donor_id"]].copy()
    obs.rename(columns={ct_col: "cell_type"}, inplace=True)
    adata.file.close()
    obs["UMAP1"] = umap[:, 0]
    obs["UMAP2"] = umap[:, 1]
    return obs


# ===================================================================
# Panel 1: Fragment size distribution (nucleosome banding)
# ===================================================================
def panel_01_fragment_size(qc_df):
    """Sample fragment sizes from fragment files to show nucleosome banding."""
    import gzip

    print("Panel 01: Fragment size distribution...")
    frag_sizes = []
    frag_files = sorted(
        [f for f in os.listdir(FRAG_DIR) if f.endswith("_fragments.tsv.gz")]
    )

    if not frag_files:
        print("  WARNING: No fragment files found, skipping panel 01")
        return

    # Sample from first 6 donors for speed
    for fname in frag_files[:6]:
        path = os.path.join(FRAG_DIR, fname)
        count = 0
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    size = int(parts[2]) - int(parts[1])
                    if 0 < size < 1000:
                        frag_sizes.append(size)
                count += 1
                if count >= 500000:  # Sample 500K fragments per donor
                    break

    if not frag_sizes:
        print("  WARNING: No fragments parsed, skipping panel 01")
        return

    fig, ax = plt.subplots(figsize=(FIG_HALF, 2.5))
    ax.hist(frag_sizes, bins=np.arange(0, 1001, 5), color=MASLD_DN,
            alpha=0.8, density=True, linewidth=0)
    ax.set_xlabel("Fragment size (bp)")
    ax.set_ylabel("Density")
    print("  [caption] Fragment size distribution")
    ax.set_xlim(0, 1000)

    # Annotate nucleosome peaks
    for pos, label in [(150, "NFR"), (200, "Mono"), (400, "Di"), (600, "Tri")]:
        ax.axvline(pos, color=MASLD_NS, linestyle="--", linewidth=0.5, alpha=0.6)
        ax.text(pos, ax.get_ylim()[1] * 0.95, label, ha="center", fontsize=6,
                color="#808080")

    save_panel(fig, "panel_01_fragment_size.pdf")


# ===================================================================
# Panel 2: TSS enrichment score distribution
# ===================================================================
def panel_02_tss_enrichment(qc_df):
    print("Panel 02: TSS enrichment distribution...")
    if "tsse" not in qc_df.columns:
        print("  WARNING: tsse not in QC data, skipping")
        return

    fig, ax = plt.subplots(figsize=(FIG_HALF, 2.5))
    tsse = qc_df["tsse"].dropna()
    ax.hist(tsse, bins=100, color=MASLD_DN, alpha=0.8, edgecolor="none")
    ax.axvline(5, color=MASLD_UP, linestyle="--", linewidth=0.8, label="QC threshold (5)")
    ax.set_xlabel("TSS enrichment score")
    ax.set_ylabel("Number of cells")
    print("  [caption] TSS enrichment score distribution")
    ax.legend(frameon=False)

    n_pass = (tsse >= 5).sum()
    n_total = len(tsse)
    ax.text(0.95, 0.85, f"{n_pass:,}/{n_total:,} pass ({100*n_pass/n_total:.1f}%)",
            transform=ax.transAxes, ha="right", fontsize=6, color="#4D4D4D")

    save_panel(fig, "panel_02_tss_enrichment.pdf")


# ===================================================================
# Panel 3: QC scatter — n_fragment vs TSS enrichment
# ===================================================================
def panel_03_qc_scatter(qc_df):
    print("Panel 03: QC scatter...")
    if "tsse" not in qc_df.columns or "n_fragment" not in qc_df.columns:
        print("  WARNING: Missing tsse or n_fragment, skipping")
        return

    fig, ax = plt.subplots(figsize=(FIG_HALF, FIG_HALF * 0.8))
    # Subsample for plotting
    n = min(30000, len(qc_df))
    sub = qc_df.sample(n=n, random_state=42)

    ax.scatter(sub["n_fragment"], sub["tsse"], s=0.3, alpha=0.3,
               c=MASLD_DN, rasterized=True)
    ax.axhline(5, color=MASLD_UP, linestyle="--", linewidth=0.5, alpha=0.7)
    ax.set_xlabel("Number of fragments")
    ax.set_ylabel("TSS enrichment score")
    print("  [caption] QC: fragments vs TSS enrichment")
    ax.set_xscale("log")

    save_panel(fig, "panel_03_qc_scatter.pdf")


# ===================================================================
# Panel 4: UMAP by cell type
# ===================================================================
def panel_04_umap_celltype(obs_df):
    print("Panel 04: UMAP by cell type...")
    fig, ax = plt.subplots(figsize=(FIG_HALF + 1, FIG_HALF * 0.8))

    # Shuffle for even overplotting
    obs_df = obs_df.sample(frac=1, random_state=42)

    for ct in CT_PALETTE:
        mask = obs_df["cell_type"] == ct
        if mask.sum() == 0:
            continue
        ax.scatter(
            obs_df.loc[mask, "UMAP1"], obs_df.loc[mask, "UMAP2"],
            s=0.1, alpha=0.3, c=CT_PALETTE[ct], label=f"{ct} ({mask.sum():,})",
            rasterized=True,
        )

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    print("  [caption] scATAC-seq: cell types (88,814 cells)")
    ax.legend(markerscale=8, fontsize=6, frameon=False, loc="center left",
              bbox_to_anchor=(1.0, 0.5), handletextpad=0.3)

    save_panel(fig, "panel_04_umap_celltype.pdf")


# ===================================================================
# Panel 5: UMAP by donor
# ===================================================================
def panel_05_umap_donor(obs_df):
    print("Panel 05: UMAP by donor...")
    fig, ax = plt.subplots(figsize=(FIG_HALF + 1, FIG_HALF * 0.8))

    donors = sorted(obs_df["donor_id"].unique())
    cmap = plt.cm.get_cmap("tab20", len(donors))
    obs_df = obs_df.sample(frac=1, random_state=42)

    for i, donor in enumerate(donors):
        mask = obs_df["donor_id"] == donor
        ax.scatter(
            obs_df.loc[mask, "UMAP1"], obs_df.loc[mask, "UMAP2"],
            s=0.1, alpha=0.3, c=[cmap(i)], label=donor, rasterized=True,
        )

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    print("  [caption] scATAC-seq: donors (batch mixing)")
    ax.legend(markerscale=8, fontsize=6, frameon=False, loc="center left",
              bbox_to_anchor=(1.0, 0.5), ncol=1, handletextpad=0.3)

    save_panel(fig, "panel_05_umap_donor.pdf")


# ===================================================================
# Panel 6: Cell type proportions by condition (MASLD vs Normal)
# ===================================================================
def panel_06_celltype_proportions(obs_df):
    print("Panel 06: Cell type proportions by condition...")
    ct_counts = (
        obs_df.groupby(["donor_id", "condition", "cell_type"])
        .size()
        .reset_index(name="count")
    )
    # Compute proportions per donor
    totals = ct_counts.groupby("donor_id")["count"].transform("sum")
    ct_counts["proportion"] = ct_counts["count"] / totals

    # Average proportions per condition × cell type
    ct_avg = (
        ct_counts.groupby(["condition", "cell_type"])["proportion"]
        .agg(["mean", "sem"])
        .reset_index()
    )

    fig, ax = plt.subplots(figsize=(FIG_HALF, 3))
    cell_types = list(CT_PALETTE.keys())
    conditions = sorted(ct_avg["condition"].unique())
    x = np.arange(len(cell_types))
    n_cond = len(conditions)
    width = 0.8 / n_cond
    # Assign colors: NORMAL=blue, MASL=pink, MASH=magenta
    cond_colors = {}
    for c in conditions:
        c_upper = c.upper()
        if "NORMAL" in c_upper or "CONTROL" in c_upper or "HEALTHY" in c_upper:
            cond_colors[c] = MASLD_CTRL
        elif "MASH" in c_upper or "NASH" in c_upper:
            cond_colors[c] = MASLD_DIS
        else:
            cond_colors[c] = "#F4A674"  # MASL/intermediate = Liang peach

    for i, cond in enumerate(conditions):
        vals = []
        errs = []
        for ct in cell_types:
            row = ct_avg[(ct_avg["condition"] == cond) & (ct_avg["cell_type"] == ct)]
            vals.append(row["mean"].values[0] if len(row) > 0 else 0)
            errs.append(row["sem"].values[0] if len(row) > 0 else 0)
        ax.bar(x + (i - n_cond/2 + 0.5) * width, vals, width, yerr=errs, label=cond,
               color=cond_colors.get(cond, MASLD_NS), edgecolor="none",
               capsize=2, error_kw={"linewidth": 0.5})

    ax.set_xticks(x)
    ax.set_xticklabels([ct.replace("_", "\n") for ct in cell_types],
                       rotation=45, ha="right", fontsize=6)
    ax.set_ylabel("Proportion")
    print("  [caption] Cell type proportions by condition")
    ax.legend(frameon=False)

    save_panel(fig, "panel_06_celltype_proportions.pdf")


# ===================================================================
# Panel 7: Marker gene activity dot plot
# ===================================================================
def panel_07_marker_dotplot(obs_df):
    """Dot plot of canonical liver marker gene activity scores by cell type."""
    import anndata as ad

    print("Panel 07: Marker gene activity dot plot...")
    ga_path = os.path.join(SNAP_DIR, "gene_activity_matrix.h5ad")
    if not os.path.exists(ga_path):
        print("  WARNING: gene_activity_matrix.h5ad not found, skipping")
        return

    # Marker genes for liver cell types (CellTypist label-transfer names)
    markers = {
        "Hepatocytes": ["ALB", "APOB", "CYP3A4", "HNF4A"],
        "Cholangiocytes": ["KRT19", "KRT7", "EPCAM"],
        "Endothelial cells": ["PECAM1", "CDH5", "VWF"],
        "Fibroblasts": ["ACTA2", "COL1A1", "DCN"],
        "Macrophages": ["CD68", "MARCO", "CLEC4F"],
        "T cells": ["CD3D", "NKG7"],
        "Plasma cells": ["JCHAIN", "MZB1"],
    }
    all_genes = []
    for gl in markers.values():
        all_genes.extend(gl)
    all_genes = list(dict.fromkeys(all_genes))  # deduplicate preserving order

    # Load into memory (backed mode has indexing bugs in this anndata version)
    adata = ad.read_h5ad(ga_path)
    var_names = list(adata.var_names)

    # Map gene names to indices
    gene_idx = {}
    for g in all_genes:
        if g in var_names:
            gene_idx[g] = var_names.index(g)

    if not gene_idx:
        print("  WARNING: No marker genes found in gene activity matrix, skipping")
        del adata
        return

    # Get cell type labels from main adata (prefer label-transferred)
    main_path = _get_processed_h5ad()
    main_adata = ad.read_h5ad(main_path, backed="r")
    ct_col = "cell_type_transferred" if "cell_type_transferred" in main_adata.obs.columns else "cell_type"
    cell_types_series = main_adata.obs[ct_col].copy()
    # Verify alignment: gene activity and processed h5ad must share obs_names
    shared = adata.obs_names.intersection(main_adata.obs_names)
    if len(shared) < 0.9 * min(adata.n_obs, main_adata.n_obs):
        print(f"  WARNING: obs_name overlap only {len(shared)}/{adata.n_obs} — using positional alignment")
    try:
        main_adata.file.close()
    except Exception:
        pass

    # Read gene activity for selected genes (slice columns)
    # For backed sparse, we need to read in chunks
    found_genes = list(gene_idx.keys())
    n_cells = adata.shape[0]

    # Read full matrix for selected genes only
    from scipy.sparse import issparse

    # Extract marker gene columns (in-memory, no chunking needed)
    expr_data = {}
    col_indices = [gene_idx[g] for g in found_genes]
    X = adata.X
    if issparse(X):
        mat = X[:, col_indices].toarray().astype(np.float32)
    else:
        mat = np.asarray(X[:, col_indices], dtype=np.float32)
    for i, g in enumerate(found_genes):
        expr_data[g] = mat[:, i]
    del mat, adata

    expr_df = pd.DataFrame(expr_data, index=range(n_cells))
    expr_df["cell_type"] = cell_types_series.values[:n_cells]

    # Compute mean expression and fraction expressed per cell type per gene
    cts = list(CT_PALETTE.keys())
    dot_data = []
    for ct in cts:
        mask = expr_df["cell_type"] == ct
        if mask.sum() == 0:
            continue
        sub = expr_df.loc[mask, found_genes]
        for g in found_genes:
            vals = sub[g].values
            mean_expr = np.mean(vals)
            frac_expr = np.mean(vals > 0)
            dot_data.append({
                "cell_type": ct, "gene": g,
                "mean_expr": mean_expr, "frac_expr": frac_expr,
            })

    dot_df = pd.DataFrame(dot_data)
    if dot_df.empty:
        print("  WARNING: No dot plot data, skipping")
        return

    # Normalize mean expression per gene for color scaling
    for g in found_genes:
        mask = dot_df["gene"] == g
        vals = dot_df.loc[mask, "mean_expr"]
        vmax = vals.max()
        if vmax > 0:
            dot_df.loc[mask, "scaled_expr"] = vals / vmax
        else:
            dot_df.loc[mask, "scaled_expr"] = 0

    # Plot
    fig, ax = plt.subplots(figsize=(FIG_FULL * 0.7, FIG_HALF * 0.7))

    # Create grid
    gene_order = found_genes
    ct_order = cts

    for i, ct in enumerate(ct_order):
        for j, g in enumerate(gene_order):
            row = dot_df[(dot_df["cell_type"] == ct) & (dot_df["gene"] == g)]
            if row.empty:
                continue
            frac = row["frac_expr"].values[0]
            scaled = row["scaled_expr"].values[0]
            size = frac * 200  # Scale dot size
            color_val = plt.cm.Reds(scaled * 0.8 + 0.1)
            ax.scatter(j, i, s=size, c=[color_val], edgecolors="#808080",
                       linewidths=0.3)

    ax.set_xticks(range(len(gene_order)))
    ax.set_xticklabels(gene_order, rotation=90, ha="center", fontsize=6)
    ax.set_yticks(range(len(ct_order)))
    ax.set_yticklabels([ct.replace("_", " ") for ct in ct_order], fontsize=6)
    print("  [caption] Gene activity scores: canonical liver markers")
    ax.set_xlim(-0.5, len(gene_order) - 0.5)
    ax.set_ylim(-0.5, len(ct_order) - 0.5)

    # Size legend
    for frac_val, label in [(0.25, "25%"), (0.5, "50%"), (0.75, "75%")]:
        ax.scatter([], [], s=frac_val * 200, c="#808080", edgecolors="#808080",
                   linewidths=0.3, label=label)
    ax.legend(title="Frac. expr.", frameon=False, fontsize=6,
              title_fontsize=6, loc="upper left", bbox_to_anchor=(1.01, 1))

    save_panel(fig, "panel_07_marker_dotplot.pdf")


# ===================================================================
# Panel 8: Label transfer confusion-style bar chart
# ===================================================================
def panel_08_label_transfer(obs_df):
    """
    Show cell type composition from label transfer.
    Since no separate confusion matrix CSV exists, show per-donor cell type
    breakdown to validate label transfer consistency across donors.
    """
    print("Panel 08: Label transfer cell type assignment per donor...")

    ct_counts = (
        obs_df.groupby(["donor_id", "cell_type"])
        .size()
        .reset_index(name="count")
    )
    totals = ct_counts.groupby("donor_id")["count"].transform("sum")
    ct_counts["proportion"] = ct_counts["count"] / totals

    # Pivot for stacked bar
    pivot = ct_counts.pivot_table(
        index="donor_id", columns="cell_type", values="proportion", fill_value=0
    )
    # Order cell types by overall abundance
    ct_order = [ct for ct in CT_PALETTE if ct in pivot.columns]
    pivot = pivot[ct_order]

    fig, ax = plt.subplots(figsize=(FIG_HALF, 3))
    bottom = np.zeros(len(pivot))
    donors = pivot.index.tolist()
    x = np.arange(len(donors))

    for ct in ct_order:
        vals = pivot[ct].values
        color = CT_PALETTE.get(ct, MASLD_NS)
        ax.bar(x, vals, bottom=bottom, color=color, edgecolor="none",
               width=0.8, label=ct.replace("_", " "))
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(donors, rotation=45, ha="right", fontsize=6)
    ax.set_ylabel("Proportion")
    ax.set_xlabel("Donor")
    print("  [caption] Cell type composition per donor (label transfer)")
    ax.legend(fontsize=6, frameon=False, loc="center left",
              bbox_to_anchor=(1.0, 0.5), ncol=1)
    ax.set_ylim(0, 1)

    save_panel(fig, "panel_08_label_transfer.pdf")


# ===================================================================
# Main
# ===================================================================
def main():
    print("=" * 60)
    print("Generating supplementary ATAC panels (Python)")
    print(f"Output: {OUT}")
    print("=" * 60)

    # Load QC data from per-donor h5ad files
    print("\nLoading per-donor QC metrics...")
    qc_df = load_per_donor_qc()
    print(f"  {len(qc_df):,} cells from {qc_df['donor_id'].nunique()} donors")
    print(f"  QC columns: {list(qc_df.columns)}")

    # Load main adata for UMAPs
    print("\nLoading main processed h5ad for UMAP coordinates...")
    obs_df = load_main_adata()
    print(f"  {len(obs_df):,} cells, {obs_df['cell_type'].nunique()} cell types")

    # Generate panels
    panel_01_fragment_size(qc_df)
    panel_02_tss_enrichment(qc_df)
    panel_03_qc_scatter(qc_df)
    panel_04_umap_celltype(obs_df)
    panel_05_umap_donor(obs_df)
    panel_06_celltype_proportions(obs_df)
    panel_07_marker_dotplot(obs_df)
    panel_08_label_transfer(obs_df)

    print("\n" + "=" * 60)
    print("Python panels complete.")
    print("=" * 60)


if __name__ == "__main__":
    main()

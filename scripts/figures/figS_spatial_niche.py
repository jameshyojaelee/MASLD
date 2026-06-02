#!/usr/bin/env python3
"""figS_spatial_niche.py — Squidpy nhood_enrichment z-score matrix, Healthy vs Steatotic.
Idiom adapted from Elison/Ren/Gaulton 2025 Fig 5L. GSE192741 Visium (n=6,546 spots),
cell2location 16-cell-type deconvolution collapsed to 11 functional groups; per-spot
label = argmax of column-z-scored group-summed abundance (recovers rare cell types
that winner-take-all hides under hepatocyte dominance).
"""
import pathlib
import sys
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import squidpy as sq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import zscore

BASE = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
H5AD = BASE / "Analysis/Spatial/results/cell2location/spatial_model/spatial_deconvolved.h5ad"
OUT_DIR = BASE / "figures/supplementary/figS_spatial_niche"
OUT_PDF = OUT_DIR / "figS_spatial_niche.pdf"
DATA_OUT = BASE / "Analysis/Spatial/results/communication"
PREFIX = "c2l_q05cell_abundance_w_sf_means_per_cluster_mu_fg_"

# Collapse the c2l 16-type label space to 11 functional groups.
# Singletons (kept as-is): Hepatocytes, Cholangiocytes, Endothelial, Fibroblasts,
# Macrophages, T cells. Grouped: Monocytes (mono+derived), Granulocytes (neutro+baso),
# B/Plasma, NK (resident+circulating/NKT), DCs (cDC1+cDC2+pDC).
GROUPS = {
    "Hepatocytes":         ["Hepatocytes"],
    "Cholangiocytes":      ["Cholangiocytes"],
    "Endothelial cells":   ["Endothelial cells"],
    "Fibroblasts":         ["Fibroblasts"],
    "Macrophages":         ["Macrophages"],
    "Monocytes":           ["Mono+mono derived cells"],
    "Granulocytes":        ["Neutrophils", "Basophils"],
    "T cells":             ["T cells"],
    "B/Plasma":            ["B cells", "Plasma cells"],
    "NK":                  ["Resident NK", "Circulating NK/NKT"],
    "DCs":                 ["cDC1s", "cDC2s", "pDCs"],
}
ORDER = list(GROUPS.keys())


def assign_argmax_z(adata):
    """Per-spot label = argmax of column-z-scored group-summed c2l abundances."""
    c2l_cols = [c for c in adata.obs.columns if c.startswith("c2l_q05")]
    short = {c: c.replace(PREFIX, "") for c in c2l_cols}
    df = adata.obs[c2l_cols].astype(float).rename(columns=short)
    grouped = pd.DataFrame(
        {g: df[[m for m in members if m in df.columns]].sum(axis=1)
         for g, members in GROUPS.items()},
        index=df.index,
    )
    z = grouped.apply(zscore, axis=0)
    label = z.idxmax(axis=1)
    adata.obs["cell_type_argmaxz"] = pd.Categorical(label.values, categories=ORDER)


def _nhood_z_one_slice(sub, n_neighs=6, n_perms=1000):
    """nhood_enrichment z-matrix for a SINGLE donor slice (graph built within it)."""
    sub = sub.copy()
    sub.obs["cell_type_argmaxz"] = sub.obs["cell_type_argmaxz"].cat.remove_unused_categories()
    # Build the spatial graph within this one slide only: spots from different
    # donors reuse the same array_row/array_col frame, so a pooled graph would
    # wire physically distinct slides together (F210). One slice -> no cross-donor
    # edges.
    sq.gr.spatial_neighbors(sub, coord_type="generic", n_neighs=n_neighs)
    sq.gr.nhood_enrichment(sub, cluster_key="cell_type_argmaxz", n_perms=n_perms, seed=42)
    res = sub.uns["cell_type_argmaxz_nhood_enrichment"]
    cats = list(sub.obs["cell_type_argmaxz"].cat.categories)
    z = pd.DataFrame(res["zscore"], index=cats, columns=cats)
    return z.reindex(index=ORDER, columns=ORDER)


def compute_zscore_matrix(adata, condition, donor_col="sample_id",
                          n_neighs=6, n_perms=1000, min_spots=50):
    """Donor-aware nhood_enrichment: compute the z-matrix PER donor slice, then
    average across the 2-3 donors in the condition (F208/F210).

    The previous version pooled all spots of a condition into one graph and ran a
    single permutation z-score over thousands of *spots* -- but the cohort is only
    2-3 donors per arm and spots within a slice are spatially autocorrelated, so
    that n was pseudoreplicated. Here the unit of inference is the donor: one
    z-matrix per donor, then the cell-wise mean. With n=2-3 donors this panel is
    descriptive (see caption), which is the honest consequence of the fix.
    """
    cond = adata[adata.obs["condition"] == condition]
    donors = [d for d in cond.obs[donor_col].astype(str).unique()]
    per_donor = []
    for d in donors:
        sub = cond[cond.obs[donor_col].astype(str) == d]
        if sub.n_obs < min_spots:
            print(f"    skip donor {d} ({sub.n_obs} spots < {min_spots})")
            continue
        try:
            per_donor.append(_nhood_z_one_slice(sub, n_neighs=n_neighs, n_perms=n_perms))
        except Exception as e:
            print(f"    donor {d} nhood failed: {str(e)[:80]}")
    if not per_donor:
        return pd.DataFrame(np.nan, index=ORDER, columns=ORDER)
    # Cell-wise mean across donors (NaN where a cell-type pair is absent in a donor).
    stacked = np.stack([m.values for m in per_donor])
    mean_z = np.nanmean(stacked, axis=0)
    print(f"    {condition}: averaged nhood z over {len(per_donor)} donor(s)")
    return pd.DataFrame(mean_z, index=ORDER, columns=ORDER)


def plot_two_panel(z_h, z_d, n_h, n_d, vlim=8):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    norm = TwoSlopeNorm(vmin=-vlim, vcenter=0, vmax=vlim)
    cmap = "RdBu_r"
    # Report DONOR n (the unit of inference), not pooled spot n (F208). n_h/n_d
    # are donor counts; the z-matrices are per-donor means (descriptive at n=2-3).
    titles = [f"Healthy (n={n_h} donors)", f"Steatotic (n={n_d} donors)"]
    for ax, mat, title in zip(axes, [z_h, z_d], titles):
        im = ax.imshow(mat.values, cmap=cmap, norm=norm, aspect="equal")
        ax.set_xticks(range(len(ORDER)))
        ax.set_yticks(range(len(ORDER)))
        ax.set_xticklabels(ORDER, rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(ORDER, fontsize=8)
        ax.set_title(title, fontsize=10)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
    cbar = fig.colorbar(im, ax=axes, shrink=0.7, label="Squidpy neighborhood enrichment z-score")
    cbar.outline.set_visible(False)
    fig.suptitle("Spatial cell-type neighborhood enrichment (GSE192741 Visium)\n"
                 "per-donor z-matrices averaged across donors (descriptive; n=2-3 donors/arm)",
                 fontsize=10)
    fig.savefig(OUT_PDF, format="pdf", bbox_inches="tight")
    plt.close(fig)


def top_pairs(z, k=5):
    cells = z.index.tolist()
    rows = []
    for i, a in enumerate(cells):
        for j, b in enumerate(cells):
            if i < j and pd.notna(z.loc[a, b]):
                rows.append((a, b, float(z.loc[a, b])))
    rows.sort(key=lambda r: r[2], reverse=True)
    return rows[:k]


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DATA_OUT.mkdir(parents=True, exist_ok=True)
    print(f"Loading {H5AD}")
    adata = ad.read_h5ad(H5AD)
    assign_argmax_z(adata)
    print("Argmax-z label counts by condition:")
    print(pd.crosstab(adata.obs["cell_type_argmaxz"], adata.obs["condition"]))

    # Donor (not spot) counts per arm — the unit of inference (F208).
    n_h = int(adata.obs.loc[adata.obs["condition"] == "Healthy", "sample_id"].nunique())
    n_d = int(adata.obs.loc[adata.obs["condition"] == "Steatotic", "sample_id"].nunique())
    z_h = compute_zscore_matrix(adata, "Healthy")
    z_d = compute_zscore_matrix(adata, "Steatotic")
    z_h.to_csv(DATA_OUT / "nhood_enrichment_grouped_Healthy.csv")
    z_d.to_csv(DATA_OUT / "nhood_enrichment_grouped_Steatotic.csv")
    delta = z_d - z_h
    delta.to_csv(DATA_OUT / "nhood_enrichment_grouped_delta.csv")

    plot_two_panel(z_h, z_d, n_h, n_d)
    print(f"Wrote {OUT_PDF}")

    print("Top 5 disease-enriched pairs (delta z, Steatotic - Healthy):")
    cells = delta.index.tolist()
    pairs = []
    for i, a in enumerate(cells):
        for j, b in enumerate(cells):
            if i < j and pd.notna(delta.loc[a, b]):
                pairs.append((a, b, float(delta.loc[a, b])))
    pairs.sort(key=lambda r: r[2], reverse=True)
    for a, b, d in pairs[:5]:
        print(f"  {a} <-> {b}: delta_z = {d:+.2f}  (H={z_h.loc[a,b]:+.2f}  D={z_d.loc[a,b]:+.2f})")


if __name__ == "__main__":
    main()

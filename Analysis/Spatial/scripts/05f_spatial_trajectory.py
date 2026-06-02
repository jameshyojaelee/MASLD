#!/usr/bin/env python3
"""
05f_spatial_trajectory.py — Disease progression trajectory in spatial context.

Computes diffusion pseudotime on spatial data, mapping disease progression
onto tissue coordinates. Identifies spatial transition zones.

Runs on the deconvolved GSE192741 'Steatotic' spots (3 slices / 2 donors:
JBO014/JBO015 = donor H35, JBO019 = donor H37). There is NO fibrosis-stage
spectrum in this input, so the trajectory is a within-Steatotic descriptive
illustration at n~2 donors, not a staged progression. (The earlier
'HRA007511 (MASLD spectrum)' note was wrong — HRA007511 is permanently blocked
[HMSMA embargo] and was never processed.)

SLURM: --partition=cpu --cpus=8 --mem=64G --time=4:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, load_deconvolved_adata,
    save_checkpoint, save_csv, print_header,
)
from spatial_stats import spearman_by_donor

# Donor (slide) column. The trajectory runs on GSE192741 'Steatotic' spots =
# 3 slices / 2 donors (JBO014+JBO015 = donor H35, JBO019 = donor H37). Spots
# are spatially autocorrelated and share donors, so the effective n is ~2
# donors, NOT thousands of spots — every test here is donor-pseudoreplicated
# and is reported as descriptive only (F063).
DONOR_COL = "sample_id"


def compute_spatial_pseudotime(adata):
    """Compute diffusion pseudotime in spatial context."""
    # PCA + neighbors + diffusion map
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3")

    sc.pp.pca(adata, n_comps=30, use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.diffmap(adata, n_comps=10)

    # Set root: spots with lowest disease signature (most "healthy-like")
    disease_markers = ["PNPLA3", "TM6SF2", "MBOAT7", "HSD17B13", "COL1A1", "ACTA2"]
    present = [g for g in disease_markers if g in adata.var_names]

    if present:
        sc.tl.score_genes(adata, gene_list=present, score_name="disease_score")
    else:
        # Fallback: use periportal score (healthy hepatocytes)
        if "periportal_score" in adata.obs.columns:
            adata.obs["disease_score"] = -adata.obs["periportal_score"]
        else:
            adata.obs["disease_score"] = 0

    # Root = spot with lowest disease score
    root_idx = adata.obs["disease_score"].idxmin()
    adata.uns["iroot"] = adata.obs.index.get_loc(root_idx)

    # Diffusion pseudotime
    sc.tl.dpt(adata)
    print(f"  DPT range: [{adata.obs['dpt_pseudotime'].min():.3f}, "
          f"{adata.obs['dpt_pseudotime'].max():.3f}]")
    return adata


def identify_transition_zones(adata, percentile=90):
    """Identify spatial locations with rapid pseudotime changes."""
    sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6)
    conn = adata.obsp["spatial_connectivities"]
    dpt = adata.obs["dpt_pseudotime"].values

    gradients = np.zeros(adata.n_obs)
    for i in range(adata.n_obs):
        neighbors = conn[i].nonzero()[1]
        if len(neighbors) > 0:
            gradients[i] = np.std(dpt[neighbors])

    adata.obs["dpt_gradient"] = gradients
    threshold = np.percentile(gradients[gradients > 0], percentile)
    adata.obs["transition_zone"] = adata.obs["dpt_gradient"] > threshold

    n_tz = adata.obs["transition_zone"].sum()
    print(f"  Transition zones: {n_tz} spots ({n_tz/adata.n_obs*100:.1f}%)")
    return adata


def gene_dynamics_along_trajectory(adata, genes, n_bins=20):
    """Compute gene expression dynamics along pseudotime."""
    dpt = adata.obs["dpt_pseudotime"].values
    bins = pd.qcut(dpt, q=n_bins, labels=False, duplicates="drop")

    dynamics = []
    for gene in genes:
        if gene not in adata.var_names:
            continue
        expr = np.asarray(adata[:, gene].X.todense()).flatten()
        for b in range(int(bins.max()) + 1):
            mask = bins == b
            if mask.sum() > 0:
                dynamics.append({
                    "gene": gene,
                    "bin": b,
                    "mean_dpt": dpt[mask].mean(),
                    "mean_expr": expr[mask].mean(),
                    "std_expr": expr[mask].std(),
                    "n_spots": mask.sum(),
                })
    return pd.DataFrame(dynamics)


def main():
    print_header("05f: Spatial Disease Trajectory")

    config = load_config()
    output_dir = RESULTS_DIR / "trajectory"
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = load_deconvolved_adata()
    print(f"  Loaded: {adata.n_obs} spots")

    # Filter to MASLD condition only (trajectory doesn't apply to healthy baseline)
    conditions = adata.obs["condition"].unique().tolist()
    masld_conditions = [c for c in conditions if "healthy" not in c.lower() and c != "Healthy"]

    if not masld_conditions:
        print("  WARNING: No MASLD samples found. Trajectory analysis requires disease data.")
        print("  Skipping trajectory analysis.")
        return

    masld_cond = masld_conditions[0]
    adata_masld = adata[adata.obs["condition"] == masld_cond].copy()
    print(f"  MASLD condition: {masld_cond} ({adata_masld.n_obs} spots)")

    # Compute pseudotime
    print("\n  Computing diffusion pseudotime...")
    adata_masld = compute_spatial_pseudotime(adata_masld)

    # Identify transition zones
    print("\n  Identifying transition zones...")
    adata_masld = identify_transition_zones(adata_masld)

    # Save pseudotime results
    traj_df = adata_masld.obs[[
        "sample_id", "dpt_pseudotime", "dpt_gradient", "transition_zone",
        "disease_score",
    ]].copy()
    if "zonation_score" in adata_masld.obs.columns:
        traj_df["zonation_score"] = adata_masld.obs["zonation_score"]
    save_csv(traj_df, "spatial_pseudotime.csv", subdir="trajectory")

    # Gene dynamics along trajectory for key genes
    print("\n  Computing gene dynamics along trajectory...")
    key_genes = [
        # Lipid metabolism
        "FASN", "SCD", "ACACA", "DGAT2", "PNPLA3",
        # Fibrosis
        "COL1A1", "COL1A2", "ACTA2", "TGFB1", "PDGFRB",
        # Inflammation
        "TNF", "IL1B", "CCL2", "CD68", "TREM2",
        # Zonation
        "CYP2E1", "GLUL", "HAL", "ASS1",
    ]
    dynamics = gene_dynamics_along_trajectory(adata_masld, key_genes)
    save_csv(dynamics, "gene_dynamics_along_trajectory.csv", subdir="trajectory")
    print(f"  Gene dynamics: {len(dynamics)} data points for "
          f"{dynamics['gene'].nunique()} genes")

    # Correlation of pseudotime with zonation — donor-blocked (F063).
    # A pooled spot-level spearmanr p would be pseudoreplicated (n=spots, ~2
    # donors). Report the mean per-donor rho; the donor-level p is descriptive
    # only at n~2 donors (spearman_by_donor returns NaN p below 3 donors).
    if "zonation_score" in adata_masld.obs.columns and DONOR_COL in adata_masld.obs.columns:
        sb = spearman_by_donor(
            adata_masld.obs["dpt_pseudotime"].values,
            adata_masld.obs["zonation_score"].values,
            adata_masld.obs[DONOR_COL].astype(str).values,
        )
        print(f"\n  Pseudotime vs zonation (donor-blocked): "
              f"mean_rho={sb['mean_rho']:.3f} over n_donors={sb['n_donors']} "
              f"(donor-level p={sb.get('pval', float('nan'))}; descriptive only)")

    # Save annotated AnnData
    save_checkpoint(adata_masld, "spatial_with_trajectory.h5ad", subdir="trajectory")

    print_header("05f: Complete")


if __name__ == "__main__":
    main()

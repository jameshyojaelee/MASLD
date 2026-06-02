#!/usr/bin/env python3
"""
05a_neighborhood_enrichment.py — Spatial cell-type neighborhood analysis.

Computes cell-type co-localization patterns and compares between conditions.
Uses squidpy for neighborhood enrichment via permutation testing.

Donor-aware (2026-06-01 rigor pass, F036): the squidpy permutation test is run
PER DONOR (sample_id) and the per-donor z-score matrices are averaged within a
condition, instead of pooling ~10^3 spatially-autocorrelated spots across the
~2 donors/condition as if they were independent (GSE192741: Healthy=JBO018/022,
Steatotic=JBO014/015/019). The per-condition z-score matrix is therefore a
donor-mean of per-donor z-scores; with ~2 donors/condition it is descriptive,
not a powered test. The n_donors per condition is logged.

SLURM: --partition=cpu --cpus=8 --mem=64G --time=4:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

np.random.seed(42)

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, load_deconvolved_adata,
    save_csv, print_header, print_step,
)

# Donor (biological replicate) column for GSE192741 — DO NOT pool across these.
DONOR_COL = "sample_id"
# Minimum spots per donor slice to attempt a nhood-enrichment permutation test.
MIN_SPOTS_PER_DONOR = 100


def _nhood_one_slice(adata_slice, n_neighs, n_perms):
    """nhood_enrichment z-score matrix for a single donor slice (or None)."""
    adata_slice = adata_slice.copy()
    # Restrict categories to those present in THIS slice so the z-score matrix
    # never carries zero-spot cell types across donors.
    adata_slice.obs["cell_type_dominant"] = (
        adata_slice.obs["cell_type_dominant"].astype("category")
        .cat.remove_unused_categories()
    )
    sq.gr.spatial_neighbors(adata_slice, coord_type="generic", n_neighs=n_neighs)
    sq.gr.nhood_enrichment(
        adata_slice, cluster_key="cell_type_dominant",
        n_perms=n_perms, seed=42,
    )
    zscore_matrix = adata_slice.uns["cell_type_dominant_nhood_enrichment"]["zscore"]
    ct_labels = adata_slice.obs["cell_type_dominant"].cat.categories
    return pd.DataFrame(zscore_matrix, index=ct_labels, columns=ct_labels)


def run_neighborhood_enrichment(adata, condition, n_neighs=6, n_perms=1000):
    """Per-donor neighborhood enrichment for one condition, averaged (F036).

    Each donor slice is tested independently (no cross-donor / cross-slice
    pooling), then the per-donor z-score matrices are averaged cell-type-pair
    wise (aligned on the union of cell types, NaN-skipping). Returns the
    donor-mean z-score matrix, or None if no donor has enough spots.
    """
    adata_cond = adata[adata.obs["condition"] == condition]
    donors = adata_cond.obs[DONOR_COL].astype(str)
    per_donor = []
    for donor in sorted(donors.unique()):
        sl = adata_cond[donors.values == donor]
        if sl.n_obs < MIN_SPOTS_PER_DONOR:
            print(f"    {donor}: {sl.n_obs} spots < {MIN_SPOTS_PER_DONOR}, skipping")
            continue
        try:
            per_donor.append(_nhood_one_slice(sl, n_neighs, n_perms))
        except Exception as e:
            print(f"    WARNING: nhood_enrichment failed for {donor}: {e}")
    if not per_donor:
        print(f"    WARNING: no donor with >= {MIN_SPOTS_PER_DONOR} spots for {condition}")
        return None

    # Average per-donor z-score matrices over the union of cell types.
    all_cts = sorted(set().union(*[m.index for m in per_donor]))
    stacked = [m.reindex(index=all_cts, columns=all_cts) for m in per_donor]
    df = pd.concat(stacked).groupby(level=0).mean().reindex(index=all_cts, columns=all_cts)
    print(f"    Averaged nhood z-scores over {len(per_donor)} donor(s)")
    return df


def differential_neighborhoods(zscore_healthy, zscore_masld):
    """Identify cell-type pairs with changed spatial co-localization."""
    shared_cts = zscore_healthy.index.intersection(zscore_masld.index)
    z_h = zscore_healthy.loc[shared_cts, shared_cts]
    z_m = zscore_masld.loc[shared_cts, shared_cts]
    delta = z_m - z_h

    pairs = []
    cts = list(shared_cts)
    for i, ct1 in enumerate(cts):
        for j, ct2 in enumerate(cts):
            if i < j:
                pairs.append({
                    "cell_type_1": ct1, "cell_type_2": ct2,
                    "z_healthy": z_h.loc[ct1, ct2],
                    "z_masld": z_m.loc[ct1, ct2],
                    "delta_z": delta.loc[ct1, ct2],
                    "abs_delta": abs(delta.loc[ct1, ct2]),
                })
    return pd.DataFrame(pairs).sort_values("abs_delta", ascending=False)


def main():
    print_header("05a: Spatial Neighborhood Enrichment")

    config = load_config()
    comm_config = config["communication"]
    output_dir = RESULTS_DIR / "communication"
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = load_deconvolved_adata()
    print(f"  Loaded: {adata.n_obs} spots")

    conditions = adata.obs["condition"].unique().tolist()
    print(f"  Conditions: {conditions}")

    zscore_matrices = {}
    for condition in conditions:
        cond_donors = adata.obs.loc[adata.obs["condition"] == condition, DONOR_COL].nunique()
        print(f"\n  Computing neighborhood enrichment: {condition} (n={cond_donors} donors)...")
        zscore = run_neighborhood_enrichment(
            adata, condition,
            n_neighs=comm_config["n_neighbors"],
            n_perms=comm_config["n_perms"],
        )
        if zscore is not None:
            zscore_matrices[condition] = zscore
            save_csv(zscore, f"nhood_enrichment_{condition}_zscore.csv",
                     subdir="communication")
            # Top co-enriched pairs
            pairs = []
            for i, ct1 in enumerate(zscore.index):
                for j, ct2 in enumerate(zscore.columns):
                    if i < j:
                        pairs.append((ct1, ct2, zscore.loc[ct1, ct2]))
            pairs.sort(key=lambda x: x[2], reverse=True)
            print(f"    Top 5 co-enriched pairs:")
            for ct1, ct2, z in pairs[:5]:
                print(f"      {ct1} ↔ {ct2}: z={z:.2f}")

    # Differential analysis
    if len(zscore_matrices) >= 2:
        cond_list = list(zscore_matrices.keys())
        # Assume first is healthy, second is disease
        cond_h = [c for c in cond_list if "healthy" in c.lower() or c == "Healthy"]
        cond_m = [c for c in cond_list if c not in cond_h]
        if cond_h and cond_m:
            print(f"\n  Differential neighborhoods: {cond_h[0]} vs {cond_m[0]}...")
            diff = differential_neighborhoods(
                zscore_matrices[cond_h[0]], zscore_matrices[cond_m[0]]
            )
            save_csv(diff, "differential_neighborhoods.csv", subdir="communication")
            print(f"  Top 10 changed interactions:")
            for _, row in diff.head(10).iterrows():
                direction = "↑" if row["delta_z"] > 0 else "↓"
                print(f"    {row['cell_type_1']} ↔ {row['cell_type_2']}: "
                      f"Δz={row['delta_z']:+.2f} {direction}")

    print_header("05a: Complete")


if __name__ == "__main__":
    main()

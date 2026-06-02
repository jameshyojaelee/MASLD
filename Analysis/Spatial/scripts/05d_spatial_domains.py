#!/usr/bin/env python3
"""
05d_spatial_domains.py — Unsupervised spatial domain identification.

Groups spots into tissue microenvironments using combined expression +
spatial proximity clustering. Characterizes domains by cell type composition,
zonation, and marker genes.

SLURM: --partition=gpu --gres=gpu:1 --cpus=8 --mem=64G --time=6:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from sklearn.preprocessing import normalize

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, load_deconvolved_adata,
    save_checkpoint, save_csv, print_header,
)


def identify_spatial_domains(adata, config):
    """Identify spatial domains via spatially-smoothed Leiden clustering."""
    alpha = config["spatial_weight"]
    n_expr_neighbors = config["n_expression_neighbors"]
    n_spatial_neighbors = config["n_spatial_neighbors"]
    target = config["target_n_domains"]

    # PCA on HVGs
    if "highly_variable" in adata.var.columns:
        sc.pp.pca(adata, n_comps=30, use_highly_variable=True)
    else:
        sc.pp.pca(adata, n_comps=30)

    # Expression-based neighbors
    sc.pp.neighbors(adata, n_neighbors=n_expr_neighbors, n_pcs=30,
                    key_added="expression_neighbors")

    # Spatial neighbors
    sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=n_spatial_neighbors,
                            key_added="spatial_neighbors")

    # Combine connectivity matrices
    expr_conn = adata.obsp["expression_neighbors_connectivities"]
    spat_conn = adata.obsp["spatial_neighbors_connectivities"]

    expr_norm = normalize(expr_conn, norm="l1", axis=1)
    spat_norm = normalize(spat_conn, norm="l1", axis=1)
    combined = (1 - alpha) * expr_norm + alpha * spat_norm

    adata.obsp["connectivities"] = combined
    adata.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "params": {"method": "combined", "alpha": alpha},
    }

    # Find resolution giving ~target domains
    # F052: pin random_state on Leiden so the domain count is reproducible
    # (flavor='igraph' still has its own RNG; random_state seeds it). The chosen
    # resolution remains a grid/proximity-to-target pick — record best_res in
    # adata.uns so the actual resolution used is persisted for methods/reporting.
    best_res = 0.5
    best_diff = 999
    res_sweep = []
    for resolution in config["resolution_range"]:
        sc.tl.leiden(adata, resolution=resolution, key_added="spatial_domain",
                    flavor="igraph", n_iterations=2, directed=False, random_state=42)
        n_clusters = adata.obs["spatial_domain"].nunique()
        res_sweep.append({"resolution": resolution, "n_domains": int(n_clusters)})
        diff = abs(n_clusters - target)
        if diff < best_diff:
            best_diff = diff
            best_res = resolution
        if diff <= 1:
            break

    sc.tl.leiden(adata, resolution=best_res, key_added="spatial_domain",
                flavor="igraph", n_iterations=2, directed=False, random_state=42)
    n_domains = adata.obs["spatial_domain"].nunique()
    adata.uns["spatial_domain_resolution"] = best_res
    # store as dict-of-arrays (h5ad-serializable); a list-of-dicts fails the .uns write
    adata.uns["spatial_domain_res_sweep"] = {
        "resolution": [float(r["resolution"]) for r in res_sweep],
        "n_domains": [int(r["n_domains"]) for r in res_sweep],
    }
    print(f"  Identified {n_domains} domains (resolution={best_res}, target={target})")
    print(f"  Resolution sweep (count vs resolution): {res_sweep}")
    return adata


def characterize_domains(adata):
    """Characterize each domain by cell types, zonation, markers."""
    # Cell type composition per domain
    if "cell_type_dominant" in adata.obs.columns:
        comp = pd.crosstab(
            adata.obs["spatial_domain"],
            adata.obs["cell_type_dominant"],
            normalize="index",
        )
    else:
        comp = pd.DataFrame()

    # Marker genes per domain
    sc.tl.rank_genes_groups(adata, groupby="spatial_domain", method="wilcoxon", n_genes=50)
    markers = sc.get.rank_genes_groups_df(adata, group=None)

    # Zonation per domain
    zon = pd.DataFrame()
    if "zonation_score" in adata.obs.columns:
        zon = adata.obs.groupby("spatial_domain")["zonation_score"].agg(["mean", "std", "count"])

    # Domain sizes
    sizes = adata.obs["spatial_domain"].value_counts().rename("n_spots")

    return comp, markers, zon, sizes


def compare_domains_between_conditions(adata):
    """Compare domain proportions between conditions."""
    if "condition" not in adata.obs.columns:
        return pd.DataFrame()

    # Per-sample domain proportions
    props = pd.crosstab(
        [adata.obs["sample_id"], adata.obs["condition"]],
        adata.obs["spatial_domain"],
        normalize="index",
    )
    props = props.reset_index()
    return props


def main():
    print_header("05d: Spatial Domain Identification")

    config = load_config()
    dom_config = config["domains"]
    output_dir = RESULTS_DIR / "domains"
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = load_deconvolved_adata()
    print(f"  Loaded: {adata.n_obs} spots")

    # F055: drop mitochondrial / ribosomal / hemoglobin genes BEFORE HVG/PCA/
    # clustering and marker calling. Otherwise MT- genes (MT-CO1/2, MT-ND1..4,
    # MT-CYB) dominate the HVG set and rank_genes_groups markers, so domains are
    # partly defined by technical mito content rather than biology. Standard
    # Visium domain workflow.
    import re
    techn_mask = adata.var_names.str.match(r"^(MT-|MT_|RPS\d|RPL\d|HB[AB])", case=False)
    n_drop = int(techn_mask.sum())
    if n_drop > 0:
        print(f"  Dropping {n_drop} MT-/ribosomal/hemoglobin genes before clustering")
        adata = adata[:, ~techn_mask].copy()
        # Invalidate any stale HVG flag so it is recomputed on the filtered set.
        if "highly_variable" in adata.var.columns:
            del adata.var["highly_variable"]

    # Ensure log-normalized expression with HVGs
    if "counts" in adata.layers:
        # HVG selection on raw counts (seurat_v3 requires raw counts)
        adata.X = adata.layers["counts"].copy()
        sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3")
        # Then normalize for PCA/clustering
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    elif "highly_variable" not in adata.var.columns:
        sc.pp.highly_variable_genes(adata, n_top_genes=3000)

    # Identify domains
    print("\n  Identifying spatial domains...")
    adata = identify_spatial_domains(adata, dom_config)

    # Characterize
    print("\n  Characterizing domains...")
    comp, markers, zon, sizes = characterize_domains(adata)

    # Save results
    save_csv(comp, "domain_composition.csv", subdir="domains")
    save_csv(markers, "domain_markers.csv", subdir="domains")
    if len(zon) > 0:
        save_csv(zon, "domain_zonation.csv", subdir="domains")

    domain_summary = adata.obs[["sample_id", "dataset", "condition", "spatial_domain"]].copy()
    save_csv(domain_summary, "spatial_domains.csv", subdir="domains")

    # Print summary
    print(f"\n  Domain summary:")
    for domain in sorted(adata.obs["spatial_domain"].unique()):
        n = sizes.get(domain, 0)
        pct = n / adata.n_obs * 100
        print(f"    Domain {domain}: {n} spots ({pct:.1f}%)")
        if len(comp) > 0 and domain in comp.index:
            top_ct = comp.loc[domain].nlargest(3)
            ct_str = ", ".join(f"{ct}={v:.0%}" for ct, v in top_ct.items())
            print(f"      Top cell types: {ct_str}")
        if len(zon) > 0 and domain in zon.index:
            print(f"      Zonation: {zon.loc[domain, 'mean']:.3f} ± {zon.loc[domain, 'std']:.3f}")

    # Compare between conditions
    print("\n  Comparing domains between conditions...")
    domain_props = compare_domains_between_conditions(adata)
    if len(domain_props) > 0:
        save_csv(domain_props, "domain_proportions_by_sample.csv", subdir="domains")

    # Save AnnData with domains
    save_checkpoint(adata, "spatial_with_domains.h5ad", subdir="domains")

    print_header("05d: Complete")


if __name__ == "__main__":
    main()

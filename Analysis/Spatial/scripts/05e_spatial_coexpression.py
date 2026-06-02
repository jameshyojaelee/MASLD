#!/usr/bin/env python3
"""
05e_spatial_coexpression.py — Spatially-resolved gene co-expression modules.

Uses Hotspot to identify gene modules co-expressed in spatially coherent
patterns. Annotates modules with bulk evidence (dream DEGs, Conserved).

SLURM: --partition=cpu --cpus=16 --mem=128G --time=8:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, load_spatial_adata, load_dream_degs,
    load_conserved, save_csv, print_header,
)

# Donor (slide) column. GSE192741 = JBO014/015/018/019/022 (n=5; 2 Healthy,
# 3 Steatotic). Used to prevent cross-slide spatial edges (F056).
DONOR_COL = "sample_id"


def offset_spatial_per_slice(adata, donor_col=DONOR_COL, gap_frac=2.0):
    """Offset each slide's spatial coordinates so they occupy disjoint blocks.

    F056: within a condition the spots of 2-3 donor slides are pooled into ONE
    generic-kNN spatial graph. Visium slide coordinates share an overlapping
    frame, so the KNN can connect spots from physically distinct slides,
    corrupting Hotspot autocorrelation with cross-donor edges. Translating each
    slide along x by a multiple of the global span guarantees inter-slide
    distances exceed any intra-slide distance, so no spatial neighbor edge can
    cross a slide boundary. Module discovery is still donor-POOLED (modules are
    not a per-donor consensus); this only removes the spurious cross-slide
    edges. Returns the modified adata (operates in place on obsm['spatial']).
    """
    if donor_col not in adata.obs.columns:
        print(f"    WARNING: '{donor_col}' missing; cannot offset slices "
              f"(cross-slide edges possible)")
        return adata
    coords = np.asarray(adata.obsm["spatial"], dtype=float).copy()
    span_x = float(coords[:, 0].max() - coords[:, 0].min()) or 1.0
    donors = adata.obs[donor_col].astype(str).values
    for k, d in enumerate(pd.unique(donors)):
        mask = donors == d
        coords[mask, 0] = coords[mask, 0] + k * gap_frac * span_x
    adata.obsm["spatial"] = coords
    return adata


def run_hotspot(adata, config):
    """Identify spatially co-expressed gene modules using Hotspot."""
    import hotspot

    # Pre-compute spatial neighbors for Hotspot
    import squidpy as sq
    sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=config["n_neighbors"])

    hs = hotspot.Hotspot(
        adata,
        model=config["hotspot_model"],
        distances_obsp_key="spatial_distances",
    )

    # Build KNN graph from the spatial distances (required before autocorrelations)
    hs.create_knn_graph(
        weighted_graph=False,
        n_neighbors=config["n_neighbors"],
    )

    # Step 1: Per-gene spatial autocorrelation
    print("    Computing autocorrelations...")
    hs_results = hs.compute_autocorrelations(jobs=config["n_jobs"])
    sig_genes = hs_results[hs_results["FDR"] < config["fdr_threshold"]].index.tolist()
    print(f"    Spatially autocorrelated genes: {len(sig_genes)}")

    if len(sig_genes) < config["min_gene_threshold"]:
        print(f"    WARNING: Too few significant genes ({len(sig_genes)} < {config['min_gene_threshold']})")
        return hs_results, pd.DataFrame(), pd.DataFrame()

    # Step 2: Pairwise local correlations
    print("    Computing local correlations...")
    lcz = hs.compute_local_correlations(sig_genes, jobs=config["n_jobs"])

    # Step 3: Module identification
    print("    Identifying modules...")
    modules = hs.create_modules(
        min_gene_threshold=config["min_gene_threshold"],
        core_only=True,
        fdr_threshold=config["fdr_threshold"],
    )

    # Step 4: Module scores per spot
    print("    Computing module scores...")
    module_scores = hs.calculate_module_scores()

    return hs_results, modules, module_scores


def annotate_modules(modules, dream_degs, conserved):
    """Annotate spatial modules with bulk evidence."""
    gene_col = "symbol" if "symbol" in dream_degs.columns else dream_degs.columns[0]
    deg_set = set(dream_degs[gene_col])
    core_set = set(conserved)

    annotations = []
    for mod_id in modules["Module"].unique():
        if mod_id == -1:
            continue
        genes = modules[modules["Module"] == mod_id].index.tolist()
        n_degs = len(set(genes) & deg_set)
        n_core = len(set(genes) & core_set)
        annotations.append({
            "module": mod_id,
            "n_genes": len(genes),
            "n_dream_degs": n_degs,
            "pct_degs": n_degs / max(len(genes), 1) * 100,
            "n_conserved": n_core,
            "pct_conserved": n_core / max(len(genes), 1) * 100,
            "top_genes": ", ".join(genes[:15]),
        })
    return pd.DataFrame(annotations)


def main():
    print_header("05e: Spatial Co-expression Module Detection")

    config = load_config()
    coex_config = config["coexpression"]
    output_dir = RESULTS_DIR / "coexpression"
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = load_spatial_adata()
    print(f"  Loaded: {adata.n_obs} spots, {adata.n_vars} genes")

    # Ensure counts for Hotspot
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()

    # Run per condition
    conditions = adata.obs["condition"].unique().tolist()

    for condition in conditions:
        print(f"\n  Processing: {condition}")
        adata_sub = adata[adata.obs["condition"] == condition].copy()
        if adata_sub.n_obs < 200:
            print(f"    Skipping (too few spots: {adata_sub.n_obs})")
            continue

        # F056: document donor structure and offset slides so the pooled spatial
        # KNN graph cannot form cross-slide edges. Modules remain donor-pooled
        # (not a per-donor consensus) — caveat for the figure/methods.
        n_donors = (adata_sub.obs[DONOR_COL].nunique()
                    if DONOR_COL in adata_sub.obs.columns else None)
        print(f"    n_donors (slides) in {condition}: {n_donors} "
              f"(modules are donor-POOLED, not consensus)")
        adata_sub = offset_spatial_per_slice(adata_sub)

        try:
            hs_results, modules, module_scores = run_hotspot(adata_sub, coex_config)

            save_csv(hs_results, f"spatial_autocorr_{condition}.csv", subdir="coexpression")

            if len(modules) > 0:
                save_csv(modules, f"modules_{condition}.csv", subdir="coexpression")
                n_modules = modules["Module"].nunique() - (1 if -1 in modules["Module"].values else 0)
                print(f"    Modules identified: {n_modules}")

                for mod_id in sorted(modules["Module"].unique()):
                    if mod_id == -1:
                        continue
                    n = (modules["Module"] == mod_id).sum()
                    top = modules[modules["Module"] == mod_id].index[:5].tolist()
                    print(f"      Module {mod_id}: {n} genes — {', '.join(top)}")

            if isinstance(module_scores, pd.DataFrame) and len(module_scores) > 0:
                save_csv(module_scores, f"module_scores_{condition}.csv", subdir="coexpression")

        except Exception as e:
            print(f"    ERROR: Hotspot failed for {condition}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Annotate modules with bulk evidence
    print("\n  Annotating modules with bulk evidence...")
    try:
        dream_degs = load_dream_degs(padj_thresh=0.1, lfc_thresh=0.8)
        conserved = load_conserved()

        for condition in conditions:
            mod_path = RESULTS_DIR / "coexpression" / f"modules_{condition}.csv"
            if mod_path.exists():
                modules = pd.read_csv(mod_path, index_col=0)
                annot = annotate_modules(modules, dream_degs, conserved)
                save_csv(annot, f"module_annotation_{condition}.csv", subdir="coexpression")
                print(f"\n  {condition} module annotations:")
                for _, row in annot.iterrows():
                    print(f"    Module {row['module']}: {row['n_genes']} genes, "
                          f"{row['pct_degs']:.0f}% DEGs, {row['n_conserved']} Conserved")
    except Exception as e:
        print(f"    WARNING: Annotation failed: {e}")

    print_header("05e: Complete")


if __name__ == "__main__":
    main()

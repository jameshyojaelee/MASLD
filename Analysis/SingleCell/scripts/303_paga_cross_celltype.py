#!/usr/bin/env python3
"""
303: PAGA Cross-Cell-Type Connectivity Analysis.

Runs partition-based graph abstraction (PAGA) on the full atlas to map
connectivity between cell sub-states across conditions, then compares
Healthy vs. NASH connectivity to identify disease-driven inter-cell-type
communication axes.

Inputs:
    - Integrated atlas: Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad

Outputs (to results_gpu_v2/pseudotime/):
    - paga_connectivity_all.csv
    - paga_connectivity_healthy.csv
    - paga_connectivity_nash.csv
    - paga_differential_connectivity.csv
    - paga_global_pseudotime.csv
    - paga_cluster_annotations.csv

Usage:
    sbatch run_pseudotime_pipeline.sh
"""

import os
import sys
import warnings
import logging
import gc

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

warnings.filterwarnings("ignore")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
SC_DIR = os.path.join(BASE, "Analysis/SingleCell")
RESULTS = os.path.join(SC_DIR, "results_gpu_v2")
H5AD = os.path.join(RESULTS, "integrated_atlas.h5ad")
PT_DIR = os.path.join(RESULTS, "pseudotime")

# Target cell types for annotation
CORE_TYPES = [
    "Hepatocytes", "Macrophages", "Fibroblasts",
    "Endothelial cells", "Cholangiocytes",
]

MAX_CELLS = 100_000  # Subsample for PAGA tractability

# ---------------------------------------------------------------------------
# GPU init
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.join(SC_DIR, "scripts"))
try:
    from gpu_utils import init_gpu, get_processor
    USE_GPU = init_gpu()
except Exception:
    USE_GPU = False
    log.info("No GPU — running CPU-only")

import scanpy as sc
import anndata as ad

if USE_GPU:
    pp, tl = get_processor(True)
    from gpu_utils import to_gpu, from_gpu
else:
    pp, tl = sc.pp, sc.tl


# ---------------------------------------------------------------------------
# Fix h5ad categoricals
# ---------------------------------------------------------------------------
def fix_h5ad_categoricals(h5ad_path):
    """Patch h5ad: add missing 'ordered' attr to categoricals."""
    with h5py.File(h5ad_path, "a") as f:
        for col in f["obs"].keys():
            if col == "_index":
                continue
            attrs = f["obs"][col].attrs
            enc = attrs.get("encoding-type", b"")
            if isinstance(enc, bytes):
                enc = enc.decode()
            if enc == "categorical" and "ordered" not in attrs:
                attrs["ordered"] = False


# ---------------------------------------------------------------------------
# Extract PAGA connectivity as DataFrame
# ---------------------------------------------------------------------------
def paga_to_dataframe(adata, group_key):
    """Extract PAGA connectivity matrix as edge list DataFrame."""
    conn = adata.uns["paga"]["connectivities"].toarray()
    labels = adata.obs[group_key].cat.categories.tolist()

    edges = []
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            if conn[i, j] > 0:
                edges.append({
                    "source": labels[i],
                    "target": labels[j],
                    "weight": conn[i, j],
                })
    return pd.DataFrame(edges)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log.info("=" * 60)
    log.info("303: PAGA Cross-Cell-Type Connectivity")
    log.info("=" * 60)

    # Load atlas (with categorical fix — does NOT modify original file)
    log.info("Loading atlas from %s ...", H5AD)
    try:
        adata = sc.read_h5ad(H5AD)
    except KeyError:
        # Fix categoricals on a temp copy
        import shutil
        tmp = os.path.join(PT_DIR, ".atlas_patched.h5ad")
        if not os.path.exists(tmp):
            shutil.copy2(H5AD, tmp)
            fix_h5ad_categoricals(tmp)
        adata = sc.read_h5ad(tmp)
    log.info("Atlas shape: %s", adata.shape)

    # Focus on core cell types
    core_mask = adata.obs["cell_type"].isin(CORE_TYPES)
    log.info("Core cell types: %d / %d cells", core_mask.sum(), len(adata))
    adata = adata[core_mask].copy()

    # Subsample if needed
    if len(adata) > MAX_CELLS:
        log.info("Subsampling from %d to %d cells...", len(adata), MAX_CELLS)
        np.random.seed(42)
        # Stratified subsampling by cell type
        keep_idx = []
        per_type = MAX_CELLS // len(CORE_TYPES)
        for ct in CORE_TYPES:
            ct_idx = np.where(adata.obs["cell_type"] == ct)[0]
            n_take = min(len(ct_idx), per_type)
            keep_idx.extend(np.random.choice(ct_idx, n_take, replace=False))
        adata = adata[sorted(keep_idx)].copy()
        log.info("After subsampling: %d cells", len(adata))

    # Data is already log-normalized — do NOT re-normalize
    log.info("Embedding (data already log-normalized)...")

    if USE_GPU:
        to_gpu(adata)

    # HVG selection — seurat flavor works on log-normalized data
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat")
    pp.scale(adata, max_value=10)
    tl.pca(adata, n_comps=50)
    pp.neighbors(adata, n_neighbors=30, n_pcs=50)
    tl.umap(adata)

    # Fine-resolution clustering for sub-state identification
    log.info("Fine-resolution clustering (res=1.5)...")
    tl.leiden(adata, resolution=1.5, key_added="leiden_fine")

    if USE_GPU:
        from_gpu(adata)

    n_clusters = adata.obs["leiden_fine"].nunique()
    log.info("Found %d fine clusters", n_clusters)

    # Annotate clusters: cell_type + cluster ID
    adata.obs["celltype_cluster"] = (
        adata.obs["cell_type"].astype(str) + "_" +
        adata.obs["leiden_fine"].astype(str)
    )

    # Dominant cell type per cluster (for labeling)
    cluster_annot = (
        adata.obs.groupby("leiden_fine")["cell_type"]
        .agg(lambda x: x.mode().iloc[0])
        .reset_index()
        .rename(columns={"cell_type": "dominant_celltype"})
    )
    cluster_annot["n_cells"] = (
        adata.obs.groupby("leiden_fine").size().values
    )

    # Condition composition per cluster
    cond_comp = pd.crosstab(
        adata.obs["leiden_fine"], adata.obs["condition"], normalize="index"
    )
    # Reset categorical index to plain strings for merge compatibility
    cond_comp.index = cond_comp.index.astype(str)
    cond_comp.columns = cond_comp.columns.astype(str)
    cluster_annot["leiden_fine"] = cluster_annot["leiden_fine"].astype(str)
    cluster_annot = cluster_annot.merge(cond_comp, left_on="leiden_fine",
                                         right_index=True, how="left")
    cluster_annot.to_csv(os.path.join(PT_DIR, "paga_cluster_annotations.csv"),
                          index=False)

    # =====================================================================
    # PAGA on full (subsampled) atlas
    # =====================================================================
    log.info("\n--- PAGA: All conditions ---")
    sc.tl.paga(adata, groups="leiden_fine")

    paga_all = paga_to_dataframe(adata, "leiden_fine")
    paga_all.to_csv(os.path.join(PT_DIR, "paga_connectivity_all.csv"), index=False)
    log.info("PAGA all: %d edges", len(paga_all))

    # =====================================================================
    # Condition-stratified PAGA: Healthy
    # =====================================================================
    log.info("\n--- PAGA: Healthy only ---")
    healthy_mask = adata.obs["condition"] == "Healthy"
    if healthy_mask.sum() > 500:
        adata_h = adata[healthy_mask].copy()
        # Re-compute neighbors on this subset
        sc.pp.neighbors(adata_h, n_neighbors=min(30, len(adata_h) // 10), n_pcs=50)
        sc.tl.paga(adata_h, groups="leiden_fine")
        paga_healthy = paga_to_dataframe(adata_h, "leiden_fine")
        paga_healthy.to_csv(os.path.join(PT_DIR, "paga_connectivity_healthy.csv"),
                             index=False)
        log.info("PAGA healthy: %d edges", len(paga_healthy))
        del adata_h
    else:
        paga_healthy = pd.DataFrame()
        log.warning("Too few Healthy cells (%d) for PAGA", healthy_mask.sum())

    # =====================================================================
    # Condition-stratified PAGA: NASH
    # =====================================================================
    log.info("\n--- PAGA: NASH only ---")
    nash_mask = adata.obs["condition"] == "NASH"
    if nash_mask.sum() > 500:
        adata_n = adata[nash_mask].copy()
        sc.pp.neighbors(adata_n, n_neighbors=min(30, len(adata_n) // 10), n_pcs=50)
        sc.tl.paga(adata_n, groups="leiden_fine")
        paga_nash = paga_to_dataframe(adata_n, "leiden_fine")
        paga_nash.to_csv(os.path.join(PT_DIR, "paga_connectivity_nash.csv"),
                          index=False)
        log.info("PAGA NASH: %d edges", len(paga_nash))
        del adata_n
    else:
        paga_nash = pd.DataFrame()
        log.warning("Too few NASH cells (%d) for PAGA", nash_mask.sum())

    # =====================================================================
    # Differential connectivity
    # =====================================================================
    if len(paga_healthy) > 0 and len(paga_nash) > 0:
        log.info("\n--- Differential PAGA connectivity ---")

        # Merge on source-target pairs
        paga_healthy_indexed = paga_healthy.set_index(["source", "target"])
        paga_nash_indexed = paga_nash.set_index(["source", "target"])

        all_edges = paga_healthy_indexed.index.union(paga_nash_indexed.index)
        diff_rows = []
        for edge in all_edges:
            w_h = paga_healthy_indexed.loc[edge, "weight"] if edge in paga_healthy_indexed.index else 0
            w_n = paga_nash_indexed.loc[edge, "weight"] if edge in paga_nash_indexed.index else 0

            # Get cell types for this edge
            src_ct = cluster_annot.loc[
                cluster_annot["leiden_fine"] == str(edge[0]), "dominant_celltype"
            ]
            tgt_ct = cluster_annot.loc[
                cluster_annot["leiden_fine"] == str(edge[1]), "dominant_celltype"
            ]
            src_ct = src_ct.iloc[0] if len(src_ct) > 0 else "Unknown"
            tgt_ct = tgt_ct.iloc[0] if len(tgt_ct) > 0 else "Unknown"

            diff_rows.append({
                "source": edge[0],
                "target": edge[1],
                "source_celltype": src_ct,
                "target_celltype": tgt_ct,
                "weight_healthy": w_h,
                "weight_nash": w_n,
                "diff": w_n - w_h,
                "cross_celltype": src_ct != tgt_ct,
            })

        diff_df = pd.DataFrame(diff_rows).sort_values("diff", ascending=False)
        diff_df.to_csv(os.path.join(PT_DIR, "paga_differential_connectivity.csv"),
                        index=False)

        # Report cross-cell-type edges gained in NASH
        cross = diff_df[diff_df["cross_celltype"]]
        gained = cross[cross["diff"] > 0.1].sort_values("diff", ascending=False)
        log.info("Cross-cell-type edges strengthened in NASH: %d", len(gained))
        if len(gained) > 0:
            log.info("Top gained:\n%s", gained.head(10).to_string())

    # =====================================================================
    # PAGA-initialized DPT (global pseudotime)
    # =====================================================================
    log.info("\n--- PAGA-initialized DPT ---")
    # Set root to the most "healthy" cluster
    cluster_health = adata.obs.groupby("leiden_fine").apply(
        lambda g: (g["condition"] == "Healthy").mean()
    )
    root_cluster = cluster_health.idxmax()
    log.info("Root cluster: %s (%.1f%% healthy)", root_cluster,
             100 * cluster_health[root_cluster])

    # Find root cell within root cluster
    root_mask = adata.obs["leiden_fine"] == root_cluster
    root_candidates = np.where(root_mask)[0]
    adata.uns["iroot"] = root_candidates[0]

    sc.tl.dpt(adata)
    dpt_vals = adata.obs["dpt_pseudotime"].values
    valid = np.isfinite(dpt_vals)
    log.info("Global DPT: %d valid cells (%.1f%%)", valid.sum(),
             100 * valid.sum() / len(dpt_vals))

    global_pt = pd.DataFrame(index=adata.obs_names)
    global_pt["paga_global_pseudotime"] = dpt_vals
    global_pt["cell_type"] = adata.obs["cell_type"].values
    global_pt["condition"] = adata.obs["condition"].values
    global_pt["leiden_fine"] = adata.obs["leiden_fine"].values
    global_pt.to_csv(os.path.join(PT_DIR, "paga_global_pseudotime.csv"))

    log.info("\n=== 303: PAGA Cross-Cell-Type COMPLETE ===")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
S3 / 360: Unconstrained Leiden resolution sweep on hepatocyte atlas.

Goal (Pachter §2.3): the original 309 fixed Leiden to 6-10 clusters via
min_frac / max_frac filter, biasing the cluster cardinality. Here we sweep
Leiden over {0.1, 0.15, ..., 2.0} with 20 random seeds per resolution,
NO min_frac / max_frac filter, NO target-cluster constraint.

For each resolution we:
  - Compute consensus matrix across the 20 seeds (fraction of pairs sharing a cluster)
  - Compute per-resolution stability (mean ARI across seed pairs)
  - Identify Progressor-like clusters: cluster whose top markers (computed via
    rank_genes_groups on the modal-seed labels) overlap with S1's pseudobulk
    Progressor markers if available, else fallback to the legacy per-cell
    Progressor markers from `subtype_markers.csv` (TODO_REPLACE_WITH_PSEUDOBULK).
  - Write per-resolution cell-level labels for downstream cross-cohort test.

Input
  /gpfs/.../results_gpu_v2/hepatocyte_subtypes/hepatocyte_atlas_annotated.h5ad
  /gpfs/.../results_gpu_v2/hepatocyte_subtypes/subtype_markers.csv  (legacy)
  /gpfs/.../results_gpu_v2/hepatocyte_subtypes/pseudobulk_markers/Progressor.csv  (S1, optional)

Outputs (worktree path)
  Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/resolution_sweep/
    leiden_labels.parquet            # cell_id x (resolution, seed) long table
    resolution_stability.csv          # resolution, mean_ari, n_clusters_modal
    progressor_clusters.csv           # resolution, seed, cluster_id, n_marker_overlap
    consensus_summary.csv             # resolution, mean_consensus_score
    progressor_per_cell.parquet       # per-cell modal Progressor membership per resolution

Environment: rapids_singlecell (GPU Leiden via rapids_singlecell.tl.leiden).
"""
from __future__ import annotations

import logging
import os
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stdout,
)
log = logging.getLogger("360_resolution_sweep")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
MAIN_REPO = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
WORKTREE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation",
)

H5AD_PATH = os.path.join(
    MAIN_REPO,
    "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_atlas_annotated.h5ad",
)
LEGACY_MARKERS_CSV = os.path.join(
    MAIN_REPO,
    "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/subtype_markers.csv",
)
PSEUDOBULK_MARKERS_CSV = os.path.join(
    WORKTREE,
    "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/pseudobulk_markers/Progressor.csv",
)
# Fallback: 313-retrofit produces published meta-subtype DE under
# .../crossmodal/bulk_sc_convergence/pseudobulk_meta_subtype/subtype_Disease-Progressor_de.csv
# Search both worktree and main project for that file.
PSEUDOBULK_313_CSVS = [
    os.path.join(
        WORKTREE,
        "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/crossmodal/bulk_sc_convergence/pseudobulk_meta_subtype/subtype_Disease-Progressor_de.csv",
    ),
    os.path.join(
        MAIN_REPO,
        "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/crossmodal/bulk_sc_convergence/pseudobulk_meta_subtype/subtype_Disease-Progressor_de.csv",
    ),
]

OUT_DIR = Path(
    os.path.join(
        WORKTREE,
        "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/resolution_sweep",
    )
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Sweep config
# ---------------------------------------------------------------------------
RESOLUTIONS = [round(0.1 + 0.05 * i, 4) for i in range(39)]  # 0.1, 0.15, ..., 2.0  -> 39 values
SEEDS = list(range(20))
TOP_N_MARKERS = 50
MIN_MARKER_OVERLAP = 5  # cluster called "Progressor-like" if >=5 of top-50 markers match

# ---------------------------------------------------------------------------
# Imports that need the env
# ---------------------------------------------------------------------------
import scanpy as sc  # noqa: E402

try:
    import rapids_singlecell as rsc  # noqa: E402
    HAS_RSC = True
except ImportError:
    HAS_RSC = False
    log.warning("rapids_singlecell not importable; falling back to CPU leidenalg.")

# Disable GPU path on Blackwell (sm_100/120) GPUs where rapids/cupy lacks
# precompiled kernels: CUDA_ERROR_NO_BINARY_FOR_GPU at moduleLoadData.
# Set FORCE_CPU_LEIDEN=1 to force the CPU path regardless of GPU detection.
if os.environ.get("FORCE_CPU_LEIDEN", "0") == "1":
    HAS_RSC = False
    log.warning("FORCE_CPU_LEIDEN=1 — using CPU leidenalg.")
else:
    try:
        import cupy
        cupy.zeros((1,), dtype="float32")
    except Exception as _e:
        HAS_RSC = False
        log.warning("cupy GPU init failed (%s); falling back to CPU leidenalg.", _e)

# ARI for stability summary
from sklearn.metrics import adjusted_rand_score  # noqa: E402


def load_progressor_markers() -> list[str]:
    if os.path.exists(PSEUDOBULK_MARKERS_CSV):
        df = pd.read_csv(PSEUDOBULK_MARKERS_CSV)
        # Heuristic column detection
        gene_col = next(
            (c for c in ("gene", "feature", "Gene", "names", "name") if c in df.columns),
            df.columns[0],
        )
        markers = df[gene_col].astype(str).head(TOP_N_MARKERS).tolist()
        log.info(
            "Loaded %d pseudobulk Progressor markers from %s",
            len(markers),
            PSEUDOBULK_MARKERS_CSV,
        )
        return markers
    # First fallback: 313 retrofit's Disease-Progressor meta-subtype DE
    for cand in PSEUDOBULK_313_CSVS:
        if os.path.exists(cand):
            df = pd.read_csv(cand)
            # Filter to up-DEG markers, sort by logFC desc.
            if "is_marker" in df.columns:
                df = df[df["is_marker"] == True]  # noqa: E712
            elif "is_marker_strict" in df.columns:
                df = df[df["is_marker_strict"] == True]  # noqa: E712
            if "logFC" in df.columns:
                df = df.sort_values("logFC", ascending=False)
            gene_col = next(
                (c for c in ("gene", "feature", "Gene", "names", "name") if c in df.columns),
                df.columns[0],
            )
            markers = df[gene_col].astype(str).head(TOP_N_MARKERS).tolist()
            log.info(
                "Loaded %d Disease-Progressor markers from 313 retrofit: %s",
                len(markers), cand,
            )
            return markers
    # Last fallback: legacy subtype_markers.csv contains rank_genes_groups output
    # with one row per (cluster, gene). Filter to clusters labeled Progressor.
    log.warning(
        "Pseudobulk markers missing at %s; falling back to legacy subtype_markers.csv "
        "(TODO_REPLACE_WITH_PSEUDOBULK).",
        PSEUDOBULK_MARKERS_CSV,
    )
    if not os.path.exists(LEGACY_MARKERS_CSV):
        log.error("No legacy markers either; cannot define Progressor.")
        return []
    df = pd.read_csv(LEGACY_MARKERS_CSV)
    # subtype_markers.csv has columns: cluster, gene, pval, logfc, ... (varies)
    col_cluster = next((c for c in df.columns if c.lower() in ("cluster", "group", "subtype")), None)
    col_gene = next((c for c in df.columns if c.lower() in ("gene", "names", "feature", "name")), None)
    if col_cluster is None or col_gene is None:
        log.error("Cannot infer cluster/gene columns from %s; cols=%s", LEGACY_MARKERS_CSV, list(df.columns))
        return []
    # Match any cluster label containing "Progressor" or canonical S2 fate (Progressor == S2)
    mask = df[col_cluster].astype(str).str.contains("Progressor", case=False, na=False)
    if mask.sum() == 0:
        # Try S2 fallback
        mask = df[col_cluster].astype(str).str.contains("S2", case=False, na=False)
    if mask.sum() == 0:
        log.warning("No Progressor / S2 markers in legacy file; returning empty list.")
        return []
    markers = (
        df[mask][col_gene].astype(str).drop_duplicates().head(TOP_N_MARKERS).tolist()
    )
    log.info("Loaded %d legacy Progressor markers (TODO_REPLACE_WITH_PSEUDOBULK).", len(markers))
    return markers


def main():
    t0 = time.time()
    log.info("Loading AnnData from %s", H5AD_PATH)
    adata = sc.read_h5ad(H5AD_PATH)
    log.info("AnnData: %s; obsm keys: %s", adata.shape, list(adata.obsm.keys()))

    if "X_scVI" not in adata.obsm:
        raise RuntimeError("X_scVI not in obsm; cannot reuse latent.")

    # Reuse the existing scVI latent — no scVI retraining (Phase 4 of rigor plan).
    # Build neighbors once on the scVI latent. rapids_singlecell can do it on GPU.
    if HAS_RSC:
        log.info("Computing GPU neighbors on X_scVI ...")
        rsc.get.anndata_to_GPU(adata)
        rsc.pp.neighbors(adata, n_neighbors=30, use_rep="X_scVI")
    else:
        log.info("Computing CPU neighbors on X_scVI ...")
        sc.pp.neighbors(adata, n_neighbors=30, use_rep="X_scVI")

    progressor_markers = load_progressor_markers()
    progressor_markers_set = set(progressor_markers)

    # Long format storage of labels: dict of (resolution, seed) -> int8 labels.
    # We avoid writing every (cell, resolution, seed) row to disk — too large.
    # Instead per-resolution we save (n_cells, n_seeds) int16 arrays as parquet
    # keyed by resolution.
    cell_index = adata.obs.index.to_numpy()
    n_cells = len(cell_index)

    stability_rows = []
    progressor_rows = []
    consensus_rows = []
    per_cell_progressor = {}  # resolution -> bool array (modal Progressor membership)

    for res in RESOLUTIONS:
        log.info("==== Resolution %s ====", res)
        labels_mat = np.full((n_cells, len(SEEDS)), -1, dtype=np.int32)
        for j, seed in enumerate(SEEDS):
            t_l = time.time()
            try:
                if HAS_RSC:
                    rsc.tl.leiden(
                        adata,
                        resolution=float(res),
                        random_state=int(seed),
                        key_added=f"leiden_tmp",
                    )
                else:
                    sc.tl.leiden(
                        adata,
                        resolution=float(res),
                        random_state=int(seed),
                        key_added=f"leiden_tmp",
                        flavor="igraph",          # 10-100× faster than leidenalg
                        directed=False,
                        n_iterations=2,
                    )
                labels_mat[:, j] = adata.obs["leiden_tmp"].astype(int).to_numpy()
            except Exception as e:  # noqa: BLE001
                log.error("Leiden failed res=%s seed=%s: %s", res, seed, e)
                continue
            log.info(
                "  seed=%d done in %.1fs; n_clusters=%d",
                seed,
                time.time() - t_l,
                len(np.unique(labels_mat[:, j])),
            )

        # ARI across seed pairs
        ari_vals = []
        for i, k in combinations(range(len(SEEDS)), 2):
            if (labels_mat[:, i] < 0).any() or (labels_mat[:, k] < 0).any():
                continue
            ari_vals.append(adjusted_rand_score(labels_mat[:, i], labels_mat[:, k]))
        mean_ari = float(np.mean(ari_vals)) if ari_vals else float("nan")

        # Modal n_clusters
        n_clust_per_seed = [
            int(len(np.unique(labels_mat[:, j]))) for j in range(len(SEEDS))
            if (labels_mat[:, j] >= 0).all()
        ]
        modal_n = int(pd.Series(n_clust_per_seed).mode().iloc[0]) if n_clust_per_seed else -1

        # Consensus score: mean over a random sub-sample of 5000 cells of pair-agreement fraction.
        rng = np.random.default_rng(42)
        if n_cells > 5000:
            idx = rng.choice(n_cells, size=5000, replace=False)
            sub = labels_mat[idx]
        else:
            sub = labels_mat
        # For each pair of cells, fraction of seeds where they are co-clustered.
        valid_seeds = [j for j in range(len(SEEDS)) if (labels_mat[:, j] >= 0).all()]
        if len(valid_seeds) >= 2:
            sub_v = sub[:, valid_seeds]
            n_seed_v = len(valid_seeds)
            # vectorized: for each seed, build co-cluster indicator, average across seeds.
            # equality matrix per seed too big; iterate seeds and accumulate.
            cs_n = sub_v.shape[0]
            agree = np.zeros((cs_n, cs_n), dtype=np.float32)
            for s_idx in range(n_seed_v):
                col = sub_v[:, s_idx]
                # outer equality via broadcasting
                eq = (col[:, None] == col[None, :]).astype(np.float32)
                agree += eq
            agree /= n_seed_v
            # off-diagonal mean
            iu = np.triu_indices(cs_n, k=1)
            mean_consensus = float(agree[iu].mean())
        else:
            mean_consensus = float("nan")

        stability_rows.append(
            {
                "resolution": res,
                "mean_ari": mean_ari,
                "modal_n_clusters": modal_n,
                "n_seeds_succeeded": int(len(valid_seeds)),
            }
        )
        consensus_rows.append(
            {"resolution": res, "mean_consensus_score": mean_consensus}
        )

        # Identify Progressor-like clusters per seed via marker overlap.
        # Reuse the per-seed labels: for each seed, compute rank_genes_groups
        # top-50 per cluster, check overlap with progressor_markers.
        # Doing this for every (resolution, seed) is expensive; do it per seed.
        cluster_progressor_calls = np.zeros(n_cells, dtype=np.float32)  # accumulate seed votes
        if progressor_markers_set:
            for j in valid_seeds:
                col = f"leiden_res{res}_seed{SEEDS[j]}"
                adata.obs[col] = pd.Categorical(labels_mat[:, j].astype(str))
                try:
                    sc.tl.rank_genes_groups(
                        adata,
                        col,
                        method="wilcoxon",
                        n_genes=TOP_N_MARKERS,
                        use_raw=False,
                    )
                except Exception as e:  # noqa: BLE001
                    log.warning("rank_genes_groups failed res=%s seed=%s: %s", res, SEEDS[j], e)
                    del adata.obs[col]
                    continue
                names = adata.uns["rank_genes_groups"]["names"]
                clusters = list(names.dtype.names)
                # For each cluster, compute marker overlap with progressor set.
                for cl in clusters:
                    top_genes = set(names[cl][:TOP_N_MARKERS].tolist())
                    overlap = len(top_genes & progressor_markers_set)
                    progressor_rows.append(
                        {
                            "resolution": res,
                            "seed": SEEDS[j],
                            "cluster": cl,
                            "n_top_markers_overlap": overlap,
                        }
                    )
                    if overlap >= MIN_MARKER_OVERLAP:
                        # this cluster is "Progressor-like" in this seed
                        in_cluster = labels_mat[:, j] == int(cl)
                        cluster_progressor_calls[in_cluster] += 1
                del adata.obs[col]
            # cell is "modal Progressor" if >=10/20 seeds called it
            modal_progressor = cluster_progressor_calls >= (len(valid_seeds) / 2.0)
        else:
            modal_progressor = np.zeros(n_cells, dtype=bool)

        per_cell_progressor[res] = modal_progressor

        # Save per-resolution labels matrix as parquet
        df_lab = pd.DataFrame(
            labels_mat,
            columns=[f"seed{s}" for s in SEEDS],
            index=cell_index,
        )
        df_lab.to_parquet(OUT_DIR / f"labels_res{res}.parquet")

        log.info(
            "Res=%s done: ARI=%.3f modal_n=%d consensus=%.3f progressor_cells=%d",
            res,
            mean_ari,
            modal_n,
            mean_consensus,
            int(modal_progressor.sum()),
        )

    # Save stability + consensus + progressor calls
    pd.DataFrame(stability_rows).to_csv(OUT_DIR / "resolution_stability.csv", index=False)
    pd.DataFrame(consensus_rows).to_csv(OUT_DIR / "consensus_summary.csv", index=False)
    pd.DataFrame(progressor_rows).to_csv(OUT_DIR / "progressor_clusters.csv", index=False)

    # Per-cell modal progressor membership across resolutions (long parquet)
    df_pcp = pd.DataFrame(per_cell_progressor, index=cell_index).astype(np.int8)
    df_pcp.index.name = "cell_id"
    df_pcp.reset_index().to_parquet(OUT_DIR / "progressor_per_cell.parquet")

    # Also write cell metadata needed by 361 (dataset, sample, stage)
    cell_meta = adata.obs[
        ["dataset", "sample", "disease_stage_coarse"]
    ].copy()
    cell_meta.index.name = "cell_id"
    cell_meta.reset_index().to_parquet(OUT_DIR / "cell_metadata.parquet")

    log.info("Total runtime: %.1f min", (time.time() - t0) / 60.0)


if __name__ == "__main__":
    main()

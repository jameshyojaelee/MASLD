#!/usr/bin/env python3
"""Spatial projection of the 22 frozen donor-level Hotspot programs.

The primary statistic is residual Moran's I after adjusting each spot score for
the matching cell2location abundance, log library size, and detected genes. A
hepatocyte-zonation-adjusted value is retained as a sensitivity analysis.
Matched-gene nulls additionally match correlation with the corresponding
cell2location lineage so cell-type-derived programs are not compared with a
generic genome-wide background. Null sets are created once per cohort and
reused across physical sections.

GSE192741 is collapsed donor-first (the two H35 sections are technical/physical
replicates) before cohort summaries and its 2-vs-2 disease delta is explicitly
descriptive. Vu remains summarized at ten physical arrays because the local
repository does not contain the barcode-to-biopsy map needed to recover the
paper's 33 biopsies from 32 patients. Spatial graphs are split into connected
tissue islands within each array before Moran statistics are calculated.
"""

from __future__ import annotations

import gc
import hashlib
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


BASE = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
ROOT = BASE / "Analysis/Multimodal_Program_Projection"
SMOKE = os.environ.get("FIG4_SPATIAL_SMOKE", "FALSE").upper() in {"TRUE", "1", "YES"}
SMOKE_DATASET = os.environ.get("FIG4_SPATIAL_SMOKE_DATASET", "GSE192741")
OUT = ROOT / ("results/smoke/spatial" if SMOKE else "results/spatial")
OUT.mkdir(parents=True, exist_ok=True)

REGISTRY_FILE = ROOT / "results/frozen_programs.tsv"
MEMBERSHIP_FILE = ROOT / "results/frozen_program_membership.tsv"
GENE_META_FILE = BASE / "data/gencode_v49_gene_metadata.tsv.gz"
N_NULL = int(os.environ.get("FIG4_SPATIAL_N_NULL", "19" if SMOKE else "9999"))
N_SENSITIVITY_NULL = int(
    os.environ.get("FIG4_SPATIAL_N_SENSITIVITY_NULL", "19" if SMOKE else "999")
)
SEED = 42
GRAPH_DISTANCE_MULTIPLIER = float(os.environ.get("FIG4_SPATIAL_GRAPH_DISTANCE_MULTIPLIER", "2.5"))
MIN_GRAPH_SPOTS = int(os.environ.get("FIG4_SPATIAL_MIN_GRAPH_SPOTS", "8"))

DATASETS = {
    "GSE192741": BASE
    / "Analysis/Spatial/results/cell2location/spatial_model/spatial_deconvolved.h5ad",
    "Vu_et_al_2025": BASE
    / "Analysis/Spatial/results/cell2location/spatial_model_vu/spatial_deconvolved_vu.h5ad",
}

LINEAGE_TO_C2L = {
    "hepatocytes": "Hepatocytes",
    "fibroblasts": "Fibroblasts",
    "macrophages": "Macrophages",
    "cholangiocytes": "Cholangiocytes",
}

PERIPORTAL = ["HAL", "SDS", "ASS1", "CPS1", "ALB", "ASL", "HAMP", "HSD17B13", "GLS2"]
PERICENTRAL = ["CYP2E1", "CYP1A2", "GLUL", "CYP3A4", "CYP2A6"]


def bh(pvalues: pd.Series) -> pd.Series:
    """Benjamini-Hochberg correction preserving missing values."""
    out = pd.Series(np.nan, index=pvalues.index, dtype=float)
    valid = pvalues.notna() & np.isfinite(pvalues)
    p = pvalues.loc[valid].to_numpy(float)
    if not len(p):
        return out
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    vals = np.empty_like(adj)
    vals[order] = np.minimum(adj, 1.0)
    out.loc[valid] = vals
    return out


def dense_counts(adata: ad.AnnData) -> sparse.csr_matrix:
    x = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if sparse.issparse(x):
        return x.tocsr().astype(np.float32)
    return sparse.csr_matrix(np.asarray(x, dtype=np.float32))


def log_normalize(counts: sparse.csr_matrix, library: np.ndarray) -> sparse.csr_matrix:
    scale = np.divide(1e4, library, out=np.zeros_like(library, dtype=float), where=library > 0)
    norm = sparse.diags(scale.astype(np.float32)) @ counts
    norm = norm.tocsr()
    norm.data = np.log1p(norm.data)
    return norm


def standardize_columns(mat: np.ndarray) -> np.ndarray:
    mean = np.nanmean(mat, axis=0, keepdims=True)
    sd = np.nanstd(mat, axis=0, ddof=1, keepdims=True)
    sd[~np.isfinite(sd) | (sd <= 0)] = 1.0
    out = (mat - mean) / sd
    out[~np.isfinite(out)] = 0.0
    return out.astype(np.float32, copy=False)


def extract_z(norm: sparse.csr_matrix, indices: np.ndarray) -> np.ndarray:
    if not len(indices):
        return np.zeros((norm.shape[0], 0), dtype=np.float32)
    arr = norm[:, indices].toarray().astype(np.float32, copy=False)
    return standardize_columns(arr)


def design_matrix(obs: pd.DataFrame, abundance: np.ndarray, zonation: np.ndarray | None = None) -> np.ndarray:
    cols = [
        np.ones(len(obs), dtype=float),
        np.asarray(abundance, dtype=float),
        np.log1p(pd.to_numeric(obs["total_counts"], errors="coerce").to_numpy(float)),
        np.log1p(pd.to_numeric(obs["n_genes_by_counts"], errors="coerce").to_numpy(float)),
    ]
    if zonation is not None:
        cols.append(np.asarray(zonation, dtype=float))
    x = np.column_stack(cols)
    for j in range(1, x.shape[1]):
        mu = np.nanmean(x[:, j])
        sd = np.nanstd(x[:, j], ddof=1)
        x[:, j] = (x[:, j] - mu) / sd if np.isfinite(sd) and sd > 0 else 0.0
    x[~np.isfinite(x)] = 0.0
    return x


def residualize(values: np.ndarray, design: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float32)
    return values - design @ (np.linalg.pinv(design) @ values)


def section_graphs(
    obs: pd.DataFrame,
    coords: np.ndarray,
    k: int = 6,
) -> tuple[dict[str, list[tuple[np.ndarray, np.ndarray, np.ndarray]]], pd.DataFrame]:
    """Build kNN graphs separately within connected tissue islands.

    Visium arrays can contain multiple physically separated biopsy pieces. A
    graph over the whole array creates artificial nearest-neighbor edges across
    empty space. We estimate the within-tissue spot spacing from each array's
    nearest-neighbor distances, remove candidate edges longer than 2.5 times
    that spacing, and compute connected components before constructing the
    final six-neighbor graphs.
    """
    graphs: dict[str, list[tuple[np.ndarray, np.ndarray, np.ndarray]]] = {}
    audit_rows: list[dict] = []
    sample_ids = obs["sample_id"].astype(str).to_numpy()
    for sample in pd.unique(sample_ids):
        idx = np.flatnonzero(sample_ids == str(sample))
        if len(idx) < MIN_GRAPH_SPOTS:
            audit_rows.append(
                {
                    "sample_id": str(sample),
                    "n_spots": len(idx),
                    "n_tissue_islands": 0,
                    "n_spots_in_graph": 0,
                    "n_spots_excluded_small_islands": len(idx),
                    "median_nearest_distance": np.nan,
                    "max_allowed_edge_distance": np.nan,
                    "max_observed_edge_distance": np.nan,
                    "n_edges": 0,
                }
            )
            continue

        local_coords = np.asarray(coords[idx], dtype=float)
        tree = cKDTree(local_coords)
        dist, neigh = tree.query(local_coords, k=min(k + 1, len(idx)))
        dist = np.atleast_2d(dist)
        neigh = np.atleast_2d(neigh)
        nearest = dist[:, 1]
        nearest = nearest[np.isfinite(nearest) & (nearest > 0)]
        if not len(nearest):
            audit_rows.append(
                {
                    "sample_id": str(sample),
                    "n_spots": len(idx),
                    "n_tissue_islands": 0,
                    "n_spots_in_graph": 0,
                    "n_spots_excluded_small_islands": len(idx),
                    "median_nearest_distance": np.nan,
                    "max_allowed_edge_distance": np.nan,
                    "max_observed_edge_distance": np.nan,
                    "n_edges": 0,
                }
            )
            continue
        median_nearest = float(np.median(nearest))
        max_edge_distance = GRAPH_DISTANCE_MULTIPLIER * median_nearest

        candidate_edges: set[tuple[int, int]] = set()
        for i, (drow, nrow) in enumerate(zip(dist, neigh, strict=True)):
            for distance, j in zip(drow[1:], nrow[1:], strict=True):
                if np.isfinite(distance) and distance <= max_edge_distance:
                    a, b = sorted((int(i), int(j)))
                    if a != b:
                        candidate_edges.add((a, b))
        if not candidate_edges:
            audit_rows.append(
                {
                    "sample_id": str(sample),
                    "n_spots": len(idx),
                    "n_tissue_islands": 0,
                    "n_spots_in_graph": 0,
                    "n_spots_excluded_small_islands": len(idx),
                    "median_nearest_distance": median_nearest,
                    "max_allowed_edge_distance": max_edge_distance,
                    "max_observed_edge_distance": np.nan,
                    "n_edges": 0,
                }
            )
            continue

        edge_arr = np.asarray(sorted(candidate_edges), dtype=np.int32)
        adjacency = sparse.coo_matrix(
            (
                np.ones(2 * len(edge_arr), dtype=np.int8),
                (
                    np.concatenate([edge_arr[:, 0], edge_arr[:, 1]]),
                    np.concatenate([edge_arr[:, 1], edge_arr[:, 0]]),
                ),
            ),
            shape=(len(idx), len(idx)),
        ).tocsr()
        _, component = connected_components(adjacency, directed=False, return_labels=True)

        sample_graphs: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
        excluded = 0
        observed_distances: list[float] = []
        for component_id in np.unique(component):
            member_local = np.flatnonzero(component == component_id)
            if len(member_local) < MIN_GRAPH_SPOTS:
                excluded += len(member_local)
                continue
            local_map = np.full(len(idx), -1, dtype=np.int32)
            local_map[member_local] = np.arange(len(member_local), dtype=np.int32)
            keep = np.isin(edge_arr[:, 0], member_local) & np.isin(edge_arr[:, 1], member_local)
            component_edges = edge_arr[keep]
            if not len(component_edges):
                excluded += len(member_local)
                continue
            left = local_map[component_edges[:, 0]]
            right = local_map[component_edges[:, 1]]
            component_idx = idx[member_local]
            sample_graphs.append((component_idx, left, right))
            observed_distances.extend(
                np.linalg.norm(
                    local_coords[component_edges[:, 0]] - local_coords[component_edges[:, 1]],
                    axis=1,
                ).tolist()
            )

        if sample_graphs:
            graphs[str(sample)] = sample_graphs
        audit_rows.append(
            {
                "sample_id": str(sample),
                "n_spots": len(idx),
                "n_tissue_islands": len(sample_graphs),
                "n_spots_in_graph": int(sum(len(x[0]) for x in sample_graphs)),
                "n_spots_excluded_small_islands": excluded,
                "median_nearest_distance": median_nearest,
                "max_allowed_edge_distance": max_edge_distance,
                "max_observed_edge_distance": max(observed_distances) if observed_distances else np.nan,
                "n_edges": int(sum(len(x[1]) for x in sample_graphs)),
            }
        )
    return graphs, pd.DataFrame(audit_rows)


def moran_columns(values: np.ndarray, graphs: dict) -> dict[str, np.ndarray]:
    """Moran I by array after spot-weighted tissue-island aggregation."""
    x = np.asarray(values, dtype=np.float32)
    if x.ndim == 1:
        x = x[:, None]
    out = {}
    for sample, sample_graphs in graphs.items():
        island_values = []
        island_weights = []
        for idx, left, right in sample_graphs:
            z = x[idx]
            z = z - np.mean(z, axis=0, keepdims=True)
            denom = np.sum(z * z, axis=0)
            numer = np.sum(z[left] * z[right], axis=0)
            val = np.divide(
                len(idx) * numer,
                len(left) * denom,
                out=np.full(x.shape[1], np.nan, dtype=float),
                where=(len(left) > 0) & (denom > 0),
            )
            island_values.append(val)
            island_weights.append(len(idx))
        out[sample] = np.average(
            np.vstack(island_values),
            axis=0,
            weights=np.asarray(island_weights, dtype=float),
        )
    return out


def donor_collapse(values_by_sample: dict[str, np.ndarray], obs: pd.DataFrame, dataset: str) -> np.ndarray:
    table = obs[["sample_id", "individual"]].drop_duplicates()
    rows = []
    for individual, part in table.groupby("individual", observed=False):
        arrays = [str(x) for x in part["sample_id"] if str(x) in values_by_sample]
        if arrays:
            rows.append(np.mean(np.vstack([values_by_sample[x] for x in arrays]), axis=0))
    if not rows:
        return np.array([], dtype=float)
    # Vu individuals are intentionally the ten physical arrays; GSE H35 is the
    # only two-section donor and is averaged here before the cohort mean.
    return np.mean(np.vstack(rows), axis=0)


def matched_pool_table(
    var_names: pd.Index,
    mean: np.ndarray,
    detection: np.ndarray,
    biotype: dict[str, str],
) -> pd.DataFrame:
    tbl = pd.DataFrame(
        {
            "gene": var_names.astype(str),
            "mean": mean,
            "detection": detection,
            "biotype": [biotype.get(str(g), "unknown") for g in var_names],
        }
    )
    tbl["mt"] = tbl["gene"].str.startswith("MT-")
    tbl["ribo"] = tbl["gene"].str.match(r"^RP[SL]")
    tbl["mean_bin"] = pd.qcut(tbl["mean"].rank(method="first"), 10, labels=False, duplicates="drop")
    tbl["detect_bin"] = pd.qcut(
        tbl["detection"].rank(method="first"), 10, labels=False, duplicates="drop"
    )
    return tbl


def abundance_correlation(norm: sparse.csr_matrix, abundance: np.ndarray) -> np.ndarray:
    """Pearson correlation of every gene with a cell2location abundance vector."""
    a = np.asarray(abundance, dtype=float)
    a = a - np.nanmean(a)
    a[~np.isfinite(a)] = 0.0
    abundance_ss = float(np.sum(a * a))
    gene_mean = np.asarray(norm.mean(axis=0)).ravel()
    gene_sq_mean = np.asarray(norm.multiply(norm).mean(axis=0)).ravel()
    gene_ss = norm.shape[0] * np.maximum(gene_sq_mean - gene_mean * gene_mean, 0.0)
    numerator = np.asarray(norm.T @ a).ravel()
    denominator = np.sqrt(gene_ss * abundance_ss)
    return np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator),
        where=denominator > 0,
    )


def build_matched_sets(
    measured: pd.DataFrame,
    pool: pd.DataFrame,
    program_gene_union: set[str],
    n_null: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, list[int], list[int]]:
    pool_idx = pool.set_index("gene", drop=False)
    candidate_map: list[np.ndarray] = []
    expression_relaxation: list[int] = []
    lineage_relaxation: list[int] = []
    excluded = set(program_gene_union)
    for gene in measured["gene_symbol"]:
        target = pool_idx.loc[gene]
        candidates = np.array([], dtype=int)
        used_expression_relax = 3
        used_lineage_relax = 3
        for lineage_relax in (0, 1, 2):
            for expression_relax in (0, 1, 2):
                mask = (
                    (pool["biotype"] == target["biotype"])
                    & (pool["mt"] == target["mt"])
                    & (pool["ribo"] == target["ribo"])
                    & ((pool["mean_bin"] - target["mean_bin"]).abs() <= expression_relax)
                    & ((pool["detect_bin"] - target["detect_bin"]).abs() <= expression_relax)
                    & (
                        (pool["lineage_corr_bin"] - target["lineage_corr_bin"]).abs()
                        <= lineage_relax
                    )
                    & (~pool["gene"].isin(excluded))
                )
                candidates = np.flatnonzero(mask.to_numpy())
                if len(candidates) >= 10:
                    used_expression_relax = expression_relax
                    used_lineage_relax = lineage_relax
                    break
            if len(candidates) >= 10:
                break
        # A rare biotype may have only one eligible strict control. Retaining
        # that control is preferable to relaxing the lineage window: its gene
        # contribution is fixed across null sets, while the other program
        # genes still generate the combinatorial matched-set distribution.
        if not len(candidates):
            mask = (
                (pool["biotype"] == target["biotype"])
                & ((pool["lineage_corr_bin"] - target["lineage_corr_bin"]).abs() <= 2)
                & (~pool["gene"].isin(excluded))
            )
            candidates = np.flatnonzero(mask.to_numpy())
            used_expression_relax = 3
            used_lineage_relax = 2
        if not len(candidates):
            mask = (pool["biotype"] == target["biotype"]) & (~pool["gene"].isin(excluded))
            candidates = np.flatnonzero(mask.to_numpy())
            used_expression_relax = 3
            used_lineage_relax = 3
        if not len(candidates):
            raise RuntimeError(f"No matched spatial control genes for {gene}")
        candidate_map.append(candidates)
        expression_relaxation.append(used_expression_relax)
        lineage_relaxation.append(used_lineage_relax)

    sets = np.empty((n_null, len(candidate_map)), dtype=np.int32)
    for i in range(n_null):
        chosen: set[int] = set()
        for j, candidates in enumerate(candidate_map):
            available = candidates[~np.isin(candidates, np.fromiter(chosen, dtype=int))]
            source = available if len(available) else candidates
            pick = int(rng.choice(source))
            sets[i, j] = pick
            chosen.add(pick)
    return sets, expression_relaxation, lineage_relaxation


def null_moran(
    norm: sparse.csr_matrix,
    matched_sets: np.ndarray,
    weights: np.ndarray,
    design: np.ndarray,
    design_zonation: np.ndarray,
    graphs: dict,
    obs: pd.DataFrame,
    dataset: str,
    chunk_size: int = 64,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    raw_out = []
    residual_out = []
    zonation_out = []
    for start in range(0, matched_sets.shape[0], chunk_size):
        block = matched_sets[start : start + chunk_size]
        uniq, inverse = np.unique(block, return_inverse=True)
        z = extract_z(norm, uniq)
        mapped = inverse.reshape(block.shape)
        scores = np.zeros((norm.shape[0], len(block)), dtype=np.float32)
        for j, weight in enumerate(weights):
            scores += z[:, mapped[:, j]] * np.float32(weight)
        resid = residualize(scores, design)
        resid_z = residualize(scores, design_zonation)
        raw_by = moran_columns(scores, graphs)
        res_by = moran_columns(resid, graphs)
        zon_by = moran_columns(resid_z, graphs)
        raw_out.append(donor_collapse(raw_by, obs, dataset))
        residual_out.append(donor_collapse(res_by, obs, dataset))
        zonation_out.append(donor_collapse(zon_by, obs, dataset))
    return (
        np.concatenate(raw_out),
        np.concatenate(residual_out),
        np.concatenate(zonation_out),
    )


def null_residual_moran(
    norm: sparse.csr_matrix,
    matched_sets: np.ndarray,
    weights: np.ndarray,
    design: np.ndarray,
    graphs: dict,
    obs: pd.DataFrame,
    dataset: str,
    chunk_size: int = 64,
) -> np.ndarray:
    """Residual Moran null for a sensitivity-specific weighting scheme."""
    residual_out = []
    for start in range(0, matched_sets.shape[0], chunk_size):
        block = matched_sets[start : start + chunk_size]
        uniq, inverse = np.unique(block, return_inverse=True)
        z = extract_z(norm, uniq)
        mapped = inverse.reshape(block.shape)
        scores = np.zeros((norm.shape[0], len(block)), dtype=np.float32)
        for j, weight in enumerate(weights):
            scores += z[:, mapped[:, j]] * np.float32(weight)
        residual_by = moran_columns(residualize(scores, design), graphs)
        residual_out.append(donor_collapse(residual_by, obs, dataset))
    return np.concatenate(residual_out)


def standardized_against_null(observed: float, null: np.ndarray) -> tuple[float, float, float, float]:
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    minimum_null = 10 if SMOKE else 100
    if len(null) < minimum_null or not np.isfinite(observed):
        return np.nan, np.nan, np.nan, np.nan
    mu = float(np.mean(null))
    sd = float(np.std(null, ddof=1))
    z = (observed - mu) / sd if sd > 0 else np.nan
    p = (1.0 + float(np.sum(null >= observed))) / (len(null) + 1.0)
    return z, p, mu, sd


def process_dataset(
    dataset: str,
    path: Path,
    registry: pd.DataFrame,
    membership: pd.DataFrame,
    biotype: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    print(f"[spatial] loading {dataset}: {path}", flush=True)
    adata = ad.read_h5ad(path)
    adata.var_names = adata.var_names.astype(str)
    obs = adata.obs.copy()
    obs["sample_id"] = obs["sample_id"].astype(str)
    obs["individual"] = obs["individual"].astype(str)
    counts = dense_counts(adata)
    library = pd.to_numeric(obs["total_counts"], errors="coerce").to_numpy(float)
    norm = log_normalize(counts, library)
    mean = np.asarray(norm.mean(axis=0)).ravel()
    detection = np.asarray((counts > 0).mean(axis=0)).ravel()
    pool = matched_pool_table(adata.var_names, mean, detection, biotype)
    gene_to_index = {g: i for i, g in enumerate(adata.var_names)}
    program_gene_union = set(membership["gene_symbol"].dropna().astype(str))

    factor_names = [str(x).replace("means_per_cluster_mu_fg_", "") for x in adata.uns["mod"]["factor_names"]]
    abundance_mat = np.asarray(adata.obsm["q05_cell_abundance_w_sf"])
    abundance_by_lineage = {
        lineage: abundance_mat[:, factor_names.index(label)]
        for lineage, label in LINEAGE_TO_C2L.items()
    }
    lineage_pools = {}
    for lineage, abundance in abundance_by_lineage.items():
        lineage_pool = pool.copy()
        lineage_pool["lineage_corr"] = abundance_correlation(norm, abundance)
        lineage_pool["lineage_corr_bin"] = pd.qcut(
            lineage_pool["lineage_corr"].rank(method="first"),
            10,
            labels=False,
            duplicates="drop",
        )
        lineage_pools[lineage] = lineage_pool

    zon_genes = [g for g in PERIPORTAL + PERICENTRAL if g in gene_to_index]
    zon_z = extract_z(norm, np.array([gene_to_index[g] for g in zon_genes], dtype=int))
    zmap = {g: zon_z[:, i] for i, g in enumerate(zon_genes)}
    pp = np.mean(np.column_stack([zmap[g] for g in PERIPORTAL if g in zmap]), axis=1)
    pc = np.mean(np.column_stack([zmap[g] for g in PERICENTRAL if g in zmap]), axis=1)
    zonation = pc - pp

    coords = np.asarray(adata.obsm["spatial"])
    graphs, graph_audit = section_graphs(obs, coords, k=6)
    graph_audit.insert(0, "dataset", dataset)
    expected_arrays = set(obs["sample_id"].astype(str))
    missing_graphs = sorted(expected_arrays.difference(graphs))
    if missing_graphs:
        raise RuntimeError(
            f"No valid tissue-island spatial graph for {dataset} array(s): "
            + ", ".join(missing_graphs)
        )
    rng = np.random.default_rng(SEED + (0 if dataset == "GSE192741" else 1000))

    result_rows = []
    section_rows = []
    match_rows = []
    null_summary_rows = []
    for row in registry.itertuples(index=False):
        pid = row.program_id
        m = membership[membership["program_id"] == pid].copy()
        m = m.groupby("gene_symbol", as_index=False)["original_l1_weight"].sum()
        m["present"] = m["gene_symbol"].isin(gene_to_index)
        m["detected"] = m["gene_symbol"].map(
            lambda g: detection[gene_to_index[g]] >= 0.01 if g in gene_to_index else False
        )
        measured = m[m["present"] & m["detected"]].copy()
        retained = float(measured["original_l1_weight"].sum())
        testable = len(measured) >= 8 and retained >= 0.20
        base = {
            "program_id": pid,
            "dataset": dataset,
            "display_order": row.display_order,
            "cell_type": row.cell_type,
            "module": row.module,
            "program_name": row.program_name,
            "n_measured": len(measured),
            "retained_l1_weight": retained,
            "testable": testable,
            "n_null": N_NULL,
        }
        if not testable:
            result_rows.append(base)
            continue

        measured["weight"] = measured["original_l1_weight"] / retained
        indices = np.array([gene_to_index[g] for g in measured["gene_symbol"]], dtype=int)
        weights = measured["weight"].to_numpy(float)
        z = extract_z(norm, indices)
        score = z @ weights
        equal_score = np.mean(z, axis=1)
        top = int(np.argmax(weights))
        leave_weights = np.delete(weights, top)
        leave_weights = leave_weights / leave_weights.sum()
        leave_score = np.delete(z, top, axis=1) @ leave_weights

        design = design_matrix(obs, abundance_by_lineage[row.cell_type])
        design_zonation = design_matrix(obs, abundance_by_lineage[row.cell_type], zonation=zonation)
        residual = residualize(score, design)
        residual_zonation = residualize(score, design_zonation)
        equal_residual = residualize(equal_score, design)
        leave_residual = residualize(leave_score, design)

        raw_by = moran_columns(score, graphs)
        residual_by = moran_columns(residual, graphs)
        zonation_by = moran_columns(residual_zonation, graphs)
        equal_by = moran_columns(equal_residual, graphs)
        leave_by = moran_columns(leave_residual, graphs)

        observed_raw = float(donor_collapse(raw_by, obs, dataset)[0])
        observed_residual = float(donor_collapse(residual_by, obs, dataset)[0])
        observed_zonation = float(donor_collapse(zonation_by, obs, dataset)[0])
        observed_equal = float(donor_collapse(equal_by, obs, dataset)[0])
        observed_leave = float(donor_collapse(leave_by, obs, dataset)[0])

        lineage_pool = lineage_pools[row.cell_type]
        matched_sets, expression_relax, lineage_relax = build_matched_sets(
            measured, lineage_pool, program_gene_union, N_NULL, rng
        )
        matched_set_sha256 = hashlib.sha256(matched_sets.tobytes()).hexdigest()
        null_raw, null_residual, null_zonation = null_moran(
            norm,
            matched_sets,
            weights,
            design,
            design_zonation,
            graphs,
            obs,
            dataset,
        )
        sensitivity_sets = matched_sets[: min(N_SENSITIVITY_NULL, len(matched_sets))]
        equal_null_residual = null_residual_moran(
            norm,
            sensitivity_sets,
            np.repeat(1.0 / len(weights), len(weights)),
            design,
            graphs,
            obs,
            dataset,
        )
        leave_null_residual = null_residual_moran(
            norm,
            np.delete(sensitivity_sets, top, axis=1),
            leave_weights,
            design,
            graphs,
            obs,
            dataset,
        )
        equal_null_mean = float(np.nanmean(equal_null_residual))
        leave_null_mean = float(np.nanmean(leave_null_residual))
        raw_z, raw_p, raw_mu, raw_sd = standardized_against_null(observed_raw, null_raw)
        res_z, res_p, res_mu, res_sd = standardized_against_null(
            observed_residual, null_residual
        )
        zon_zscore, zon_p, zon_mu, zon_sd = standardized_against_null(
            observed_zonation, null_zonation
        )
        for statistic, null_values in (
            ("raw_moran_i", null_raw),
            ("residual_moran_i", null_residual),
            ("zonation_moran_i", null_zonation),
        ):
            finite_null = np.asarray(null_values, dtype=float)
            finite_null = finite_null[np.isfinite(finite_null)]
            quantiles = np.quantile(finite_null, [0.001, 0.01, 0.05, 0.50, 0.95, 0.99, 0.999])
            null_summary_rows.append(
                {
                    "program_id": pid,
                    "dataset": dataset,
                    "statistic": statistic,
                    "n_null": len(finite_null),
                    "matched_set_sha256": matched_set_sha256,
                    "null_mean": float(np.mean(finite_null)),
                    "null_sd": float(np.std(finite_null, ddof=1)),
                    "q001": quantiles[0],
                    "q010": quantiles[1],
                    "q050": quantiles[2],
                    "q500": quantiles[3],
                    "q950": quantiles[4],
                    "q990": quantiles[5],
                    "q999": quantiles[6],
                }
            )

        sample_mean = pd.DataFrame(
            {
                "sample_id": obs["sample_id"].to_numpy(),
                "individual": obs["individual"].to_numpy(),
                "condition": obs["condition"].astype(str).to_numpy(),
                "residual_score": residual,
            }
        ).groupby(["sample_id", "individual", "condition"], as_index=False)["residual_score"].mean()
        donor_mean = sample_mean.groupby(["individual", "condition"], as_index=False)["residual_score"].mean()
        disease_delta = np.nan
        if dataset == "GSE192741":
            healthy = donor_mean.loc[donor_mean["condition"] == "Healthy", "residual_score"]
            steatotic = donor_mean.loc[donor_mean["condition"] == "Steatotic", "residual_score"]
            if len(healthy) == 2 and len(steatotic) == 2:
                disease_delta = float(steatotic.mean() - healthy.mean())

        result_rows.append(
            {
                **base,
                "raw_moran_i": observed_raw,
                "raw_null_mean": raw_mu,
                "raw_null_sd": raw_sd,
                "raw_moran_z": raw_z,
                "raw_pvalue": raw_p,
                "residual_moran_i": observed_residual,
                "residual_null_mean": res_mu,
                "residual_null_sd": res_sd,
                "residual_moran_z": res_z,
                "residual_pvalue": res_p,
                "zonation_moran_i": observed_zonation,
                "zonation_null_mean": zon_mu,
                "zonation_null_sd": zon_sd,
                "zonation_moran_z": zon_zscore,
                "zonation_pvalue": zon_p,
                "equal_residual_moran_i": observed_equal,
                "equal_residual_null_mean": equal_null_mean,
                "leave_top_residual_moran_i": observed_leave,
                "leave_top_residual_null_mean": leave_null_mean,
                "disease_delta_descriptive": disease_delta,
                "matched_set_sha256": matched_set_sha256,
                "target_lineage_corr_mean": float(
                    np.mean(lineage_pool["lineage_corr"].to_numpy()[indices])
                ),
            }
        )
        for sample in graphs:
            meta_row = obs.loc[obs["sample_id"] == sample].iloc[0]
            section_rows.append(
                {
                    "program_id": pid,
                    "dataset": dataset,
                    "sample_id": sample,
                    "individual": str(meta_row["individual"]),
                    "condition": str(meta_row["condition"]),
                    "raw_moran_i": float(raw_by[sample][0]),
                    "residual_moran_i": float(residual_by[sample][0]),
                    "zonation_moran_i": float(zonation_by[sample][0]),
                }
            )
        for gene_index, (gene, expression_level, _lineage_level) in enumerate(
            zip(measured["gene_symbol"], expression_relax, lineage_relax, strict=True)
        ):
            target = lineage_pool.loc[lineage_pool["gene"] == gene].iloc[0]
            control_rows = lineage_pool.iloc[matched_sets[:, gene_index]]
            lineage_bin_difference = np.abs(
                control_rows["lineage_corr_bin"].to_numpy(float)
                - float(target["lineage_corr_bin"])
            )
            match_rows.append(
                {
                    "program_id": pid,
                    "dataset": dataset,
                    "gene_symbol": gene,
                    "expression_match_relaxation": expression_level,
                    # The realized sampled controls are authoritative. A pool
                    # smaller than the preferred ten can remain strictly
                    # matched even when the search initializer was not reset.
                    "lineage_match_relaxation": int(
                        np.max(lineage_bin_difference)
                    ),
                    "target_lineage_correlation": target["lineage_corr"],
                    "target_lineage_correlation_bin": target["lineage_corr_bin"],
                    "mean_control_lineage_correlation": float(control_rows["lineage_corr"].mean()),
                    "n_unique_control_genes": int(
                        np.unique(matched_sets[:, gene_index]).size
                    ),
                    "mean_absolute_lineage_bin_difference": float(
                        np.mean(lineage_bin_difference)
                    ),
                    "max_absolute_lineage_bin_difference": float(
                        np.max(lineage_bin_difference)
                    ),
                }
            )
        print(
            f"[spatial] {dataset} {pid}: n={len(measured)} weight={retained:.3f} "
            f"residual Z={res_z:.2f} p={res_p:.4g}",
            flush=True,
        )

    results = pd.DataFrame(result_rows)
    results["residual_padj"] = bh(results["residual_pvalue"])
    results["raw_padj"] = bh(results["raw_pvalue"])
    results["zonation_padj"] = bh(results["zonation_pvalue"])
    results["sensitivity_sign_agree"] = (
        results["testable"].fillna(False)
        & np.isfinite(results["residual_moran_z"])
        & (
            np.sign(results["residual_moran_z"])
            == np.sign(
                results["equal_residual_moran_i"]
                - results["equal_residual_null_mean"]
            )
        )
        & (
            np.sign(results["residual_moran_z"])
            == np.sign(
                results["leave_top_residual_moran_i"]
                - results["leave_top_residual_null_mean"]
            )
        )
    )
    results["robust"] = (
        results["testable"].fillna(False)
        & (results["residual_padj"] < 0.05)
        & results["sensitivity_sign_agree"].fillna(False)
    )

    del adata, counts, norm
    gc.collect()
    return (
        results,
        pd.DataFrame(section_rows),
        pd.DataFrame(match_rows),
        graph_audit,
        pd.DataFrame(null_summary_rows),
    )


def main() -> None:
    required = [REGISTRY_FILE, MEMBERSHIP_FILE, GENE_META_FILE, *DATASETS.values()]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing spatial input(s): " + ", ".join(missing))
    if N_NULL < 9999:
        print(
            f"[spatial] WARNING: FIG4_SPATIAL_N_NULL={N_NULL} is a smoke-test override; "
            "canonical production requires >=9999",
            flush=True,
        )

    registry = pd.read_csv(REGISTRY_FILE, sep="\t")
    if SMOKE:
        registry = registry.iloc[:1].copy()
    membership = pd.read_csv(MEMBERSHIP_FILE, sep="\t")
    membership = membership[
        membership["mapped_symbol"].fillna(False)
        & membership["gene_symbol"].notna()
    ].copy()
    gene_meta = pd.read_csv(
        GENE_META_FILE,
        sep="\t",
        usecols=["gene_name", "gene_biotype"],
        compression="gzip",
    ).drop_duplicates("gene_name")
    biotype = dict(zip(gene_meta["gene_name"].astype(str), gene_meta["gene_biotype"].astype(str)))

    all_results = []
    all_sections = []
    all_matches = []
    all_graph_audits = []
    all_null_summaries = []
    if SMOKE:
        if SMOKE_DATASET not in DATASETS:
            raise ValueError(
                f"Unknown FIG4_SPATIAL_SMOKE_DATASET={SMOKE_DATASET}; expected one of {list(DATASETS)}"
            )
        datasets_to_run = {SMOKE_DATASET: DATASETS[SMOKE_DATASET]}
    else:
        datasets_to_run = DATASETS
    for dataset, path in datasets_to_run.items():
        res, sec, match, graph_audit, null_summary = process_dataset(
            dataset, path, registry, membership, biotype
        )
        all_results.append(res)
        all_sections.append(sec)
        all_matches.append(match)
        all_graph_audits.append(graph_audit)
        all_null_summaries.append(null_summary)

    results = pd.concat(all_results, ignore_index=True)
    sections = pd.concat(all_sections, ignore_index=True)
    matches = pd.concat(all_matches, ignore_index=True)
    graph_audits = pd.concat(all_graph_audits, ignore_index=True)
    null_summaries = pd.concat(all_null_summaries, ignore_index=True)
    results.to_csv(OUT / "spatial_program_results.tsv", sep="\t", index=False)
    sections.to_csv(OUT / "spatial_section_results.tsv", sep="\t", index=False)
    matches.to_csv(OUT / "spatial_matching_audit.tsv", sep="\t", index=False)
    graph_audits.to_csv(OUT / "spatial_graph_audit.tsv", sep="\t", index=False)
    null_summaries.to_csv(OUT / "spatial_null_summary.tsv", sep="\t", index=False)
    print("[spatial] testable/robust by dataset", flush=True)
    print(results.groupby("dataset")[["testable", "robust"]].sum(), flush=True)


if __name__ == "__main__":
    main()

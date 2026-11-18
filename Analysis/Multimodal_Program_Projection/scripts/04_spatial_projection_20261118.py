#!/usr/bin/env python3
"""Visium Spatial Projection Engine (20261118 candidate, N08).

Key methodological updates over 04_spatial_projection.py:
1. Within-section nonlinear covariate residualization:
   For each physical `sample_id` section s separately, residualize spot log1p-CPM10k
   gene expression on
     X_s = [1, log1p(q05_s), q05_s, log1p(total_counts_s), log1p(n_genes_s)]
   using ridge regression with unpenalized intercept (Lambda = diag(0, 1e-6, ..., 1e-6)),
   then standardize residuals to zero mean and unit variance within section s.
2. Single-gene mean Moran's I statistic:
   Precompute single-gene Moran's I_s(z_tilde_{s,g}) across all detected genes within each
   section's connected tissue islands, collapse across sections per donor/unit (I_donor,
   cohort mean bar_I), and evaluate the weighted single-gene mean Moran's I:
     I_P = sum_{g in P} w_g * bar_I_g
   against 9,999 matched-gene draws across all 22 v1 AND all 117 frozen v2 Hotspot
   programs (GENCODE v49 primary-assembly membership). Eliminates composite-score
   off-diagonal shot-noise cancellation and reduces runtime to seconds.
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
OUT = Path(
    os.environ.get(
        "SPATIAL_OUT_20261118",
        str(ROOT / "candidates/spatial_20261118"),
    )
)

V1_REGISTRY_FILE = ROOT / "results/frozen_programs.tsv"
V1_MEMBERSHIP_FILE = ROOT / "results/frozen_program_membership.tsv"
V2_REGISTRY_FILE = (
    BASE
    / "figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv"
)
V2_MEMBERSHIP_FILE = (
    ROOT
    / "candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/program_membership_v2.tsv"
)
GENE_META_FILE = BASE / "data/gencode_v49_gene_metadata.tsv.gz"

N_NULL = int(os.environ.get("FIG4_SPATIAL_N_NULL", "9999"))
N_SENSITIVITY_NULL = int(os.environ.get("FIG4_SPATIAL_N_SENSITIVITY_NULL", "999"))
SEED = 42
RIDGE_LAMBDA = 1e-6
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
    "endothelial": "Endothelial cells",
    "tcells": "T cells",
}


def nan_weighted_dot(vals_1d: np.ndarray, w_1d: np.ndarray) -> float:
    ok = np.isfinite(vals_1d)
    sw = float(w_1d[ok].sum())
    return float(np.dot(w_1d[ok], vals_1d[ok]) / sw) if sw > 0 else np.nan


def nan_weighted_mat(vals_2d: np.ndarray, w_1d: np.ndarray) -> np.ndarray:
    ok = np.isfinite(vals_2d)
    w_ok = np.where(ok, w_1d[None, :], 0.0)
    sw = w_ok.sum(axis=1)
    num = np.where(ok, vals_2d * w_1d[None, :], 0.0).sum(axis=1)
    return np.divide(num, sw, out=np.full(vals_2d.shape[0], np.nan, dtype=np.float64), where=sw > 0)

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


def section_graphs(
    obs: pd.DataFrame,
    coords: np.ndarray,
    k: int = 6,
) -> tuple[dict[str, list[tuple[np.ndarray, np.ndarray, np.ndarray]]], pd.DataFrame]:
    """Build kNN graphs separately within connected tissue islands per physical section."""
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


def _standardize_within_section_design(cols: list[np.ndarray]) -> np.ndarray:
    x = np.column_stack(cols).astype(np.float64, copy=False)
    for j in range(1, x.shape[1]):
        col = x[:, j]
        mu = np.nanmean(col)
        sd = np.nanstd(col, ddof=1)
        x[:, j] = (col - mu) / sd if (np.isfinite(sd) and sd > 0) else 0.0
    x[~np.isfinite(x)] = 0.0
    return x


def _standardize_matrix_cols(mat: np.ndarray) -> np.ndarray:
    mu = np.mean(mat, axis=0, keepdims=True)
    centered = mat - mu
    sd = np.std(centered, axis=0, ddof=1, keepdims=True)
    valid = np.isfinite(sd) & (sd > 1e-12)
    out = np.divide(centered, np.where(valid, sd, 1.0))
    out[:, ~valid.ravel()] = 0.0
    return out


def _section_moran_from_z(
    z_sec: np.ndarray,
    idx_sec_global: np.ndarray,
    sample_graphs: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
) -> np.ndarray:
    """Compute vector of single-gene Moran's I (length G) for one physical section."""
    g2l = np.full(int(idx_sec_global.max()) + 1, -1, dtype=np.int32)
    g2l[idx_sec_global] = np.arange(len(idx_sec_global), dtype=np.int32)
    isl_vals = []
    isl_wts = []
    n_genes = z_sec.shape[1]
    for isl_global, left, right in sample_graphs:
        loc = g2l[isl_global]
        zi = z_sec[loc, :]
        zi = zi - zi.mean(axis=0, keepdims=True)
        denom = np.sum(zi * zi, axis=0)
        numer = np.sum(zi[left, :] * zi[right, :], axis=0)
        val = np.divide(
            len(loc) * numer,
            len(left) * denom,
            out=np.full(n_genes, np.nan, dtype=np.float64),
            where=(len(left) > 0) & (denom > 0),
        )
        isl_vals.append(val)
        isl_wts.append(len(loc))
    vals_stack = np.vstack(isl_vals)
    wts = np.asarray(isl_wts, dtype=np.float64)
    # Spot-weighted mean across islands (nan-aware per gene)
    valid = np.isfinite(vals_stack)
    w_mat = np.where(valid, wts[:, None], 0.0)
    w_sum = w_mat.sum(axis=0)
    v_sum = np.where(valid, vals_stack * w_mat, 0.0).sum(axis=0)
    return np.divide(v_sum, w_sum, out=np.full(n_genes, np.nan, dtype=np.float64), where=w_sum > 0)


def compute_within_section_single_gene_moran(
    norm: sparse.csr_matrix,
    obs: pd.DataFrame,
    graphs: dict[str, list[tuple[np.ndarray, np.ndarray, np.ndarray]]],
    q05_matrix: np.ndarray,
    extra_cov: np.ndarray | None = None,
    lam: float = RIDGE_LAMBDA,
    compute_raw: bool = True,
) -> dict[str, object]:
    """Residualize each section on [1, log1p(q05_s), q05_s, log1p(total_counts_s), log1p(n_genes_s)]
    via ridge (lambda = 1e-6 on non-intercept columns), standardize within section, and
    precompute single-gene Moran's I across all detected genes.
    """
    q05 = np.asarray(q05_matrix, dtype=np.float64)
    if q05.ndim == 1:
        q05 = q05[:, None]
    if np.any(q05 < -1e-12) or not np.isfinite(q05).all():
        raise ValueError("q05_matrix must be finite and nonnegative for log1p(q05_s) design")
    q05 = np.maximum(q05, 0.0)

    total_counts = np.maximum(pd.to_numeric(obs["total_counts"], errors="coerce").to_numpy(np.float64), 0.0)
    n_genes = np.maximum(pd.to_numeric(obs["n_genes_by_counts"], errors="coerce").to_numpy(np.float64), 0.0)
    sample_ids = obs["sample_id"].astype(str).to_numpy()
    ind_by_sample = obs[["sample_id", "individual"]].drop_duplicates().set_index("sample_id")["individual"].astype(str).to_dict()

    n_genes_total = norm.shape[1]
    sec_resid: dict[str, np.ndarray] = {}
    sec_raw: dict[str, np.ndarray] = {}
    arrays = sorted(graphs.keys())

    for sample in arrays:
        idx_s = np.flatnonzero(sample_ids == sample)
        Y_s = norm[idx_s, :].toarray().astype(np.float64, copy=False)

        if compute_raw:
            Z_raw = _standardize_matrix_cols(Y_s)
            sec_raw[sample] = _section_moran_from_z(Z_raw, idx_s, graphs[sample])

        cols = [np.ones(len(idx_s), dtype=np.float64)]
        for c_idx in range(q05.shape[1]):
            qs = q05[idx_s, c_idx]
            cols.append(np.log1p(qs))
            cols.append(qs)
        cols.append(np.log1p(total_counts[idx_s]))
        cols.append(np.log1p(n_genes[idx_s]))
        if extra_cov is not None:
            ec = np.asarray(extra_cov, dtype=np.float64)
            if ec.ndim == 1:
                ec = ec[:, None]
            for j in range(ec.shape[1]):
                cols.append(ec[idx_s, j])

        X_s = _standardize_within_section_design(cols)
        p = X_s.shape[1]
        penalty = np.zeros(p, dtype=np.float64)
        penalty[1:] = lam
        xtx_reg = X_s.T @ X_s + np.diag(penalty)
        beta_s = np.linalg.solve(xtx_reg, X_s.T @ Y_s)
        R_s = Y_s - X_s @ beta_s
        Z_s = _standardize_matrix_cols(R_s)
        sec_resid[sample] = _section_moran_from_z(Z_s, idx_s, graphs[sample])

    # Collapse per donor/unit, then cohort mean across donors
    donors = sorted({ind_by_sample[s] for s in arrays})
    donor_mat_resid = np.empty((len(donors), n_genes_total), dtype=np.float64)
    donor_mat_raw = np.empty((len(donors), n_genes_total), dtype=np.float64) if compute_raw else None
    for d_idx, donor in enumerate(donors):
        d_arrays = [s for s in arrays if ind_by_sample[s] == donor]
        donor_mat_resid[d_idx, :] = np.nanmean(np.vstack([sec_resid[s] for s in d_arrays]), axis=0)
        if compute_raw:
            donor_mat_raw[d_idx, :] = np.nanmean(np.vstack([sec_raw[s] for s in d_arrays]), axis=0)

    bar_I_resid = np.nanmean(donor_mat_resid, axis=0)
    bar_I_raw = np.nanmean(donor_mat_raw, axis=0) if compute_raw else None

    return {
        "arrays": arrays,
        "donors": donors,
        "sec_resid": sec_resid,
        "sec_raw": sec_raw if compute_raw else None,
        "donor_mat_resid": donor_mat_resid,
        "donor_mat_raw": donor_mat_raw,
        "bar_I_resid": bar_I_resid,
        "bar_I_raw": bar_I_raw,
    }


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
    tbl["detect_bin"] = pd.qcut(tbl["detection"].rank(method="first"), 10, labels=False, duplicates="drop")
    return tbl


def abundance_correlation(norm: sparse.csr_matrix, abundance: np.ndarray) -> np.ndarray:
    a = np.asarray(abundance, dtype=float)
    a = a - np.nanmean(a)
    a[~np.isfinite(a)] = 0.0
    abundance_ss = float(np.sum(a * a))
    gene_mean = np.asarray(norm.mean(axis=0)).ravel()
    gene_sq_mean = np.asarray(norm.multiply(norm).mean(axis=0)).ravel()
    gene_ss = norm.shape[0] * np.maximum(gene_sq_mean - gene_mean * gene_mean, 0.0)
    numerator = np.asarray(norm.T @ a).ravel()
    denominator = np.sqrt(gene_ss * abundance_ss)
    return np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)


def build_matched_sets(
    measured: pd.DataFrame,
    pool: pd.DataFrame,
    program_gene_union: set[str],
    n_null: int,
    rng: np.random.Generator,
    finite_mask: np.ndarray | None = None,
) -> tuple[np.ndarray, list[int], list[int]]:
    pool_idx = pool.set_index("gene", drop=False)
    candidate_map: list[np.ndarray] = []
    expression_relaxation: list[int] = []
    lineage_relaxation: list[int] = []
    excluded = set(program_gene_union)
    valid_pool = np.ones(len(pool), dtype=bool) if finite_mask is None else np.asarray(finite_mask, dtype=bool)

    for gene in measured["gene_symbol"]:
        target = pool_idx.loc[gene]
        candidates = np.array([], dtype=int)
        used_expression_relax = 3
        used_lineage_relax = 3
        for lineage_relax in (0, 1, 2):
            for expression_relax in (0, 1, 2):
                mask = (
                    valid_pool
                    & (pool["biotype"].to_numpy() == target["biotype"])
                    & (pool["mt"].to_numpy() == target["mt"])
                    & (pool["ribo"].to_numpy() == target["ribo"])
                    & ((pool["mean_bin"] - target["mean_bin"]).abs().to_numpy() <= expression_relax)
                    & ((pool["detect_bin"] - target["detect_bin"]).abs().to_numpy() <= expression_relax)
                    & ((pool["lineage_corr_bin"] - target["lineage_corr_bin"]).abs().to_numpy() <= lineage_relax)
                    & (~pool["gene"].isin(excluded).to_numpy())
                )
                candidates = np.flatnonzero(mask)
                if len(candidates) >= 10:
                    used_expression_relax = expression_relax
                    used_lineage_relax = lineage_relax
                    break
            if len(candidates) >= 10:
                break
        if not len(candidates):
            mask = (
                valid_pool
                & (pool["biotype"].to_numpy() == target["biotype"])
                & ((pool["lineage_corr_bin"] - target["lineage_corr_bin"]).abs().to_numpy() <= 2)
                & (~pool["gene"].isin(excluded).to_numpy())
            )
            candidates = np.flatnonzero(mask)
            used_expression_relax = 3
            used_lineage_relax = 2
        if not len(candidates):
            mask = valid_pool & (pool["biotype"].to_numpy() == target["biotype"]) & (~pool["gene"].isin(excluded).to_numpy())
            candidates = np.flatnonzero(mask)
            used_expression_relax = 3
            used_lineage_relax = 3
        if not len(candidates):
            mask = valid_pool & (~pool["gene"].isin(excluded).to_numpy())
            candidates = np.flatnonzero(mask)
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
            pick = int(candidates[rng.integers(len(candidates))])
            if len(candidates) > len(chosen):
                while pick in chosen:
                    pick = int(candidates[rng.integers(len(candidates))])
            chosen.add(pick)
            sets[i, j] = pick
    return sets, expression_relaxation, lineage_relaxation


def standardized_against_null(observed: float, null: np.ndarray) -> tuple[float, float, float, float]:
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if len(null) < 100 or not np.isfinite(observed):
        return np.nan, np.nan, np.nan, np.nan
    mu = float(np.mean(null))
    sd = float(np.std(null, ddof=1))
    z = (observed - mu) / sd if sd > 0 else np.nan
    p = (1.0 + float(np.sum(null >= observed))) / (len(null) + 1.0)
    return z, p, mu, sd


def _score_weighted_mean_moran(
    bar_I_resid: np.ndarray,
    bar_I_raw: np.ndarray | None,
    bar_I_zon: np.ndarray | None,
    gene_indices: np.ndarray,
    weights: np.ndarray,
    sets: np.ndarray,
    n_sens: int,
) -> dict[str, object]:
    """Compute weighted single-gene mean Moran's I & matched-gene null stats in O(n_null * k)."""
    obs_res = float(np.dot(weights, bar_I_resid[gene_indices]))
    null_res = bar_I_resid[sets] @ weights
    res_z, res_p, res_mu, res_sd = standardized_against_null(obs_res, null_res)

    top = int(np.argmax(weights))
    eq_w = np.full(len(weights), 1.0 / len(weights), dtype=np.float64)
    lv_idx = np.delete(gene_indices, top)
    lv_w = np.delete(weights, top)
    lv_w = lv_w / lv_w.sum()

    sens = sets[: min(n_sens, len(sets))]
    obs_eq = float(np.mean(bar_I_resid[gene_indices]))
    null_eq = bar_I_resid[sens] @ eq_w
    obs_lv = float(np.dot(lv_w, bar_I_resid[lv_idx]))
    null_lv = bar_I_resid[np.delete(sens, top, axis=1)] @ lv_w

    out: dict[str, object] = {
        "residual_moran_i": obs_res,
        "residual_null_mean": res_mu,
        "residual_null_sd": res_sd,
        "residual_moran_z": res_z,
        "residual_pvalue": res_p,
        "equal_residual_moran_i": obs_eq,
        "equal_residual_null_mean": float(np.mean(null_eq)),
        "leave_top_residual_moran_i": obs_lv,
        "leave_top_residual_null_mean": float(np.mean(null_lv)),
        "sensitivity_sign_agree": bool(
            np.isfinite(res_z)
            and np.sign(res_z) == np.sign(obs_eq - float(np.mean(null_eq)))
            and np.sign(res_z) == np.sign(obs_lv - float(np.mean(null_lv)))
        ),
        "null_res": null_res,
    }
    if bar_I_raw is not None:
        obs_raw = float(np.dot(weights, bar_I_raw[gene_indices]))
        null_raw = bar_I_raw[sets] @ weights
        rz, rp, rmu, rsd = standardized_against_null(obs_raw, null_raw)
        out.update({"raw_moran_i": obs_raw, "raw_null_mean": rmu, "raw_null_sd": rsd, "raw_moran_z": rz, "raw_pvalue": rp})
    if bar_I_zon is not None:
        obs_zon = float(np.dot(weights, bar_I_zon[gene_indices]))
        null_zon = bar_I_zon[sets] @ weights
        zz, zp, zmu, zsd = standardized_against_null(obs_zon, null_zon)
        out.update({"zonation_moran_i": obs_zon, "zonation_null_mean": zmu, "zonation_null_sd": zsd, "zonation_moran_z": zz, "zonation_pvalue": zp})
    return out


def _evaluate_program_family(
    dataset: str,
    adata_path: Path,
    registry: pd.DataFrame,
    membership: pd.DataFrame,
    biotype: dict[str, str],
    family_tag: str,
    seed_offset: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    print(f"[spatial_20261118:{family_tag}] loading {dataset}: {adata_path}", flush=True)
    adata = ad.read_h5ad(adata_path)
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

    factor_names = [str(x).replace("means_per_cluster_mu_fg_", "") for x in adata.uns["mod"]["factor_names"]]
    q05_all = np.asarray(adata.obsm["q05_cell_abundance_w_sf"], dtype=float)
    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    graphs, graph_audit = section_graphs(obs, coords, k=6)
    graph_audit.insert(0, "dataset", dataset)

    g2i = {str(g): i for i, g in enumerate(adata.var_names)}
    pp_idx = [g2i[g] for g in PERIPORTAL if g in g2i]
    pc_idx = [g2i[g] for g in PERICENTRAL if g in g2i]
    zon_cov = None
    if pp_idx and pc_idx:
        pp_z = _standardize_matrix_cols(norm[:, pp_idx].toarray().astype(np.float64)).mean(axis=1)
        pc_z = _standardize_matrix_cols(norm[:, pc_idx].toarray().astype(np.float64)).mean(axis=1)
        zon_cov = pp_z - pc_z

    # Precompute within-section single-gene Moran matrices per required cell_type lineage
    lineages_needed = sorted(set(registry["cell_type"].astype(str)))
    lin_cache = {}
    for lin in lineages_needed:
        c2l_name = LINEAGE_TO_C2L[lin]
        cov_1d = q05_all[:, factor_names.index(c2l_name)]
        res_lin = compute_within_section_single_gene_moran(norm, obs, graphs, cov_1d, extra_cov=None, compute_raw=True)
        zon_lin = (
            compute_within_section_single_gene_moran(norm, obs, graphs, cov_1d, extra_cov=zon_cov, compute_raw=False)
            if zon_cov is not None else None
        )
        pool_lin = pool.copy()
        pool_lin["lineage_corr"] = abundance_correlation(norm, cov_1d)
        pool_lin["lineage_corr_bin"] = pd.qcut(pool_lin["lineage_corr"].rank(method="first"), 10, labels=False, duplicates="drop")
        finite_m = np.isfinite(res_lin["bar_I_resid"]) & np.isfinite(res_lin["bar_I_raw"])
        if zon_lin is not None:
            finite_m = finite_m & np.isfinite(zon_lin["bar_I_resid"])
        lin_cache[lin] = (res_lin, zon_lin, pool_lin, finite_m)

    union = set(membership["gene_symbol"].dropna().astype(str))
    rng = np.random.default_rng(SEED + seed_offset)

    res_rows = []
    sec_rows = []
    match_rows = []
    for row in registry.itertuples(index=False):
        pid = str(row.program_id)
        lin = str(row.cell_type)
        res_lin, zon_lin, pool_lin, finite_m = lin_cache[lin]
        bar_res = res_lin["bar_I_resid"]
        bar_raw = res_lin["bar_I_raw"]
        bar_zon = zon_lin["bar_I_resid"] if zon_lin is not None else None

        m = membership[membership["program_id"] == pid].groupby("gene_symbol", as_index=False)["original_l1_weight"].sum()
        m["present"] = m["gene_symbol"].isin(g2i)
        m["detected"] = m["gene_symbol"].map(lambda g: (detection[g2i[g]] >= 0.01 and bool(finite_m[g2i[g]])) if g in g2i else False)
        measured = m[m["present"] & m["detected"]].copy()
        retained = float(measured["original_l1_weight"].sum())
        testable = len(measured) >= 8 and retained >= 0.20

        base = {
            "program_id": pid,
            "dataset": dataset,
            "display_order": int(row.display_order) if hasattr(row, "display_order") else np.nan,
            "cell_type": lin,
            "module": int(row.module),
            "program_name": str(row.program_name) if hasattr(row, "program_name") else str(row.module_name),
            "n_measured": len(measured),
            "retained_l1_weight": retained,
            "testable": testable,
            "n_null": N_NULL,
            "statistic_definition": "within_section_nonlinear_ridge_weighted_single_gene_mean_moran_i",
        }
        if not testable:
            res_rows.append(base)
            continue

        measured["weight"] = measured["original_l1_weight"] / retained
        idx = np.array([g2i[g] for g in measured["gene_symbol"]], dtype=int)
        w = measured["weight"].to_numpy(np.float64)

        sets, er, lr = build_matched_sets(measured, pool_lin, union, N_NULL, rng, finite_mask=finite_m)
        match_rows.append({
            "dataset": dataset, "program_id": pid, "n_measured": len(measured),
            "n_genes_expression_relaxed": int(np.sum(np.asarray(er) > 0)),
            "n_genes_lineage_relaxed": int(np.sum(np.asarray(lr) > 0)),
            "matched_set_sha256": hashlib.sha256(sets.tobytes()).hexdigest(),
        })
        stats = _score_weighted_mean_moran(bar_res, bar_raw, bar_zon, idx, w, sets, N_SENSITIVITY_NULL)
        null_res = stats.pop("null_res")
        res_rows.append({**base, **stats, "matched_set_sha256": match_rows[-1]["matched_set_sha256"]})

        # Per-section breakdown
        ind_by_s = obs[["sample_id", "individual"]].drop_duplicates().set_index("sample_id")["individual"].astype(str).to_dict()
        for s in res_lin["arrays"]:
            sec_i = nan_weighted_dot(res_lin["sec_resid"][s][idx], w)
            sec_null = nan_weighted_mat(res_lin["sec_resid"][s][sets], w)
            sec_null_mean = float(np.nanmean(sec_null))
            sec_rows.append({
                "dataset": dataset, "program_id": pid, "sample_id": s, "individual": ind_by_s[s],
                "residual_moran_i": sec_i, "null_mean": sec_null_mean,
                "centered_excess": sec_i - sec_null_mean,
            })

    res_df = pd.DataFrame(res_rows)
    res_df["residual_padj"] = bh(res_df["residual_pvalue"])
    if "raw_pvalue" in res_df.columns:
        res_df["raw_padj"] = bh(res_df["raw_pvalue"])
    if "zonation_pvalue" in res_df.columns:
        res_df["zonation_padj"] = bh(res_df["zonation_pvalue"])
    res_df["robust"] = (
        res_df["testable"].astype("boolean").fillna(False).to_numpy(bool)
        & (res_df["residual_padj"] < 0.05).fillna(False).to_numpy(bool)
        & res_df["sensitivity_sign_agree"].astype("boolean").fillna(False).to_numpy(bool)
    )

    del adata, counts, norm
    gc.collect()
    return res_df, pd.DataFrame(sec_rows), pd.DataFrame(match_rows), graph_audit


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gene_meta = pd.read_csv(GENE_META_FILE, sep="\t", usecols=["gene_name", "gene_biotype", "chromosome"], compression="gzip")
    primary_chr = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M"]}
    gene_meta_prim = gene_meta[gene_meta["chromosome"].isin(primary_chr)].drop_duplicates("gene_name")
    biotype = dict(zip(gene_meta_prim["gene_name"].astype(str), gene_meta_prim["gene_biotype"].astype(str)))

    # 1. V1 22-program family
    v1_reg = pd.read_csv(V1_REGISTRY_FILE, sep="\t")
    v1_mem = pd.read_csv(V1_MEMBERSHIP_FILE, sep="\t")
    v1_mem = v1_mem[v1_mem["mapped_symbol"].fillna(False) & v1_mem["gene_symbol"].notna()].copy()

    all_r1, all_s1, all_m1, all_g = [], [], [], []
    for ds_idx, (ds, path) in enumerate(DATASETS.items()):
        r, s, m, g = _evaluate_program_family(ds, path, v1_reg, v1_mem, biotype, "v1_22", 1000 * ds_idx)
        all_r1.append(r); all_s1.append(s); all_m1.append(m); all_g.append(g)
    pd.concat(all_r1, ignore_index=True).to_csv(OUT / "spatial_program_results.tsv", sep="\t", index=False)
    pd.concat(all_s1, ignore_index=True).to_csv(OUT / "spatial_section_results.tsv", sep="\t", index=False)
    pd.concat(all_m1, ignore_index=True).to_csv(OUT / "spatial_matching_audit.tsv", sep="\t", index=False)
    pd.concat(all_g, ignore_index=True).to_csv(OUT / "spatial_graph_audit.tsv", sep="\t", index=False)

    # 2. V2 117-program family (GENCODE v49 primary-assembly deterministic mapping)
    v2_reg = pd.read_csv(V2_REGISTRY_FILE, sep="\t")
    v2_reg["program_id"] = v2_reg["program_uid"]
    v2_reg["program_name"] = v2_reg["module_name"]
    v2_reg = v2_reg.sort_values(["cell_type", "module"]).reset_index(drop=True)
    v2_reg["display_order"] = np.arange(1, len(v2_reg) + 1)

    v2_mem = pd.read_csv(V2_MEMBERSHIP_FILE, sep="\t")
    v2_mem["program_id"] = v2_mem["program_uid"]
    prim_counts = gene_meta_prim.groupby("gene_name").size().to_dict()
    sym = v2_mem["mapped_symbol"].fillna("").astype(str).str.strip()
    canon = v2_mem["canonical_gene"].fillna("").astype(str).str.strip()
    v2_mem["gene_symbol"] = np.where(sym != "", sym, canon)
    in_prim = (
        (v2_mem["mapped_symbol_status"] == "gencode_v49_unique_symbol_confirmed")
        | (v2_mem["mapped_symbol_status"] == "gencode_v49_unambiguous_ensembl_to_symbol")
        | ((v2_mem["mapped_symbol_status"] == "unmapped_or_ambiguous") & (v2_mem["gene_symbol"].map(lambda g: prim_counts.get(g, 0) == 1)))
    )
    v2_mem = v2_mem[in_prim & (v2_mem["gene_symbol"] != "")].copy()

    all_r2, all_s2, all_m2 = [], [], []
    for ds_idx, (ds, path) in enumerate(DATASETS.items()):
        r, s, m, _ = _evaluate_program_family(ds, path, v2_reg, v2_mem, biotype, "v2_117", 1000 * ds_idx)
        all_r2.append(r); all_s2.append(s); all_m2.append(m)
    res117 = pd.concat(all_r2, ignore_index=True)
    res117.to_csv(OUT / "spatial_program_results_117.tsv", sep="\t", index=False)
    pd.concat(all_s2, ignore_index=True).to_csv(OUT / "spatial_section_results_117.tsv", sep="\t", index=False)
    pd.concat(all_m2, ignore_index=True).to_csv(OUT / "spatial_matching_audit_117.tsv", sep="\t", index=False)

    print("[spatial_20261118] V1 22-program testable/robust:")
    print(pd.concat(all_r1, ignore_index=True).groupby("dataset")[["testable", "robust"]].sum())
    print("[spatial_20261118] V2 117-program testable/robust:")
    print(res117.groupby("dataset")[["testable", "robust"]].sum())


if __name__ == "__main__":
    main()

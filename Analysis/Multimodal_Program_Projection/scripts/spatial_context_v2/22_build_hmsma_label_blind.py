#!/usr/bin/env python3
"""Build label-blind HMSMA two-program descriptive spatial organization.

This branch deliberately produces no population P value. Physical arrays are
unresolved technical units until an authoritative array-to-donor key exists.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from spatial_resource_lib import (
    EXPECTED_REGISTRY_SHA256,
    PROGRAM_RELEASE_ID,
    RESOURCE_RELEASE_ID,
    SpatialResourceError,
    load_frozen_programs,
    sha256_file,
    write_tsv,
)


OFFSETS = ((0, -2), (0, 2), (-1, -1), (-1, 1), (1, -1), (1, 1))


def moran_hex(values: np.ndarray, rows: np.ndarray, columns: np.ndarray) -> tuple[float, int, int]:
    """Binary six-neighbor Moran's I on the fixed Visium hex lattice."""
    positions = {(int(row), int(column)): index for index, (row, column) in enumerate(zip(rows, columns))}
    edges = []
    degree = np.zeros(len(values), dtype=int)
    for index, (row, column) in enumerate(zip(rows, columns)):
        for drow, dcolumn in OFFSETS:
            neighbor = positions.get((int(row + drow), int(column + dcolumn)))
            if neighbor is not None:
                edges.append((index, neighbor))
                degree[index] += 1
    eligible = degree > 0
    if eligible.sum() < 3 or not edges:
        return float("nan"), int(eligible.sum()), 0
    retained_indices = np.flatnonzero(eligible)
    retained_values = values[retained_indices]
    centered = retained_values - np.mean(retained_values)
    old_to_new = {old: new for new, old in enumerate(retained_indices)}
    retained_edges = [(old_to_new[left], old_to_new[right]) for left, right in edges if eligible[left] and eligible[right]]
    denominator = float(np.dot(centered, centered))
    if denominator <= 0:
        return float("nan"), int(eligible.sum()), len(edges)
    numerator = sum(centered[left] * centered[right] for left, right in retained_edges)
    statistic = len(retained_values) / len(retained_edges) * numerator / denominator
    return float(statistic), int(eligible.sum()), len(retained_edges)


def gene_symbols(adata: ad.AnnData) -> list[str]:
    if "gene_symbol" in adata.var:
        return adata.var["gene_symbol"].astype(str).tolist()
    return [str(value).split(".")[0] for value in adata.var_names]


def score_program(matrix, total_counts: np.ndarray, detected_genes: np.ndarray, symbols: list[str], weights: dict[str, float]):
    symbol_to_columns: dict[str, list[int]] = {}
    for index, symbol in enumerate(symbols):
        symbol_to_columns.setdefault(symbol, []).append(index)
    measured = []
    for gene in sorted(set(weights) & set(symbol_to_columns)):
        gene_matrix = matrix[:, symbol_to_columns[gene]]
        if sp.issparse(gene_matrix):
            detected_fraction = float(np.asarray((gene_matrix > 0).sum(axis=1)).ravel().astype(bool).mean())
        else:
            detected_fraction = float(np.any(gene_matrix > 0, axis=1).mean())
        if detected_fraction >= 0.01:
            measured.append(gene)
    if not measured:
        raise SpatialResourceError("HMSMA program has no genes detected in at least 1% of array spots")
    columns = [column for gene in measured for column in symbol_to_columns[gene]]
    column_weights = np.asarray([weights[gene] / len(symbol_to_columns[gene]) for gene in measured for _ in symbol_to_columns[gene]])
    normalized = matrix[:, columns].astype(np.float64)
    scale = np.divide(10_000.0, total_counts, out=np.zeros_like(total_counts, dtype=float), where=total_counts > 0)
    normalized = sp.diags(scale) @ normalized if sp.issparse(normalized) else normalized * scale[:, None]
    if sp.issparse(normalized):
        normalized.data = np.log1p(normalized.data)
        score = np.asarray(normalized @ column_weights).ravel()
    else:
        score = np.log1p(normalized) @ column_weights
    design = np.column_stack([
        np.ones(len(score)),
        (np.log1p(total_counts) - np.mean(np.log1p(total_counts))) / max(np.std(np.log1p(total_counts)), 1e-12),
        (detected_genes - np.mean(detected_genes)) / max(np.std(detected_genes), 1e-12),
    ])
    fitted = design @ np.linalg.lstsq(design, score, rcond=None)[0]
    residual = score - fitted
    residual = (residual - np.mean(residual)) / max(np.std(residual), 1e-12)
    return residual, measured, float(sum(weights[gene] for gene in measured))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--candidate-root", type=Path, required=True)
    args = parser.parse_args()
    project = args.project_root.resolve()
    candidate = args.candidate_root.resolve()
    input_root = project / "Analysis/Spatial/results/preprocessed/HRA007511_starsolo"
    output = candidate / "hmsma_label_blind"
    if (output / "READY").exists():
        raise SpatialResourceError(f"immutable HMSMA candidate already exists: {output}")
    summary = json.loads((input_root / "build_summary.json").read_text(encoding="utf-8"))
    if summary != {
        "n_samples_built": 35,
        "n_samples_missing": 0,
        "samples_missing": [],
        "n_spots_total": 130097,
        "n_genes": 86369,
        "min_umi": 500,
        "phenotype_key_available": False,
    }:
        raise SpatialResourceError(f"HMSMA build summary drift or phenotype gate changed: {summary}")

    hotspot = project / "Analysis/Multimodal_Program_Projection/candidates" / PROGRAM_RELEASE_ID / "hotspot"
    programs, weights = load_frozen_programs(hotspot)
    selected = {uid: program for uid, program in programs.items() if program["external_test_eligible"].upper() == "TRUE"}
    files = sorted(input_root.glob("HRA_*.h5ad"))
    if len(files) != 35:
        raise SpatialResourceError(f"expected 35 HMSMA array files, found {len(files)}")

    array_rows = []
    map_frames = []
    spot_counts = []
    cached_maps = {}
    for path in files:
        adata = ad.read_h5ad(path)
        if "array_row" not in adata.obs or "array_col" not in adata.obs:
            raise SpatialResourceError(f"HMSMA array lacks registered array coordinates: {path.name}")
        total = np.asarray(adata.X.sum(axis=1)).ravel()
        detected = np.asarray((adata.X > 0).sum(axis=1)).ravel()
        symbols = gene_symbols(adata)
        spot_counts.append((path.stem, adata.n_obs))
        cached_maps[path.stem] = []
        for uid, program in selected.items():
            residual, measured, retained = score_program(adata.X, total, detected, symbols, weights[uid])
            statistic, n_graph, n_directed_edges = moran_hex(
                residual,
                adata.obs["array_row"].to_numpy(),
                adata.obs["array_col"].to_numpy(),
            )
            array_rows.append({
                "release_id": RESOURCE_RELEASE_ID,
                "program_release_id": PROGRAM_RELEASE_ID,
                "registry_sha256": EXPECTED_REGISTRY_SHA256,
                "dataset_id": "HRA007511_HMSMA",
                "assay_id": "visium_rna",
                "array_id": path.stem,
                "program_uid": uid,
                "program_label": program["module_name"],
                "membership_sha256": program["membership_sha256"],
                "n_spots": adata.n_obs,
                "n_graph_eligible_spots": n_graph,
                "n_directed_edges": n_directed_edges,
                "n_genes_measured": len(measured),
                "retained_l1_weight": f"{retained:.17g}",
                "residual_moran_i": f"{statistic:.17g}",
                "score_unit": "frozen_positive_weight_log1p_CPT_residual_z",
                "adjustment": "within_array_log_library_size_and_detected_gene_count",
                "inferential_pvalue_authorized": "FALSE",
            })
            cached_maps[path.stem].append((uid, program, residual, adata.obs.copy(), len(measured), retained))

    ordered_counts = sorted(count for _, count in spot_counts)
    median_count = float(np.median(ordered_counts))
    selected_array = min(spot_counts, key=lambda item: (abs(item[1] - median_count), item[0]))[0]
    for uid, program, residual, obs, n_measured, retained in cached_maps[selected_array]:
        map_frames.append(pd.DataFrame({
            "release_id": RESOURCE_RELEASE_ID,
            "dataset_id": "HRA007511_HMSMA",
            "array_id": selected_array,
            "selection_rule": "spot_count_closest_to_dataset_median_then_lexical_tie_break",
            "program_uid": uid,
            "program_label": program["module_name"],
            "spot_id": obs.index.astype(str),
            "array_row": obs["array_row"].to_numpy(),
            "array_col": obs["array_col"].to_numpy(),
            "residual_program_score_z": residual,
            "n_genes_measured": n_measured,
            "retained_l1_weight": retained,
        }))

    output.mkdir(parents=True, exist_ok=True)
    array_columns = tuple(array_rows[0])
    write_tsv(output / "per_array_program_organization.tsv", array_columns, array_rows)
    pd.concat(map_frames, ignore_index=True).to_parquet(output / "selected_array_program_maps.parquet", index=False)
    write_tsv(
        output / "metadata_gate.tsv",
        ("required_field", "status", "claim_if_missing"),
        [
            {"required_field": field, "status": "unresolved", "claim_if_missing": "prohibited"}
            for field in (
                "array_to_donor_or_clinical_id", "repeat_section_relationship", "technical_run_aggregation",
                "control_MASLD_MASH_label", "NAS", "fibrosis", "age_sex_BMI_batch",
                "array_to_HE_registration", "spatial_to_scRNA_donor_overlap", "MALDI_registration",
            )
        ],
    )
    write_tsv(
        output / "READY",
        ("release_id", "status", "n_arrays", "n_spots_technical", "n_programs", "population_pvalues", "selected_map_array", "per_array_sha256", "map_sha256"),
        [{
            "release_id": RESOURCE_RELEASE_ID,
            "status": "metadata_pending_label_blind_descriptive_only",
            "n_arrays": 35,
            "n_spots_technical": 130097,
            "n_programs": 2,
            "population_pvalues": "prohibited",
            "selected_map_array": selected_array,
            "per_array_sha256": sha256_file(output / "per_array_program_organization.tsv"),
            "map_sha256": sha256_file(output / "selected_array_program_maps.parquet"),
        }],
    )
    print(json.dumps({"arrays": 35, "technical_spots": 130097, "selected_map_array": selected_array}, indent=2))


if __name__ == "__main__":
    main()

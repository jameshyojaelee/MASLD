"""Geometry-preserving continual mapping into a frozen reference latent."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, DeterministicGzipTextWriter, sha256_path
from .control_adapter import _donor_balanced_centroid
from .embedding import load_embedding, matched_rows
from .execution import require_execution_ownership, source_identity, verify_execution_record
from .firewall import validate_program_firewall


def fit_scaled_orthogonal_bridge(
    source: np.ndarray, target: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    source = np.asarray(source, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    if source.shape != target.shape or source.ndim != 2 or source.shape[1] < 2:
        raise ContractError("orthogonal bridge requires matched two-dimensional matrices")
    if not np.isfinite(source).all() or not np.isfinite(target).all():
        raise ContractError("orthogonal bridge received non-finite values")
    source_mean = source.mean(axis=0)
    target_mean = target.mean(axis=0)
    centered_source = source - source_mean
    centered_target = target - target_mean
    u, singular, vt = np.linalg.svd(centered_source.T @ centered_target, full_matrices=False)
    rotation = u @ vt
    denominator = float(np.square(centered_source).sum())
    scale = float(singular.sum() / denominator) if denominator > 0 else 0.0
    if scale <= 0 or not np.isfinite(scale):
        raise ContractError("orthogonal bridge scale is not positive and finite")
    identity_error = float(np.max(np.abs(rotation.T @ rotation - np.eye(rotation.shape[1]))))
    if identity_error > 1e-10:
        raise ContractError("orthogonal bridge rotation is not orthonormal")
    return source_mean, target_mean, rotation, scale


def apply_scaled_orthogonal_bridge(
    values: np.ndarray, source_mean: np.ndarray, target_mean: np.ndarray,
    rotation: np.ndarray, scale: float,
) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    return ((values - source_mean) @ rotation * scale + target_mean).astype(np.float32)


def _balanced_fit_positions(
    cells, positions: np.ndarray, model_kind: str, config: dict[str, Any], seed: int | None = None,
):
    selected = cells.iloc[positions]
    if model_kind == "all_lineage":
        keys = selected["donor_id"].astype(str) + "|" + selected["audit_cell_type"].astype(str)
        cap = int(config["sampling"]["all_lineage_cap_per_donor_label"])
    else:
        keys = selected["donor_id"].astype(str)
        cap = int(config["sampling"]["lineage_cap_per_donor"])
    rng = np.random.default_rng(int(config["screen"]["seed"] if seed is None else seed))
    keep = []
    key_values = keys.to_numpy()
    for key in sorted(set(key_values)):
        local = np.flatnonzero(key_values == key)
        if len(local) > cap:
            local = np.sort(rng.choice(local, size=cap, replace=False))
        keep.extend(positions[local])
    result = np.asarray(sorted(keep), dtype=np.int64)
    if not len(result):
        raise ContractError("orthogonal bridge fit roster is empty")
    return result


def apply_dataset_control_offsets(
    latent: np.ndarray, cells, query_datasets: list[str], weight: float,
    minimum_control_donors: int,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    if not np.isfinite(weight) or not 0 <= weight <= 1:
        raise ContractError("orthogonal bridge offset weight is invalid")
    reference = cells["strict_reference"].to_numpy(dtype=bool)
    primary = cells["primary_query"].to_numpy(dtype=bool)
    controls = primary & cells["query_control"].to_numpy(dtype=bool) & ~reference
    reference_centroid, n_reference = _donor_balanced_centroid(latent, cells, reference)
    result = np.asarray(latent).copy()
    summaries = []
    datasets = cells["dataset"].astype(str).to_numpy()
    for dataset in query_datasets:
        control_mask = controls & (datasets == dataset)
        control_centroid, n_control = _donor_balanced_centroid(latent, cells, control_mask)
        if n_control < minimum_control_donors:
            raise ContractError(f"underpowered bridge control group: {dataset}")
        full_offset = reference_centroid - control_centroid
        offset = weight * full_offset
        apply_mask = primary & (datasets == dataset)
        before = np.asarray(latent[apply_mask], dtype=np.float64)
        result[apply_mask] = (before + offset).astype(result.dtype, copy=False)
        after = np.asarray(result[apply_mask], dtype=np.float64)
        centered_error = float(np.max(np.abs(
            (after - after.mean(axis=0)) - (before - before.mean(axis=0))
        ), initial=0.0))
        summaries.append({
            "dataset": dataset, "preparations_tied": sorted(set(cells.loc[apply_mask, "preparation"].astype(str))),
            "n_control_donors": n_control, "n_reference_donors": n_reference,
            "n_applied_cells": int(apply_mask.sum()), "control_offset_weight": weight,
            "full_offset": list(map(float, full_offset)), "applied_offset": list(map(float, offset)),
            "maximum_centered_coordinate_error": centered_error,
        })
    if not np.array_equal(result[reference], latent[reference]):
        raise ContractError("orthogonal bridge calibration changed frozen reference coordinates")
    return result, summaries


def build_orthogonal_bridge_grid(
    config: dict[str, Any], policy_value: str | Path,
    raw_pca_embedding_value: str | Path,
    architecture_embedding_value: str | Path, architecture_execution_record_value: str | Path,
    reference_embedding_value: str | Path, reference_execution_record_value: str | Path,
    output_root_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    pipeline_root = Path(config["_config_path"]).resolve().parent
    policy_path = Path(policy_value).resolve()
    if policy_path != pipeline_root / "reference" / "orthogonal_bridge_policy_v8.json":
        raise ContractError("orthogonal bridge requires the source-controlled v8 policy")
    with policy_path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != "masld-cl-orthogonal-bridge-policy-v8"
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("case_stage_program_hero_gene_and_cas13_outcomes_locked") is not True
    ):
        raise ContractError("orthogonal bridge policy is invalid")
    for key in ("parent_rejection", "comparator_diagnostic"):
        source = pipeline_root / policy[key]["relative_path"]
        if sha256_path(source) != policy[key]["sha256"]:
            raise ContractError(f"orthogonal bridge policy source changed: {key}")

    raw_path = Path(raw_pca_embedding_value).resolve()
    raw_info, raw_latent, raw_cells = load_embedding(raw_path)
    if raw_info.get("method") != "uncorrected_scalesc_pca_30d":
        raise ContractError("orthogonal bridge raw representation is invalid")
    architecture_path = Path(architecture_embedding_value).resolve()
    architecture_info, _, cells = load_embedding(architecture_path)
    architecture_run_path = architecture_path.parent / "update_manifest.json"
    with architecture_run_path.open() as handle:
        architecture_run = json.load(handle)
    model_kind = architecture_run.get("model_kind")
    if (
        architecture_run.get("method") != "architecture_surgery"
        or architecture_run.get("embedding") != architecture_info
        or model_kind not in {"all_lineage", *config["lineages"]}
        or set(architecture_run.get("query_datasets", [])) != set(policy["primary_query_datasets"])
    ):
        raise ContractError("orthogonal bridge architecture source is invalid")
    architecture_record = verify_execution_record(
        architecture_execution_record_value, pipeline_root, config["_config_sha256"],
        allow_historical_source=True,
    )
    require_execution_ownership(
        architecture_record, [architecture_path, architecture_run_path],
        role="orthogonal bridge architecture source",
    )
    reference_path = Path(reference_embedding_value).resolve()
    reference_info, reference_latent, reference_cells = load_embedding(reference_path)
    reference_run_path = reference_path.parent / "reference_manifest.json"
    with reference_run_path.open() as handle:
        reference_run = json.load(handle)
    if (
        reference_run.get("schema_version") != "masld-cl-reference-v1"
        or reference_run.get("model_kind") != model_kind
        or reference_run.get("embedding") != reference_info
    ):
        raise ContractError("orthogonal bridge reference source is invalid")
    reference_record = verify_execution_record(
        reference_execution_record_value, pipeline_root, config["_config_sha256"],
        allow_historical_source=True,
    )
    require_execution_ownership(
        reference_record, [reference_path, reference_run_path], role="orthogonal bridge reference source"
    )

    raw_left, candidate_right = matched_rows(raw_cells, cells)
    if len(candidate_right) != len(cells) or not np.array_equal(candidate_right, np.arange(len(cells))):
        order = np.argsort(candidate_right)
        raw_left, candidate_right = raw_left[order], candidate_right[order]
    if not np.array_equal(candidate_right, np.arange(len(cells))):
        raise ContractError("raw PCA does not contain the exact architecture cell roster")
    raw_values = np.asarray(raw_latent[raw_left], dtype=np.float32)
    reference_positions = np.flatnonzero(cells["strict_reference"].to_numpy(dtype=bool))
    ref_left, candidate_ref_right = matched_rows(
        reference_cells, cells.iloc[reference_positions].reset_index(drop=True)
    )
    candidate_reference_positions = reference_positions[candidate_ref_right]
    fit_positions = _balanced_fit_positions(cells, candidate_reference_positions, model_kind, config)
    ref_index_by_cell = {
        value: index for index, value in enumerate(reference_cells["cell_id"].astype(str))
    }
    fit_ref_positions = np.asarray(
        [ref_index_by_cell[value] for value in cells.iloc[fit_positions]["cell_id"].astype(str)],
        dtype=np.int64,
    )
    source_mean, target_mean, rotation, scale = fit_scaled_orthogonal_bridge(
        raw_values[fit_positions], np.asarray(reference_latent)[fit_ref_positions]
    )
    bridged = apply_scaled_orthogonal_bridge(
        raw_values, source_mean, target_mean, rotation, scale
    )
    bridged[candidate_reference_positions] = np.asarray(reference_latent)[ref_left].astype(np.float32)
    cells = cells.copy()
    prediction_column = cells.columns.get_loc("predicted_cell_type")
    cells.iloc[candidate_reference_positions, prediction_column] = (
        reference_cells.iloc[ref_left]["predicted_cell_type"].astype(str).to_numpy()
    )

    output_root = Path(output_root_value)
    output_root.mkdir(parents=True, exist_ok=False)
    outputs = []
    frozen_source = source_identity(pipeline_root)
    for weight in map(float, policy["control_offset_weight_grid"]):
        adapted, offsets = apply_dataset_control_offsets(
            bridged, cells, policy["primary_query_datasets"], weight,
            int(policy["minimum_control_donors_per_dataset"]),
        )
        output = output_root / f"w{int(round(weight * 100)):03d}"
        output.mkdir()
        latent_path = output / "embedding_latent.npy"
        cells_path = output / architecture_info["cells_file"]
        np.save(latent_path, adapted)
        with DeterministicGzipTextWriter(cells_path) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(cells.columns)
            writer.writerows(cells.itertuples(index=False, name=None))
        embedding = dict(architecture_info)
        embedding.update({
            "method": "orthogonal_bridge_continual_adapter",
            "n_latent": int(adapted.shape[1]), "latent_file": latent_path.name,
            "cells_file": cells_path.name, "latent_sha256": sha256_path(latent_path),
            "cells_sha256": sha256_path(cells_path),
        })
        write_json_exclusive(output / "embedding_manifest.json", embedding)
        manifest = {
            "schema_version": "masld-cl-orthogonal-bridge-v8",
            "config_sha256": config["_config_sha256"], "policy_realpath": str(policy_path),
            "policy_sha256": sha256_path(policy_path), "source_identity": frozen_source,
            "model_kind": model_kind, "control_offset_weight": weight,
            "control_only": True, "outcomes_unlocked": False,
            "reference_coordinates_and_predictions_frozen_to_preupdate": True,
            "query_predictions_frozen_to_architecture_surgery": True,
            "raw_pca_embedding_realpath": str(raw_path),
            "raw_pca_embedding_sha256": sha256_path(raw_path),
            "architecture_embedding_realpath": str(architecture_path),
            "architecture_embedding_sha256": sha256_path(architecture_path),
            "architecture_run_realpath": str(architecture_run_path),
            "architecture_run_sha256": sha256_path(architecture_run_path),
            "architecture_execution_record_realpath": str(Path(architecture_execution_record_value).resolve()),
            "architecture_execution_record_sha256": sha256_path(architecture_execution_record_value),
            "reference_embedding_realpath": str(reference_path),
            "reference_embedding_sha256": sha256_path(reference_path),
            "reference_run_realpath": str(reference_run_path),
            "reference_run_sha256": sha256_path(reference_run_path),
            "reference_execution_record_realpath": str(Path(reference_execution_record_value).resolve()),
            "reference_execution_record_sha256": sha256_path(reference_execution_record_value),
            "bridge": {
                "fit_cells": int(len(fit_positions)), "scale": scale,
                "source_mean": source_mean.tolist(), "target_mean": target_mean.tolist(),
                "rotation": rotation.tolist(),
                "orthonormality_max_abs_error": float(np.max(np.abs(rotation.T @ rotation - np.eye(rotation.shape[1])))),
            },
            "offsets": offsets, "embedding": embedding,
        }
        write_json_exclusive(output / "orthogonal_bridge_manifest.json", manifest)
        outputs.append({
            "control_offset_weight": weight,
            "embedding_manifest": str((output / "embedding_manifest.json").resolve()),
            "bridge_manifest": str((output / "orthogonal_bridge_manifest.json").resolve()),
            "bridge_manifest_sha256": sha256_path(output / "orthogonal_bridge_manifest.json"),
        })
    grid = {
        "schema_version": "masld-cl-orthogonal-bridge-grid-v8",
        "config_sha256": config["_config_sha256"], "model_kind": model_kind,
        "policy_sha256": sha256_path(policy_path), "source_identity": frozen_source,
        "outputs": outputs,
    }
    write_json_exclusive(output_root / "grid_manifest.json", grid)
    return grid

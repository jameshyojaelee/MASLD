"""Outcome-blind rigid alignment and exact strict7 restoration for V22."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, DeterministicGzipTextWriter, sha256_path
from .embedding import load_embedding, matched_rows
from .firewall import validate_program_firewall
from .sampling import capped_group_indices


class ExternalHealthyAnchorError(ContractError):
    """Raised when the outcome-blind V23 export changes identity."""


POLICY_SCHEMA = "masld-cl-external-healthy-anchor-policy-v23"


def load_external_healthy_anchor_policy(
    config: dict[str, Any], value: str | Path,
) -> tuple[Path, dict[str, Any], dict[str, Any]]:
    path = Path(value).resolve()
    expected = (
        Path(config["_config_path"]).resolve().parent
        / "reference" / "external_healthy_anchor_policy_v23.json"
    )
    if path != expected:
        raise ExternalHealthyAnchorError("V23 requires its source-controlled policy")
    with path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != POLICY_SCHEMA
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("status") != "single_outcome_blind_stop_go_export"
        or policy.get("case_stage_program_hero_gene_umap_cas13_used") is not False
        or policy.get("evaluation", {}).get("reuse_v22_thresholds_exactly") is not True
        or policy.get("evaluation", {}).get("no_additional_hyperparameter_or_candidate_grid") is not True
    ):
        raise ExternalHealthyAnchorError("V23 policy identity differs")
    parent_spec = policy["parent_decision"]
    parent_path = (path.parent / parent_spec["path"]).resolve()
    with parent_path.open() as handle:
        parent = json.load(handle)
    failed = [row["gate"] for row in parent.get("gates", []) if not row.get("pass")]
    if (
        sha256_path(parent_path) != parent_spec["sha256"]
        or parent.get("decision") != parent_spec["required_decision"]
        or failed != parent_spec["required_failed_gates"]
        or parent.get("case_stage_program_hero_gene_umap_cas13_read") is not False
    ):
        raise ExternalHealthyAnchorError("V23 parent decision differs")
    v22_spec = policy["parent_v22_policy"]
    v22_path = (path.parent / v22_spec["path"]).resolve()
    if sha256_path(v22_path) != v22_spec["sha256"]:
        raise ExternalHealthyAnchorError("V23 parent V22 policy changed")
    return path, policy, parent


def _policy_embedding_path(policy_path: Path, spec: dict[str, Any]) -> Path:
    path = (policy_path.parent / spec["path"]).resolve()
    if sha256_path(path) != spec["sha256"]:
        raise ExternalHealthyAnchorError("V23 source embedding changed")
    return path


def _balanced_anchor_positions(cells, maximum: int, seed: int) -> np.ndarray:
    reference_positions = np.flatnonzero(
        cells["strict_reference"].to_numpy(dtype=bool)
    )
    reference = cells.iloc[reference_positions]
    groups = list(zip(
        reference["donor_id"].astype(str),
        reference["audit_cell_type"].astype(str),
    ))
    cap = max(1, int(np.ceil(maximum / len(set(groups)))))
    local = capped_group_indices(groups, cap, seed)
    positions = reference_positions[local]
    if len(positions) > maximum:
        positions = np.random.default_rng(seed + 1).choice(
            positions, maximum, replace=False
        )
    return np.sort(positions.astype(np.int64))


def _proper_rigid_alignment(source: np.ndarray, target: np.ndarray):
    source = np.asarray(source, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    if source.shape != target.shape or source.ndim != 2 or len(source) < source.shape[1]:
        raise ExternalHealthyAnchorError("V23 rigid-alignment anchors differ")
    source_mean = source.mean(axis=0)
    target_mean = target.mean(axis=0)
    u, singular, vt = np.linalg.svd(
        (source - source_mean).T @ (target - target_mean), full_matrices=False
    )
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    orthogonality_error = float(np.max(np.abs(
        rotation.T @ rotation - np.eye(rotation.shape[0])
    )))
    if np.linalg.det(rotation) <= 0 or orthogonality_error > 1e-10:
        raise ExternalHealthyAnchorError("V23 transform is not a proper rotation")
    return source_mean, target_mean, rotation, singular, orthogonality_error


def build_external_healthy_anchor_export(
    config: dict[str, Any], policy_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    policy_path, policy, parent = load_external_healthy_anchor_policy(
        config, policy_value
    )
    candidate_spec = policy["candidate_embedding"]
    common_spec = policy["strict7_embedding"]
    candidate_path = _policy_embedding_path(policy_path, candidate_spec)
    common_path = _policy_embedding_path(policy_path, common_spec)
    candidate_manifest_path = candidate_path.parent / "continual_reference_manifest.json"
    common_manifest_path = common_path.parent / "reference_manifest.json"
    if (
        sha256_path(candidate_manifest_path) != candidate_spec["training_manifest_sha256"]
        or sha256_path(common_manifest_path) != common_spec["training_manifest_sha256"]
    ):
        raise ExternalHealthyAnchorError("V23 source training manifest changed")
    candidate = load_embedding(candidate_path)
    common = load_embedding(common_path)
    left, right = matched_rows(common[2], candidate[2])
    if len(left) != len(common[2]) or len(right) != len(candidate[2]):
        raise ExternalHealthyAnchorError("V23 source cell rosters differ")
    cells = candidate[2].copy().reset_index(drop=True)
    common_latent = np.empty_like(np.asarray(candidate[1]))
    common_latent[right] = np.asarray(common[1])[left]
    common_cells = common[2].iloc[left].reset_index(drop=True)
    candidate_cells_sorted = cells.iloc[right].reset_index(drop=True)
    if not np.array_equal(
        common_cells["audit_cell_type"].astype(str).to_numpy(),
        candidate_cells_sorted["audit_cell_type"].astype(str).to_numpy(),
    ):
        raise ExternalHealthyAnchorError("V23 frozen audit labels differ")
    common_predictions = np.empty(len(cells), dtype=object)
    common_predictions[right] = common_cells["predicted_cell_type"].astype(str).to_numpy()
    reference = cells["strict_reference"].to_numpy(dtype=bool)
    anchor_policy = policy["anchor"]
    anchor_positions = _balanced_anchor_positions(
        cells, anchor_policy["maximum_cells"], anchor_policy["seed"]
    )
    source_mean, target_mean, rotation, singular, ortho_error = _proper_rigid_alignment(
        np.asarray(candidate[1])[anchor_positions], common_latent[anchor_positions]
    )
    transformed = (
        (np.asarray(candidate[1], dtype=np.float64) - source_mean) @ rotation
        + target_mean
    ).astype(np.float32)
    transformed[reference] = common_latent[reference]
    prediction_column = cells.columns.get_loc("predicted_cell_type")
    query_predictions = cells.loc[~reference, "predicted_cell_type"].astype(str).to_numpy()
    cells.iloc[np.flatnonzero(reference), prediction_column] = common_predictions[reference]
    if not np.array_equal(
        cells.loc[~reference, "predicted_cell_type"].astype(str).to_numpy(),
        query_predictions,
    ):
        raise ExternalHealthyAnchorError("V23 changed query predictions")
    if not np.array_equal(transformed[reference], common_latent[reference]):
        raise ExternalHealthyAnchorError("V23 did not restore strict7 coordinates exactly")
    query_positions = np.flatnonzero(~reference)
    audit_positions = np.random.default_rng(anchor_policy["seed"] + 2).choice(
        query_positions, min(1000, len(query_positions)), replace=False
    )
    before = np.asarray(candidate[1])[audit_positions]
    after = transformed[audit_positions]
    before_dist = np.linalg.norm(np.diff(before.astype(np.float64), axis=0), axis=1)
    after_dist = np.linalg.norm(np.diff(after.astype(np.float64), axis=0), axis=1)
    distance_error = float(np.max(np.abs(before_dist - after_dist), initial=0.0))
    if distance_error > 1e-4:
        raise ExternalHealthyAnchorError("V23 rigid transform changed query distances")

    output = Path(output_value)
    output.mkdir(parents=True, exist_ok=False)
    latent_path = output / "embedding_latent.npy"
    cells_path = output / candidate[0]["cells_file"]
    np.save(latent_path, transformed)
    with DeterministicGzipTextWriter(cells_path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(cells.columns)
        writer.writerows(cells.itertuples(index=False, name=None))
    embedding = dict(candidate[0])
    embedding.update({
        "method": "proper_rigid_strict7_anchored_continual_export",
        "latent_file": latent_path.name,
        "cells_file": cells_path.name,
        "latent_sha256": sha256_path(latent_path),
        "cells_sha256": sha256_path(cells_path),
        "reference_coordinates_and_predictions_frozen": True,
    })
    write_json_exclusive(output / "embedding_manifest.json", embedding)
    anchor_digest = hashlib.sha256()
    for value in cells.iloc[anchor_positions]["cell_id"].astype(str):
        anchor_digest.update(value.encode("utf-8"))
        anchor_digest.update(b"\0")
    manifest = {
        "schema_version": "masld-cl-external-healthy-anchor-export-v23",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "parent_decision": {
            "path": str((policy_path.parent / policy["parent_decision"]["path"]).resolve()),
            "sha256": policy["parent_decision"]["sha256"],
            "decision": parent["decision"],
        },
        "candidate_embedding": {"path": str(candidate_path), "sha256": sha256_path(candidate_path)},
        "strict7_embedding": {"path": str(common_path), "sha256": sha256_path(common_path)},
        "anchor": {
            "n_cells": int(len(anchor_positions)),
            "cell_order_sha256": anchor_digest.hexdigest(),
            "source_mean": list(map(float, source_mean)),
            "target_mean": list(map(float, target_mean)),
            "rotation": rotation.tolist(),
            "rotation_determinant": float(np.linalg.det(rotation)),
            "orthogonality_maximum_error": ortho_error,
            "singular_values": list(map(float, singular)),
            "maximum_query_pair_distance_error": distance_error,
        },
        "strict7_coordinates_and_predictions_exact": True,
        "query_predictions_changed": False,
        "query_transform": "proper_rigid_rotation_plus_translation_no_scaling",
        "control_only": True,
        "case_stage_program_hero_gene_umap_cas13_read": False,
        "sensitivity_only": True,
        "embedding": embedding,
    }
    write_json_exclusive(output / "anchor_export_manifest.json", manifest)
    return manifest

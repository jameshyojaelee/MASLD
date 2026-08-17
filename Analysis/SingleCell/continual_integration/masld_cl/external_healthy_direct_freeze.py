"""Exact strict7 restoration with byte-identical V22 query representation."""

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


class ExternalHealthyDirectFreezeError(ContractError):
    """Raised when the no-fit V24 export changes identity."""


POLICY_SCHEMA = "masld-cl-external-healthy-direct-freeze-policy-v24"


def load_external_healthy_direct_freeze_policy(
    config: dict[str, Any], value: str | Path,
) -> tuple[Path, dict[str, Any], dict[str, Any]]:
    path = Path(value).resolve()
    expected = (
        Path(config["_config_path"]).resolve().parent
        / "reference" / "external_healthy_direct_freeze_policy_v24.json"
    )
    if path != expected:
        raise ExternalHealthyDirectFreezeError("V24 requires its source-controlled policy")
    with path.open() as handle:
        policy = json.load(handle)
    export = policy.get("export", {})
    if (
        policy.get("schema_version") != POLICY_SCHEMA
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("status") != "single_no_fit_stop_go_export"
        or policy.get("case_stage_program_hero_gene_umap_cas13_used") is not False
        or export.get("fitted_transform") is not False
        or export.get("hyperparameters") != []
        or policy.get("evaluation", {}).get("reuse_v22_thresholds_exactly") is not True
        or policy.get("evaluation", {}).get("no_additional_candidate_grid") is not True
    ):
        raise ExternalHealthyDirectFreezeError("V24 policy identity differs")
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
        raise ExternalHealthyDirectFreezeError("V24 parent decision differs")
    return path, policy, parent


def _source_path(policy_path: Path, spec: dict[str, Any]) -> Path:
    path = (policy_path.parent / spec["path"]).resolve()
    if sha256_path(path) != spec["sha256"]:
        raise ExternalHealthyDirectFreezeError("V24 source embedding changed")
    return path


def build_external_healthy_direct_freeze_export(
    config: dict[str, Any], policy_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    policy_path, policy, parent = load_external_healthy_direct_freeze_policy(
        config, policy_value
    )
    candidate_spec = policy["candidate_embedding"]
    common_spec = policy["strict7_embedding"]
    candidate_path = _source_path(policy_path, candidate_spec)
    common_path = _source_path(policy_path, common_spec)
    candidate_manifest_path = candidate_path.parent / "continual_reference_manifest.json"
    common_manifest_path = common_path.parent / "reference_manifest.json"
    if (
        sha256_path(candidate_manifest_path) != candidate_spec["training_manifest_sha256"]
        or sha256_path(common_manifest_path) != common_spec["training_manifest_sha256"]
    ):
        raise ExternalHealthyDirectFreezeError("V24 source training manifest changed")
    candidate = load_embedding(candidate_path)
    common = load_embedding(common_path)
    left, right = matched_rows(common[2], candidate[2])
    if len(left) != len(common[2]) or len(right) != len(candidate[2]):
        raise ExternalHealthyDirectFreezeError("V24 source cell rosters differ")
    cells = candidate[2].copy().reset_index(drop=True)
    common_latent = np.empty_like(np.asarray(candidate[1]))
    common_latent[right] = np.asarray(common[1])[left]
    common_cells = common[2].iloc[left].reset_index(drop=True)
    candidate_cells_sorted = cells.iloc[right].reset_index(drop=True)
    if not np.array_equal(
        common_cells["audit_cell_type"].astype(str).to_numpy(),
        candidate_cells_sorted["audit_cell_type"].astype(str).to_numpy(),
    ):
        raise ExternalHealthyDirectFreezeError("V24 frozen audit labels differ")
    common_predictions = np.empty(len(cells), dtype=object)
    common_predictions[right] = common_cells["predicted_cell_type"].astype(str).to_numpy()
    reference = cells["strict_reference"].to_numpy(dtype=bool)
    query = ~reference
    candidate_latent = np.asarray(candidate[1])
    latent = candidate_latent.copy()
    query_coordinates = latent[query].copy()
    query_coordinates_sha256 = hashlib.sha256(query_coordinates.tobytes()).hexdigest()
    query_predictions = cells.loc[query, "predicted_cell_type"].astype(str).to_numpy()
    latent[reference] = common_latent[reference]
    prediction_column = cells.columns.get_loc("predicted_cell_type")
    cells.iloc[np.flatnonzero(reference), prediction_column] = common_predictions[reference]
    if (
        not np.array_equal(latent[query], query_coordinates)
        or not np.array_equal(
            cells.loc[query, "predicted_cell_type"].astype(str).to_numpy(),
            query_predictions,
        )
        or not np.array_equal(latent[reference], common_latent[reference])
    ):
        raise ExternalHealthyDirectFreezeError("V24 exact freeze invariant failed")

    output = Path(output_value)
    output.mkdir(parents=True, exist_ok=False)
    latent_path = output / "embedding_latent.npy"
    cells_path = output / candidate[0]["cells_file"]
    np.save(latent_path, latent)
    with DeterministicGzipTextWriter(cells_path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(cells.columns)
        writer.writerows(cells.itertuples(index=False, name=None))
    embedding = dict(candidate[0])
    embedding.update({
        "method": "direct_strict7_frozen_continual_export",
        "latent_file": latent_path.name,
        "cells_file": cells_path.name,
        "latent_sha256": sha256_path(latent_path),
        "cells_sha256": sha256_path(cells_path),
        "reference_coordinates_and_predictions_frozen": True,
    })
    write_json_exclusive(output / "embedding_manifest.json", embedding)
    manifest = {
        "schema_version": "masld-cl-external-healthy-direct-freeze-export-v24",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "parent_decision": {
            "path": str((policy_path.parent / policy["parent_decision"]["path"]).resolve()),
            "sha256": policy["parent_decision"]["sha256"],
            "decision": parent["decision"],
        },
        "candidate_embedding": {"path": str(candidate_path), "sha256": sha256_path(candidate_path)},
        "strict7_embedding": {"path": str(common_path), "sha256": sha256_path(common_path)},
        "strict7_coordinates_and_predictions_exact": True,
        "query_coordinates_byte_identical_to_v22": True,
        "query_predictions_changed": False,
        "query_coordinates_sha256": query_coordinates_sha256,
        "fitted_transform": False,
        "control_only": True,
        "case_stage_program_hero_gene_umap_cas13_read": False,
        "sensitivity_only": True,
        "embedding": embedding,
    }
    write_json_exclusive(output / "anchor_export_manifest.json", manifest)
    return manifest

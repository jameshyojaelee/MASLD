"""Formal control-only pilot for the control-calibrated continual adapter."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .benchmark import _bundle_centroids, _paired_improvement, _validate_harmony_info
from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .control_evaluation import _neighbor_jaccard_loss, _reference_f1_change_ci
from .embedding import load_embedding, matched_rows
from .execution import source_identity, verify_source_identity_payload
from .firewall import assert_control_only_metric_names, validate_program_firewall


def evaluate_control_adapter_pilot(
    config: dict[str, Any], reference_embedding_value: str | Path,
    harmony_embedding_value: str | Path, architecture_embedding_value: str | Path,
    adapter_embedding_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    pipeline_root = Path(config["_config_path"]).resolve().parent
    adapter_path = Path(adapter_embedding_value).resolve()
    adapter_manifest_path = adapter_path.parent / "control_adapter_manifest.json"
    with adapter_manifest_path.open() as handle:
        adapter_manifest = json.load(handle)
    if (
        adapter_manifest.get("schema_version") != "masld-cl-control-adapter-v7"
        or adapter_manifest.get("config_sha256") != config["_config_sha256"]
        or adapter_manifest.get("control_only") is not True
        or adapter_manifest.get("outcomes_unlocked") is not False
        or adapter_manifest.get("adapter_source_identity") != source_identity(pipeline_root)
    ):
        raise ContractError("control adapter manifest or source identity is invalid")

    reference = load_embedding(reference_embedding_value)
    harmony = load_embedding(harmony_embedding_value)
    architecture = load_embedding(architecture_embedding_value)
    adapter = load_embedding(adapter_path)
    model_kind = adapter_manifest.get("model_kind")
    if (
        model_kind not in {"all_lineage", *config["lineages"]}
        or architecture[0].get("model_kind") != model_kind
        or reference[0].get("model_kind") != model_kind
        or adapter[0].get("model_kind") != model_kind
    ):
        raise ContractError("control adapter model kinds differ")
    _validate_harmony_info(harmony[0], config)
    if adapter_manifest.get("embedding") != adapter[0]:
        raise ContractError("adapter manifest does not own its embedding")
    left_base, right_adapter = matched_rows(architecture[2], adapter[2])
    if len(left_base) != len(architecture[2]) or not np.array_equal(
        architecture[2].iloc[left_base]["cell_id"].to_numpy(),
        adapter[2].iloc[right_adapter]["cell_id"].to_numpy(),
    ):
        raise ContractError("adapter changed the architecture cell roster")
    architecture_reference = architecture[2].iloc[left_base]["strict_reference"].to_numpy(dtype=bool)
    if not np.array_equal(
        architecture[2].iloc[left_base].loc[~architecture_reference, "predicted_cell_type"].astype(str).to_numpy(),
        adapter[2].iloc[right_adapter].loc[~architecture_reference, "predicted_cell_type"].astype(str).to_numpy(),
    ):
        raise ContractError("adapter changed architecture query predictions")

    reference_cells = reference[2].reset_index(drop=True)
    adapter_reference = adapter[2].loc[adapter[2]["strict_reference"]].reset_index(drop=True)
    left, right = matched_rows(reference_cells, adapter_reference)
    if not np.array_equal(
        np.asarray(reference[1])[left], np.asarray(adapter[1])[
            np.flatnonzero(adapter[2]["strict_reference"].to_numpy(dtype=bool))[right]
        ],
    ):
        raise ContractError("adapter did not preserve frozen reference coordinates")
    if not np.array_equal(
        reference_cells.iloc[left]["predicted_cell_type"].astype(str).to_numpy(),
        adapter_reference.iloc[right]["predicted_cell_type"].astype(str).to_numpy(),
    ):
        raise ContractError("adapter did not preserve frozen reference predictions")
    f1 = _reference_f1_change_ci(
        reference_cells.iloc[left].reset_index(drop=True),
        adapter_reference.iloc[right].reset_index(drop=True),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"],
    )
    adapter_reference_indices = np.flatnonzero(
        adapter[2]["strict_reference"].to_numpy(dtype=bool)
    )[right]
    jaccard = _neighbor_jaccard_loss(
        np.asarray(reference[1])[left], np.asarray(adapter[1])[adapter_reference_indices],
        reference_cells.iloc[left].reset_index(drop=True),
        k=config["evaluation"]["reference_neighborhood_k"],
        maximum_cells=config["evaluation"]["reference_neighborhood_max_cells"],
        seed=config["screen"]["seed"],
    )
    lineage = None if model_kind == "all_lineage" else model_kind
    harmony_result = _paired_improvement(
        _bundle_centroids(adapter, lineage), _bundle_centroids(harmony, lineage),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 4701,
    )
    architecture_result = _paired_improvement(
        _bundle_centroids(adapter, lineage), _bundle_centroids(architecture, lineage),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 4801,
    )
    protocol_harmony = _paired_improvement(
        _bundle_centroids(adapter, lineage, preparation="unsorted"),
        _bundle_centroids(harmony, lineage, preparation="unsorted"),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 4901,
    )
    protocol_architecture = _paired_improvement(
        _bundle_centroids(adapter, lineage, preparation="unsorted"),
        _bundle_centroids(architecture, lineage, preparation="unsorted"),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 5001,
    )
    metrics = {
        "reference_macro_f1_change_ci_low": f1["ci_low"],
        "reference_neighborhood_jaccard_loss": jaccard,
        "alignment_improvement_vs_harmony": harmony_result["improvement"],
        "alignment_improvement_vs_architecture_surgery": architecture_result["improvement"],
        "unsorted_alignment_improvement_vs_harmony": protocol_harmony["improvement"],
        "unsorted_alignment_improvement_vs_architecture_surgery": protocol_architecture["improvement"],
    }
    assert_control_only_metric_names(list(metrics))
    threshold = config["gates"]["minimum_control_alignment_improvement"]
    maximum_centered_error = max(
        value["maximum_centered_coordinate_error"]
        for value in adapter_manifest["offsets"]
    )
    gates = {
        "reference_macro_f1": metrics["reference_macro_f1_change_ci_low"]
        > -config["gates"]["reference_macro_f1_margin"],
        "reference_neighborhood": jaccard
        <= config["gates"]["reference_neighborhood_jaccard_loss"],
        "harmony_improvement": harmony_result["improvement"] >= threshold
        and harmony_result["ci_low"] > 0,
        "architecture_improvement": architecture_result["improvement"] >= threshold
        and architecture_result["ci_low"] > 0,
        "uniform_offset_invariant": maximum_centered_error <= 1e-5,
        "unsorted_harmony_no_material_worsening": protocol_harmony["improvement"] >= -0.10,
        "unsorted_architecture_no_material_worsening": protocol_architecture["improvement"] >= -0.10,
    }
    result = {
        "schema_version": "masld-cl-control-adapter-pilot-v7",
        "config_sha256": config["_config_sha256"],
        "control_only": True,
        "model_kind": model_kind,
        "global_reference_weight": adapter_manifest["global_reference_weight"],
        "outcomes_unlocked": False,
        "metrics": metrics,
        "bootstrap": {
            "harmony": harmony_result,
            "architecture_surgery": architecture_result,
            "unsorted_harmony": protocol_harmony,
            "unsorted_architecture_surgery": protocol_architecture,
        },
        "maximum_centered_coordinate_error": maximum_centered_error,
        "gates": gates,
        "passed": all(gates.values()),
        "sources": {
            "reference_embedding": str(Path(reference_embedding_value).resolve()),
            "harmony_embedding": str(Path(harmony_embedding_value).resolve()),
            "architecture_embedding": str(Path(architecture_embedding_value).resolve()),
            "adapter_embedding": str(adapter_path),
            "adapter_manifest": str(adapter_manifest_path),
            "adapter_manifest_sha256": sha256_path(adapter_manifest_path),
        },
    }
    write_json_exclusive(output_value, result)
    return result


def evaluate_orthogonal_bridge_pilot(
    config: dict[str, Any], reference_embedding_value: str | Path,
    harmony_embedding_value: str | Path, architecture_embedding_value: str | Path,
    candidate_embedding_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    """Control-only V8 evaluation, kept separate from the rejected V7 adapter."""
    validate_program_firewall(config)
    pipeline_root = Path(config["_config_path"]).resolve().parent
    candidate_path = Path(candidate_embedding_value).resolve()
    bridge_manifest_path = candidate_path.parent / "orthogonal_bridge_manifest.json"
    with bridge_manifest_path.open() as handle:
        bridge_manifest = json.load(handle)
    bridge_schema = bridge_manifest.get("schema_version")
    if (
        bridge_schema not in {
            "masld-cl-orthogonal-bridge-v8", "masld-cl-orthogonal-bridge-retarget-v9"
        }
        or bridge_manifest.get("config_sha256") != config["_config_sha256"]
        or bridge_manifest.get("control_only") is not True
        or bridge_manifest.get("outcomes_unlocked") is not False
    ):
        raise ContractError("orthogonal bridge manifest or source identity is invalid")
    verify_source_identity_payload(bridge_manifest.get("source_identity"))
    reference = load_embedding(reference_embedding_value)
    harmony = load_embedding(harmony_embedding_value)
    architecture = load_embedding(architecture_embedding_value)
    candidate = load_embedding(candidate_path)
    model_kind = bridge_manifest.get("model_kind")
    if (
        model_kind not in {"all_lineage", *config["lineages"]}
        or architecture[0].get("model_kind") != model_kind
        or reference[0].get("model_kind") != model_kind
        or candidate[0].get("model_kind") != model_kind
    ):
        raise ContractError("orthogonal bridge model kinds differ")
    _validate_harmony_info(harmony[0], config)
    if bridge_manifest.get("embedding") != candidate[0]:
        raise ContractError("orthogonal bridge manifest does not own its embedding")
    left_base, right_candidate = matched_rows(architecture[2], candidate[2])
    if len(left_base) != len(architecture[2]) or not np.array_equal(
        architecture[2].iloc[left_base]["cell_id"].to_numpy(),
        candidate[2].iloc[right_candidate]["cell_id"].to_numpy(),
    ):
        raise ContractError("orthogonal bridge changed the architecture cell roster")
    architecture_reference = architecture[2].iloc[left_base]["strict_reference"].to_numpy(dtype=bool)
    if not np.array_equal(
        architecture[2].iloc[left_base].loc[~architecture_reference, "predicted_cell_type"].astype(str).to_numpy(),
        candidate[2].iloc[right_candidate].loc[~architecture_reference, "predicted_cell_type"].astype(str).to_numpy(),
    ):
        raise ContractError("orthogonal bridge changed architecture query predictions")
    reference_cells = reference[2].reset_index(drop=True)
    candidate_reference = candidate[2].loc[candidate[2]["strict_reference"]].reset_index(drop=True)
    left, right = matched_rows(reference_cells, candidate_reference)
    candidate_reference_indices = np.flatnonzero(
        candidate[2]["strict_reference"].to_numpy(dtype=bool)
    )[right]
    if not np.array_equal(
        np.asarray(reference[1])[left], np.asarray(candidate[1])[candidate_reference_indices]
    ):
        raise ContractError("orthogonal bridge did not preserve frozen reference coordinates")
    if not np.array_equal(
        reference_cells.iloc[left]["predicted_cell_type"].astype(str).to_numpy(),
        candidate_reference.iloc[right]["predicted_cell_type"].astype(str).to_numpy(),
    ):
        raise ContractError("orthogonal bridge did not preserve frozen reference predictions")
    f1 = _reference_f1_change_ci(
        reference_cells.iloc[left].reset_index(drop=True),
        candidate_reference.iloc[right].reset_index(drop=True),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"],
    )
    jaccard = _neighbor_jaccard_loss(
        np.asarray(reference[1])[left], np.asarray(candidate[1])[candidate_reference_indices],
        reference_cells.iloc[left].reset_index(drop=True),
        k=config["evaluation"]["reference_neighborhood_k"],
        maximum_cells=config["evaluation"]["reference_neighborhood_max_cells"],
        seed=config["screen"]["seed"],
    )
    lineage = None if model_kind == "all_lineage" else model_kind
    harmony_result = _paired_improvement(
        _bundle_centroids(candidate, lineage), _bundle_centroids(harmony, lineage),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 5701,
    )
    architecture_result = _paired_improvement(
        _bundle_centroids(candidate, lineage), _bundle_centroids(architecture, lineage),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 5801,
    )
    protocol_harmony = _paired_improvement(
        _bundle_centroids(candidate, lineage, preparation="unsorted"),
        _bundle_centroids(harmony, lineage, preparation="unsorted"),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 5901,
    )
    protocol_architecture = _paired_improvement(
        _bundle_centroids(candidate, lineage, preparation="unsorted"),
        _bundle_centroids(architecture, lineage, preparation="unsorted"),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 6001,
    )
    metrics = {
        "reference_macro_f1_change_ci_low": f1["ci_low"],
        "reference_neighborhood_jaccard_loss": jaccard,
        "alignment_improvement_vs_harmony": harmony_result["improvement"],
        "alignment_improvement_vs_architecture_surgery": architecture_result["improvement"],
        "unsorted_alignment_improvement_vs_harmony": protocol_harmony["improvement"],
        "unsorted_alignment_improvement_vs_architecture_surgery": protocol_architecture["improvement"],
    }
    assert_control_only_metric_names(list(metrics))
    threshold = config["gates"]["minimum_control_alignment_improvement"]
    maximum_centered_error = max(
        value["maximum_centered_coordinate_error"] for value in bridge_manifest["offsets"]
    )
    gates = {
        "reference_macro_f1": metrics["reference_macro_f1_change_ci_low"]
        > -config["gates"]["reference_macro_f1_margin"],
        "reference_neighborhood": jaccard <= config["gates"]["reference_neighborhood_jaccard_loss"],
        "harmony_improvement": harmony_result["improvement"] >= threshold and harmony_result["ci_low"] > 0,
        "architecture_improvement": architecture_result["improvement"] >= threshold and architecture_result["ci_low"] > 0,
        "uniform_offset_invariant": maximum_centered_error <= 1e-5,
        "unsorted_harmony_no_material_worsening": protocol_harmony["improvement"] >= -0.10,
        "unsorted_architecture_no_material_worsening": protocol_architecture["improvement"] >= -0.10,
    }
    result = {
        "schema_version": (
            "masld-cl-orthogonal-bridge-retarget-pilot-v9"
            if bridge_schema == "masld-cl-orthogonal-bridge-retarget-v9"
            else "masld-cl-orthogonal-bridge-pilot-v8"
        ),
        "config_sha256": config["_config_sha256"], "control_only": True,
        "model_kind": model_kind,
        "control_offset_weight": bridge_manifest["control_offset_weight"],
        "global_reference_weight": bridge_manifest.get("global_reference_weight"),
        "outcomes_unlocked": False, "metrics": metrics,
        "bootstrap": {
            "harmony": harmony_result, "architecture_surgery": architecture_result,
            "unsorted_harmony": protocol_harmony,
            "unsorted_architecture_surgery": protocol_architecture,
        },
        "maximum_centered_coordinate_error": maximum_centered_error,
        "gates": gates, "passed": all(gates.values()),
        "sources": {
            "reference_embedding": str(Path(reference_embedding_value).resolve()),
            "harmony_embedding": str(Path(harmony_embedding_value).resolve()),
            "architecture_embedding": str(Path(architecture_embedding_value).resolve()),
            "candidate_embedding": str(candidate_path),
            "bridge_manifest": str(bridge_manifest_path),
            "bridge_manifest_sha256": sha256_path(bridge_manifest_path),
        },
    }
    write_json_exclusive(output_value, result)
    return result

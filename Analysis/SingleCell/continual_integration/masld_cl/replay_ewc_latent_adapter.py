"""Control-only adapter that retains the corrected replay-plus-EWC latent."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .backbone_diagnostic import _bundle_centroids, _paired_improvement
from .config import canonical_json_bytes, write_json_exclusive
from .contracts import ContractError, DeterministicGzipTextWriter, sha256_path
from .control_adapter import _donor_balanced_centroid
from .control_evaluation import (
    _neighbor_jaccard_loss,
    _reference_f1_change_ci,
    _shift_and_standard_error,
)
from .embedding import load_embedding, matched_rows
from .firewall import validate_program_firewall
from .orthogonal_bridge import (
    _balanced_fit_positions,
    apply_scaled_orthogonal_bridge,
    fit_scaled_orthogonal_bridge,
)


POLICY_SCHEMA = "masld-cl-replay-ewc-latent-adapter-policy-v27"


def _load_policy(config: dict[str, Any], value: str | Path) -> tuple[Path, dict[str, Any]]:
    path = Path(value).resolve()
    expected = (
        Path(config["_config_path"]).resolve().parent
        / "reference" / "replay_ewc_latent_adapter_policy_v27.json"
    )
    if path != expected:
        raise ContractError("V27 requires its source-controlled policy")
    with path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != POLICY_SCHEMA
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("diagnostic_basis", {}).get("query_labels_may_select_v27") is not False
        or policy.get("firewall", {}).get(
            "no_disease_or_held_study_metric_may_be_read_before_the_v27_control_lock"
        ) is not True
    ):
        raise ContractError("V27 policy identity or firewall differs")
    for section, pairs in {
        "parent": (("embedding", "embedding_sha256"), ("update_manifest", "update_manifest_sha256")),
        "reference": (("embedding", "embedding_sha256"), ("reference_manifest", "reference_manifest_sha256")),
        "diagnostic_basis": (("v26_sensitivity", "v26_sensitivity_sha256"),),
    }.items():
        for key, hash_key in pairs:
            source = (path.parent / policy[section][key]).resolve()
            if sha256_path(source) != policy[section][hash_key]:
                raise ContractError(f"V27 policy source changed: {section}.{key}")
    return path, policy


def _write_embedding(
    output: Path, parent_info: dict[str, Any], latent: np.ndarray, cells,
) -> dict[str, Any]:
    latent_path = output / "embedding_latent.npy"
    cells_path = output / parent_info["cells_file"]
    np.save(latent_path, latent)
    with DeterministicGzipTextWriter(cells_path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(cells.columns)
        writer.writerows(cells.itertuples(index=False, name=None))
    embedding = dict(parent_info)
    embedding.update({
        "method": "reference_frozen_replay_ewc_latent_adapter",
        "latent_file": latent_path.name,
        "cells_file": cells_path.name,
        "latent_sha256": sha256_path(latent_path),
        "cells_sha256": sha256_path(cells_path),
    })
    write_json_exclusive(output / "embedding_manifest.json", embedding)
    return embedding


def build_replay_ewc_latent_adapter_grid(
    config: dict[str, Any], policy_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    policy_path, policy = _load_policy(config, policy_value)
    parent_path = (policy_path.parent / policy["parent"]["embedding"]).resolve()
    parent_run_path = (policy_path.parent / policy["parent"]["update_manifest"]).resolve()
    reference_path = (policy_path.parent / policy["reference"]["embedding"]).resolve()
    reference_run_path = (
        policy_path.parent / policy["reference"]["reference_manifest"]
    ).resolve()
    parent_info, parent_latent, cells = load_embedding(parent_path)
    reference_info, reference_latent, reference_cells = load_embedding(reference_path)
    with parent_run_path.open() as handle:
        parent_run = json.load(handle)
    with reference_run_path.open() as handle:
        reference_run = json.load(handle)
    parent_spec = policy["parent"]
    if (
        parent_run.get("method") != parent_spec["required_method"]
        or float(parent_run.get("ewc_lambda", -1)) != parent_spec["required_ewc_lambda"]
        or float(parent_run.get("replay_fraction", -1)) != parent_spec["required_replay_fraction"]
        or parent_run.get("model_kind") != "all_lineage"
        or parent_run.get("embedding") != parent_info
        or reference_run.get("model_kind") != "all_lineage"
        or reference_run.get("embedding") != reference_info
    ):
        raise ContractError("V27 replay-plus-EWC parent or reference differs")

    reference_positions = np.flatnonzero(cells["strict_reference"].to_numpy(dtype=bool))
    ref_left, parent_ref_right = matched_rows(
        reference_cells, cells.iloc[reference_positions].reset_index(drop=True)
    )
    parent_reference_positions = reference_positions[parent_ref_right]
    if len(parent_reference_positions) != 216957:
        raise ContractError("V27 strict-reference roster differs")
    fit_positions = _balanced_fit_positions(
        cells, parent_reference_positions, "all_lineage", config
    )
    ref_by_cell = {
        cell: index for index, cell in enumerate(reference_cells["cell_id"].astype(str))
    }
    fit_reference_positions = np.asarray(
        [ref_by_cell[cell] for cell in cells.iloc[fit_positions]["cell_id"].astype(str)],
        dtype=np.int64,
    )
    source_mean, target_mean, rotation, scale = fit_scaled_orthogonal_bridge(
        np.asarray(parent_latent)[fit_positions],
        np.asarray(reference_latent)[fit_reference_positions],
    )
    bridged = apply_scaled_orthogonal_bridge(
        np.asarray(parent_latent), source_mean, target_mean, rotation, scale
    )
    bridged[parent_reference_positions] = np.asarray(reference_latent)[ref_left]
    cells = cells.copy()
    prediction_column = cells.columns.get_loc("predicted_cell_type")
    cells.iloc[parent_reference_positions, prediction_column] = (
        reference_cells.iloc[ref_left]["predicted_cell_type"].astype(str).to_numpy()
    )

    is_reference = cells["strict_reference"].to_numpy(dtype=bool)
    is_primary = cells["primary_query"].to_numpy(dtype=bool)
    is_control = is_primary & cells["query_control"].to_numpy(dtype=bool) & ~is_reference
    datasets = cells["dataset"].astype(str).to_numpy()
    preparations = cells["preparation"].astype(str).to_numpy()
    global_centroid, n_global = _donor_balanced_centroid(bridged, cells, is_reference)
    compatible = is_reference & (
        preparations == policy["control_calibration"]["compatible_reference_preparation"]
    )
    compatible_centroid, n_compatible = _donor_balanced_centroid(
        bridged, cells, compatible
    )
    if n_compatible < 3:
        raise ContractError("V27 compatible reference group is underpowered")

    output_root = Path(output_value).resolve()
    output_root.mkdir(parents=True, exist_ok=False)
    outputs = []
    for global_weight in map(
        float, policy["control_calibration"]["global_reference_weight_grid"]
    ):
        target = (
            (1.0 - global_weight) * compatible_centroid
            + global_weight * global_centroid
        )
        for offset_weight in map(
            float, policy["control_calibration"]["control_offset_weight_grid"]
        ):
            adapted = np.asarray(bridged).copy()
            offsets = []
            for dataset in policy["control_calibration"]["primary_query_datasets"]:
                control_mask = is_control & (datasets == dataset)
                control_centroid, n_control = _donor_balanced_centroid(
                    bridged, cells, control_mask
                )
                if n_control < policy["control_calibration"]["minimum_control_donors_per_dataset"]:
                    raise ContractError(f"V27 control group is underpowered: {dataset}")
                offset = offset_weight * (target - control_centroid)
                apply_mask = is_primary & (datasets == dataset)
                before = np.asarray(bridged[apply_mask], dtype=np.float64)
                adapted[apply_mask] = (before + offset).astype(adapted.dtype, copy=False)
                after = np.asarray(adapted[apply_mask], dtype=np.float64)
                centered_error = float(np.max(np.abs(
                    (after - after.mean(axis=0)) - (before - before.mean(axis=0))
                ), initial=0.0))
                offsets.append({
                    "dataset": dataset,
                    "n_control_donors": n_control,
                    "offset": list(map(float, offset)),
                    "maximum_centered_coordinate_error": centered_error,
                })
            if not np.array_equal(adapted[is_reference], bridged[is_reference]):
                raise ContractError("V27 calibration changed frozen reference coordinates")
            setting = (
                f"g{int(round(global_weight * 100)):03d}"
                f"_w{int(round(offset_weight * 100)):03d}"
            )
            output = output_root / setting
            output.mkdir()
            embedding = _write_embedding(output, parent_info, adapted, cells)
            manifest = {
                "schema_version": "masld-cl-replay-ewc-latent-adapter-v27",
                "config_sha256": config["_config_sha256"],
                "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
                "control_only": True,
                "outcomes_unlocked": False,
                "model_kind": "all_lineage",
                "global_reference_weight": global_weight,
                "control_offset_weight": offset_weight,
                "n_global_reference_donors": n_global,
                "n_compatible_reference_donors": n_compatible,
                "reference_coordinates_and_predictions_bitwise_frozen": True,
                "query_predictions_source": str(parent_path),
                "parent_embedding": {"path": str(parent_path), "sha256": sha256_path(parent_path)},
                "parent_run": {"path": str(parent_run_path), "sha256": sha256_path(parent_run_path)},
                "reference_embedding": {"path": str(reference_path), "sha256": sha256_path(reference_path)},
                "reference_run": {"path": str(reference_run_path), "sha256": sha256_path(reference_run_path)},
                "bridge": {
                    "fit_cells": int(len(fit_positions)),
                    "scale": scale,
                    "source_mean": source_mean.tolist(),
                    "target_mean": target_mean.tolist(),
                    "rotation": rotation.tolist(),
                    "orthonormality_max_abs_error": float(np.max(np.abs(
                        rotation.T @ rotation - np.eye(rotation.shape[1])
                    ))),
                },
                "offsets": offsets,
                "embedding": embedding,
            }
            manifest_path = output / "adapter_manifest.json"
            write_json_exclusive(manifest_path, manifest)
            outputs.append({
                "setting": setting,
                "global_reference_weight": global_weight,
                "control_offset_weight": offset_weight,
                "embedding": str((output / "embedding_manifest.json").resolve()),
                "manifest": str(manifest_path.resolve()),
            })
    grid = {
        "schema_version": "masld-cl-replay-ewc-latent-adapter-grid-v27",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "control_only": True,
        "outcomes_unlocked": False,
        "outputs": outputs,
    }
    write_json_exclusive(output_root / "grid_manifest.json", grid)
    return grid


def evaluate_and_select_replay_ewc_latent_adapter(
    config: dict[str, Any], policy_value: str | Path, grid_value: str | Path,
    reference_value: str | Path, harmony_value: str | Path,
    architecture_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    policy_path, policy = _load_policy(config, policy_value)
    grid_path = Path(grid_value).resolve()
    with grid_path.open() as handle:
        grid = json.load(handle)
    if (
        grid.get("schema_version") != "masld-cl-replay-ewc-latent-adapter-grid-v27"
        or grid.get("policy", {}).get("sha256") != sha256_path(policy_path)
        or grid.get("outcomes_unlocked") is not False
        or len(grid.get("outputs", [])) != 15
    ):
        raise ContractError("V27 grid identity differs")
    reference_path = Path(reference_value).resolve()
    harmony_path = Path(harmony_value).resolve()
    architecture_path = Path(architecture_value).resolve()
    reference = load_embedding(reference_path)
    harmony = load_embedding(harmony_path)
    architecture = load_embedding(architecture_path)
    if harmony[0].get("method") != "incumbent_scalesc_pca_harmony":
        raise ContractError("V27 Harmony comparator differs")

    results = []
    for index, source in enumerate(grid["outputs"]):
        embedding_path = Path(source["embedding"]).resolve()
        manifest_path = Path(source["manifest"]).resolve()
        candidate = load_embedding(embedding_path)
        with manifest_path.open() as handle:
            manifest = json.load(handle)
        if (
            manifest.get("schema_version") != "masld-cl-replay-ewc-latent-adapter-v27"
            or manifest.get("policy", {}).get("sha256") != sha256_path(policy_path)
            or manifest.get("outcomes_unlocked") is not False
            or manifest.get("embedding") != candidate[0]
        ):
            raise ContractError(f"V27 candidate manifest differs: {source['setting']}")
        positions = np.flatnonzero(candidate[2]["strict_reference"].to_numpy(dtype=bool))
        candidate_reference = candidate[2].iloc[positions].reset_index(drop=True)
        left, right = matched_rows(reference[2], candidate_reference)
        reference_cells = reference[2].iloc[left].reset_index(drop=True)
        candidate_cells = candidate_reference.iloc[right].reset_index(drop=True)
        reference_latent = np.asarray(reference[1])[left]
        candidate_latent = np.asarray(candidate[1])[positions[right]]
        exact_coordinates = np.array_equal(reference_latent, candidate_latent)
        exact_predictions = np.array_equal(
            reference_cells["predicted_cell_type"].astype(str).to_numpy(),
            candidate_cells["predicted_cell_type"].astype(str).to_numpy(),
        )
        if not exact_coordinates or not exact_predictions:
            raise ContractError(f"V27 reference export is not exact: {source['setting']}")
        f1 = _reference_f1_change_ci(
            reference_cells, candidate_cells, config["bootstrap"]["replicates"],
            config["bootstrap"]["seed"] + 2701 + index,
        )
        jaccard = _neighbor_jaccard_loss(
            reference_latent, candidate_latent, reference_cells,
            k=config["evaluation"]["reference_neighborhood_k"],
            maximum_cells=config["evaluation"]["reference_neighborhood_max_cells"],
            seed=config["screen"]["seed"] + 2701 + index,
        )
        harmony_result = _paired_improvement(
            _bundle_centroids(candidate), _bundle_centroids(harmony),
            config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 2801 + index,
        )
        architecture_result = _paired_improvement(
            _bundle_centroids(candidate), _bundle_centroids(architecture),
            config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 2901 + index,
        )
        shift, shift_se = _shift_and_standard_error(
            candidate[2], candidate[1], config["bootstrap"]["replicates"],
            config["bootstrap"]["seed"] + 3001 + index,
        )
        gates = {
            "reference_macro_f1": f1["ci_low"]
            > policy["selection"]["required_reference_macro_f1_change_ci_low"],
            "reference_neighborhood": jaccard
            <= policy["selection"]["maximum_reference_neighborhood_jaccard_loss"],
            "alignment_vs_harmony": harmony_result["improvement"]
            >= policy["selection"]["minimum_alignment_improvement_vs_harmony"]
            and harmony_result["ci_low"]
            > policy["selection"]["alignment_ci_low_must_exceed"],
            "alignment_vs_architecture_surgery": architecture_result["improvement"]
            >= config["gates"]["minimum_control_alignment_improvement"]
            and architecture_result["ci_low"] > 0,
        }
        results.append({
            **source,
            "reference_macro_f1_change": f1,
            "reference_neighborhood_jaccard_loss": jaccard,
            "reference_coordinates_bitwise_exact": exact_coordinates,
            "reference_predictions_bitwise_exact": exact_predictions,
            "shift_control": shift,
            "shift_control_standard_error": shift_se,
            "alignment_vs_harmony": harmony_result,
            "alignment_vs_architecture_surgery": architecture_result,
            "gates": gates,
            "eligible": all(gates.values()),
            "embedding_sha256": sha256_path(embedding_path),
            "manifest_sha256": sha256_path(manifest_path),
        })
    survivors = [row for row in results if row["eligible"]]
    selected = None
    threshold = None
    if survivors:
        best = min(survivors, key=lambda row: row["shift_control"])
        threshold = best["shift_control"] + best["shift_control_standard_error"]
        within_one_se = [row for row in survivors if row["shift_control"] <= threshold]
        within_one_se.sort(key=lambda row: (
            row["control_offset_weight"],
            abs(row["global_reference_weight"] - 0.75),
            row["global_reference_weight"],
            row["shift_control"],
        ))
        selected = within_one_se[0]
    decision = {
        "schema_version": "masld-cl-replay-ewc-latent-adapter-selection-v27",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "grid": {"path": str(grid_path), "sha256": sha256_path(grid_path)},
        "control_only": True,
        "outcome_variables_read": [],
        "query_audit_labels_used_to_construct_or_select": False,
        "selection_frozen": selected is not None,
        "outcomes_unlocked": selected is not None,
        "selected": selected,
        "one_standard_error_threshold": threshold,
        "eligible_settings": [row["setting"] for row in survivors],
        "results": results,
        "comparators": {
            "reference": {"path": str(reference_path), "sha256": sha256_path(reference_path)},
            "harmony": {"path": str(harmony_path), "sha256": sha256_path(harmony_path)},
            "architecture_surgery": {
                "path": str(architecture_path), "sha256": sha256_path(architecture_path)
            },
        },
    }
    decision["lock_sha256"] = hashlib.sha256(canonical_json_bytes(decision)).hexdigest()
    write_json_exclusive(output_value, decision)
    return decision

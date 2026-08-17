"""Class-routed geometry adapter guided by the replay-plus-EWC latent."""

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
from .latent_knn_label_audit import _fit_predict_knn, load_latent_knn_label_audit_policy
from .orthogonal_bridge import (
    _balanced_fit_positions,
    apply_scaled_orthogonal_bridge,
    fit_scaled_orthogonal_bridge,
)
from .training import _capped_indices


def _load_policy(config: dict[str, Any], value: str | Path) -> tuple[Path, dict[str, Any]]:
    path = Path(value).resolve()
    reference_root = Path(config["_config_path"]).resolve().parent / "reference"
    profiles = {
        reference_root / "class_routed_replay_ewc_adapter_policy_v28.json": (
            "masld-cl-class-routed-replay-ewc-adapter-policy-v28",
            "no_disease_or_held_study_metric_may_be_read_before_the_v28_control_lock",
        ),
        reference_root / "compartment_routed_replay_ewc_adapter_policy_v29.json": (
            "masld-cl-compartment-routed-replay-ewc-adapter-policy-v29",
            "no_disease_or_held_study_metric_may_be_read_before_the_v29_control_lock",
        ),
        reference_root / "immune_compartment_replay_ewc_adapter_policy_v30.json": (
            "masld-cl-immune-compartment-replay-ewc-adapter-policy-v30",
            "no_disease_or_held_study_metric_may_be_read_before_the_v30_control_lock",
        ),
    }
    profile = profiles.get(path)
    if profile is None:
        raise ContractError("class-routed adapter requires a source-controlled policy")
    with path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != profile[0]
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("routing_classifier", {}).get("query_labels_available_to_fit_or_route") is not False
        or policy.get("firewall", {}).get(
            profile[1]
        ) is not True
    ):
        raise ContractError("class-routed adapter policy identity or firewall differs")
    sources = (
        ("parent", "embedding", "embedding_sha256"),
        ("reference", "embedding", "embedding_sha256"),
        ("raw_pca", "embedding", "embedding_sha256"),
        ("diagnostic_basis", "v26", "v26_sha256"),
        ("diagnostic_basis", "v27_outcomes", "v27_outcomes_sha256"),
    )
    for section, key, hash_key in sources:
        source = (path.parent / policy[section][key]).resolve()
        if sha256_path(source) != policy[section][hash_key]:
            raise ContractError(f"V28 policy source changed: {section}.{key}")
    if "v28_outcomes" in policy["diagnostic_basis"]:
        source = (path.parent / policy["diagnostic_basis"]["v28_outcomes"]).resolve()
        if sha256_path(source) != policy["diagnostic_basis"]["v28_outcomes_sha256"]:
            raise ContractError("V29 policy source changed: diagnostic_basis.v28_outcomes")
    if "v29_outcomes" in policy["diagnostic_basis"]:
        source = (path.parent / policy["diagnostic_basis"]["v29_outcomes"]).resolve()
        if sha256_path(source) != policy["diagnostic_basis"]["v29_outcomes_sha256"]:
            raise ContractError("V30 policy source changed: diagnostic_basis.v29_outcomes")
    classifier_path = (path.parent / policy["routing_classifier"]["policy"]).resolve()
    if sha256_path(classifier_path) != policy["routing_classifier"]["policy_sha256"]:
        raise ContractError("V28 routing classifier policy changed")
    return path, policy


def _ordered_latent(bundle_cells, source_cells, source_latent) -> np.ndarray:
    left, right = matched_rows(source_cells, bundle_cells)
    order = np.argsort(right)
    left, right = left[order], right[order]
    if not np.array_equal(right, np.arange(len(bundle_cells))):
        raise ContractError("V28 cell rosters differ")
    return np.asarray(source_latent)[left]


def _write_embedding(
    output: Path, parent_info: dict[str, Any], latent, cells, method: str,
):
    latent_path = output / "embedding_latent.npy"
    cells_path = output / parent_info["cells_file"]
    np.save(latent_path, latent)
    with DeterministicGzipTextWriter(cells_path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(cells.columns)
        writer.writerows(cells.itertuples(index=False, name=None))
    info = dict(parent_info)
    info.update({
        "method": method,
        "latent_file": latent_path.name,
        "cells_file": cells_path.name,
        "latent_sha256": sha256_path(latent_path),
        "cells_sha256": sha256_path(cells_path),
    })
    write_json_exclusive(output / "embedding_manifest.json", info)
    return info


def build_class_routed_adapter_grid(
    config: dict[str, Any], policy_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    policy_path, policy = _load_policy(config, policy_value)
    version = policy["schema_version"].rsplit("v", 1)[-1]
    is_v29 = version == "29"
    is_v30 = version == "30"
    method = (
        "immune_compartment_replay_ewc_geometry_adapter"
        if is_v30 else
        "compartment_routed_replay_ewc_geometry_adapter"
        if is_v29 else "class_routed_replay_ewc_geometry_adapter"
    )
    candidate_schema = (
        "masld-cl-immune-compartment-replay-ewc-adapter-v30"
        if is_v30 else
        "masld-cl-compartment-routed-replay-ewc-adapter-v29"
        if is_v29 else "masld-cl-class-routed-replay-ewc-adapter-v28"
    )
    grid_schema = (
        "masld-cl-immune-compartment-replay-ewc-adapter-grid-v30"
        if is_v30 else
        "masld-cl-compartment-routed-replay-ewc-adapter-grid-v29"
        if is_v29 else "masld-cl-class-routed-replay-ewc-adapter-grid-v28"
    )
    parent_path = (policy_path.parent / policy["parent"]["embedding"]).resolve()
    reference_path = (policy_path.parent / policy["reference"]["embedding"]).resolve()
    raw_path = (policy_path.parent / policy["raw_pca"]["embedding"]).resolve()
    parent_info, parent_latent, cells = load_embedding(parent_path)
    _, reference_latent, reference_cells = load_embedding(reference_path)
    _, raw_latent, raw_cells = load_embedding(raw_path)
    if parent_info.get("model_kind") != "all_lineage":
        raise ContractError("V28 parent is not all-lineage")
    raw_values = _ordered_latent(cells, raw_cells, raw_latent)
    reference_positions = np.flatnonzero(cells["strict_reference"].to_numpy(dtype=bool))
    query_positions = np.flatnonzero(
        cells["analysis_eligible"].to_numpy(dtype=bool)
        & ~cells["strict_reference"].to_numpy(dtype=bool)
    )
    ref_left, parent_ref_right = matched_rows(
        reference_cells, cells.iloc[reference_positions].reset_index(drop=True)
    )
    parent_reference_positions = reference_positions[parent_ref_right]
    reference_values = np.asarray(reference_latent)[ref_left]
    if len(parent_reference_positions) != 216957 or len(query_positions) != 687559:
        raise ContractError("V28 reference or query roster differs")

    classifier_path = (policy_path.parent / policy["routing_classifier"]["policy"]).resolve()
    _, classifier_policy = load_latent_knn_label_audit_policy(config, classifier_path)
    reference_pool = _capped_indices(
        type("ADataView", (), {"obs": cells, "n_obs": len(cells)})(),
        reference_positions, "all_lineage", config,
        classifier_policy["classifier"]["reference_cap_seed"],
    )
    query_label_routes, route_identity = _fit_predict_knn(
        np.asarray(parent_latent)[reference_pool],
        cells.iloc[reference_pool]["audit_cell_type"].astype(str).to_numpy(),
        np.asarray(parent_latent)[query_positions], classifier_policy["classifier"],
    )
    reference_label_routes = (
        cells.iloc[parent_reference_positions]["audit_cell_type"].astype(str).to_numpy()
    )
    routing_groups = policy.get("routing_groups")
    if routing_groups:
        lookup = {
            label: group for group, labels in routing_groups.items() for label in labels
        }
        expected_labels = set(reference_label_routes)
        if set(lookup) != expected_labels or len(lookup) != sum(map(len, routing_groups.values())):
            raise ContractError("V29 routing-group registry differs from reference labels")
        reference_routes = np.asarray([lookup[label] for label in reference_label_routes], dtype=object)
        query_routes = np.asarray([lookup[str(label)] for label in query_label_routes], dtype=object)
    else:
        reference_routes = reference_label_routes
        query_routes = query_label_routes
    routes = np.empty(len(cells), dtype=object)
    routes[:] = None
    routes[parent_reference_positions] = reference_routes
    routes[query_positions] = query_routes
    routing_labels = np.empty(len(cells), dtype=object)
    routing_labels[:] = None
    routing_labels[parent_reference_positions] = reference_label_routes
    routing_labels[query_positions] = query_label_routes
    if any(value is None for value in routes):
        raise ContractError("V28 routing left cells unassigned")
    classes = sorted(set(map(str, routes[parent_reference_positions])))
    if set(map(str, query_routes)) - set(classes):
        raise ContractError("V28 query route is absent from the reference")

    bridged = np.empty_like(raw_values, dtype=np.float32)
    bridge_summaries = []
    ref_by_cell = {
        cell: index for index, cell in enumerate(reference_cells["cell_id"].astype(str))
    }
    for class_name in classes:
        class_reference = parent_reference_positions[routes[parent_reference_positions] == class_name]
        fit_positions = _balanced_fit_positions(
            cells, class_reference, "all_lineage", config
        )
        fit_reference_positions = np.asarray(
            [ref_by_cell[cell] for cell in cells.iloc[fit_positions]["cell_id"].astype(str)],
            dtype=np.int64,
        )
        source_mean, target_mean, rotation, scale = fit_scaled_orthogonal_bridge(
            raw_values[fit_positions], np.asarray(reference_latent)[fit_reference_positions]
        )
        apply_positions = np.flatnonzero(routes == class_name)
        bridged[apply_positions] = apply_scaled_orthogonal_bridge(
            raw_values[apply_positions], source_mean, target_mean, rotation, scale
        )
        bridge_summaries.append({
            "class": class_name, "fit_cells": int(len(fit_positions)),
            "reference_cells": int(len(class_reference)),
            "query_cells": int(np.sum(routes[query_positions] == class_name)),
            "scale": scale, "source_mean": source_mean.tolist(),
            "target_mean": target_mean.tolist(), "rotation": rotation.tolist(),
        })
    bridged[parent_reference_positions] = reference_values
    cells = cells.copy()
    prediction_column = cells.columns.get_loc("predicted_cell_type")
    cells.iloc[parent_reference_positions, prediction_column] = (
        reference_cells.iloc[ref_left]["predicted_cell_type"].astype(str).to_numpy()
    )
    cells["routing_class"] = routes
    cells["routing_label"] = routing_labels

    output_root = Path(output_value).resolve()
    output_root.mkdir(parents=True, exist_ok=False)
    outputs = []
    primary = cells["primary_query"].to_numpy(dtype=bool)
    controls = primary & cells["query_control"].to_numpy(dtype=bool)
    datasets = cells["dataset"].astype(str).to_numpy()
    minimum = int(policy["control_calibration"]["minimum_control_donors_per_dataset_class"])
    for weight in map(float, policy["control_calibration"]["control_offset_weight_grid"]):
        adapted = bridged.copy()
        offsets = []
        for class_name in classes:
            reference_mask = cells["strict_reference"].to_numpy(dtype=bool) & (routes == class_name)
            target, n_reference = _donor_balanced_centroid(bridged, cells, reference_mask)
            for dataset in policy["control_calibration"]["primary_query_datasets"]:
                control_mask = controls & (datasets == dataset) & (routes == class_name)
                donors = cells.loc[control_mask, "donor_id"].astype(str).nunique()
                apply_mask = primary & (datasets == dataset) & (routes == class_name)
                if donors < minimum:
                    offsets.append({
                        "class": class_name, "dataset": dataset, "translated": False,
                        "n_control_donors": int(donors), "n_reference_donors": n_reference,
                        "n_applied_cells": int(apply_mask.sum()),
                    })
                    continue
                control_centroid, _ = _donor_balanced_centroid(
                    bridged, cells, control_mask
                )
                offset = weight * (target - control_centroid)
                before = np.asarray(bridged[apply_mask], dtype=np.float64)
                adapted[apply_mask] = (before + offset).astype(np.float32)
                after = np.asarray(adapted[apply_mask], dtype=np.float64)
                centered_error = float(np.max(np.abs(
                    (after - after.mean(axis=0)) - (before - before.mean(axis=0))
                ), initial=0.0))
                offsets.append({
                    "class": class_name, "dataset": dataset, "translated": True,
                    "n_control_donors": int(donors), "n_reference_donors": n_reference,
                    "n_applied_cells": int(apply_mask.sum()), "offset": list(map(float, offset)),
                    "maximum_centered_coordinate_error": centered_error,
                })
        if not np.array_equal(adapted[parent_reference_positions], reference_values):
            raise ContractError("V28 changed frozen reference coordinates")
        setting = f"w{int(round(weight * 100)):03d}"
        output = output_root / setting
        output.mkdir()
        embedding = _write_embedding(output, parent_info, adapted, cells, method)
        manifest = {
            "schema_version": candidate_schema,
            "config_sha256": config["_config_sha256"],
            "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
            "control_only": True, "outcomes_unlocked": False,
            "control_offset_weight": weight, "classes": classes,
            "query_labels_available_to_fit_or_route": False,
            "routing_identity": route_identity,
            "routing_groups": routing_groups,
            "bridges": bridge_summaries, "offsets": offsets,
            "embedding": embedding,
        }
        manifest_path = output / "adapter_manifest.json"
        write_json_exclusive(manifest_path, manifest)
        outputs.append({
            "setting": setting, "control_offset_weight": weight,
            "embedding": str((output / "embedding_manifest.json").resolve()),
            "manifest": str(manifest_path.resolve()),
        })
    grid = {
        "schema_version": grid_schema,
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "control_only": True, "outcomes_unlocked": False, "outputs": outputs,
    }
    write_json_exclusive(output_root / "grid_manifest.json", grid)
    return grid


def evaluate_and_select_class_routed_adapter(
    config: dict[str, Any], policy_value: str | Path, grid_value: str | Path,
    reference_value: str | Path, harmony_value: str | Path,
    architecture_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    validate_program_firewall(config)
    policy_path, policy = _load_policy(config, policy_value)
    version = policy["schema_version"].rsplit("v", 1)[-1]
    is_v29 = version == "29"
    is_v30 = version == "30"
    grid_schema = (
        "masld-cl-immune-compartment-replay-ewc-adapter-grid-v30"
        if is_v30 else
        "masld-cl-compartment-routed-replay-ewc-adapter-grid-v29"
        if is_v29 else "masld-cl-class-routed-replay-ewc-adapter-grid-v28"
    )
    candidate_schema = (
        "masld-cl-immune-compartment-replay-ewc-adapter-v30"
        if is_v30 else
        "masld-cl-compartment-routed-replay-ewc-adapter-v29"
        if is_v29 else "masld-cl-class-routed-replay-ewc-adapter-v28"
    )
    selection_schema = (
        "masld-cl-immune-compartment-replay-ewc-adapter-selection-v30"
        if is_v30 else
        "masld-cl-compartment-routed-replay-ewc-adapter-selection-v29"
        if is_v29 else "masld-cl-class-routed-replay-ewc-adapter-selection-v28"
    )
    grid_path = Path(grid_value).resolve()
    with grid_path.open() as handle:
        grid = json.load(handle)
    if (
        grid.get("schema_version") != grid_schema
        or grid.get("policy", {}).get("sha256") != sha256_path(policy_path)
        or grid.get("outcomes_unlocked") is not False
        or len(grid.get("outputs", [])) != 3
    ):
        raise ContractError("V28 grid identity differs")
    reference_path = Path(reference_value).resolve()
    harmony_path = Path(harmony_value).resolve()
    architecture_path = Path(architecture_value).resolve()
    reference = load_embedding(reference_path)
    harmony = load_embedding(harmony_path)
    architecture = load_embedding(architecture_path)
    results = []
    for index, source in enumerate(grid["outputs"]):
        embedding_path = Path(source["embedding"]).resolve()
        manifest_path = Path(source["manifest"]).resolve()
        candidate = load_embedding(embedding_path)
        with manifest_path.open() as handle:
            manifest = json.load(handle)
        if (
            manifest.get("schema_version") != candidate_schema
            or manifest.get("policy", {}).get("sha256") != sha256_path(policy_path)
            or manifest.get("outcomes_unlocked") is not False
            or manifest.get("embedding") != candidate[0]
        ):
            raise ContractError(f"V28 candidate differs: {source['setting']}")
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
            raise ContractError(f"V28 reference export is not exact: {source['setting']}")
        f1 = _reference_f1_change_ci(
            reference_cells, candidate_cells, config["bootstrap"]["replicates"],
            config["bootstrap"]["seed"] + 3101 + index,
        )
        jaccard = _neighbor_jaccard_loss(
            reference_latent, candidate_latent, reference_cells,
            k=config["evaluation"]["reference_neighborhood_k"],
            maximum_cells=config["evaluation"]["reference_neighborhood_max_cells"],
            seed=config["screen"]["seed"] + 3101 + index,
        )
        harmony_result = _paired_improvement(
            _bundle_centroids(candidate), _bundle_centroids(harmony),
            config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 3201 + index,
        )
        architecture_result = _paired_improvement(
            _bundle_centroids(candidate), _bundle_centroids(architecture),
            config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + 3301 + index,
        )
        shift, shift_se = _shift_and_standard_error(
            candidate[2], candidate[1], config["bootstrap"]["replicates"],
            config["bootstrap"]["seed"] + 3401 + index,
        )
        gates = {
            "reference_macro_f1": f1["ci_low"]
            > policy["selection"]["required_reference_macro_f1_change_ci_low"],
            "reference_neighborhood": jaccard
            <= policy["selection"]["maximum_reference_neighborhood_jaccard_loss"],
            "alignment_vs_harmony": harmony_result["improvement"]
            >= policy["selection"]["minimum_alignment_improvement_vs_harmony"]
            and harmony_result["ci_low"] > policy["selection"]["alignment_ci_low_must_exceed"],
            "alignment_vs_architecture_surgery": architecture_result["improvement"]
            >= config["gates"]["minimum_control_alignment_improvement"]
            and architecture_result["ci_low"] > 0,
        }
        results.append({
            **source, "reference_macro_f1_change": f1,
            "reference_neighborhood_jaccard_loss": jaccard,
            "reference_coordinates_bitwise_exact": exact_coordinates,
            "reference_predictions_bitwise_exact": exact_predictions,
            "shift_control": shift, "shift_control_standard_error": shift_se,
            "alignment_vs_harmony": harmony_result,
            "alignment_vs_architecture_surgery": architecture_result,
            "gates": gates, "eligible": all(gates.values()),
            "embedding_sha256": sha256_path(embedding_path),
            "manifest_sha256": sha256_path(manifest_path),
        })
    survivors = [row for row in results if row["eligible"]]
    selected = None
    threshold = None
    if survivors:
        best = min(survivors, key=lambda row: row["shift_control"])
        threshold = best["shift_control"] + best["shift_control_standard_error"]
        selected = min(
            (row for row in survivors if row["shift_control"] <= threshold),
            key=lambda row: (row["control_offset_weight"], row["shift_control"]),
        )
    decision = {
        "schema_version": selection_schema,
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "grid": {"path": str(grid_path), "sha256": sha256_path(grid_path)},
        "control_only": True, "outcome_variables_read": [],
        "query_audit_labels_used_to_fit_route_construct_or_select": False,
        "selection_frozen": selected is not None, "outcomes_unlocked": selected is not None,
        "selected": selected, "one_standard_error_threshold": threshold,
        "eligible_settings": [row["setting"] for row in survivors], "results": results,
        "comparators": {
            "reference": {"path": str(reference_path), "sha256": sha256_path(reference_path)},
            "harmony": {"path": str(harmony_path), "sha256": sha256_path(harmony_path)},
            "architecture_surgery": {"path": str(architecture_path), "sha256": sha256_path(architecture_path)},
        },
    }
    decision["lock_sha256"] = hashlib.sha256(canonical_json_bytes(decision)).hexdigest()
    write_json_exclusive(output_value, decision)
    return decision

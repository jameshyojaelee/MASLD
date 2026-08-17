"""Prediction and acquisition-order invariants for the locked V9 bridge."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .benchmark import _label_scores
from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .embedding import load_embedding, matched_rows
from .firewall import load_retargeted_bridge_selection_lock, validate_program_firewall


def _selected_results(selection: dict[str, Any]) -> dict[str, dict[str, Any]]:
    weight = float(selection["selected_global_reference_weight"])
    result = {}
    for source in selection["control_results"]:
        with Path(source["path"]).open() as handle:
            row = json.load(handle)
        if float(row["global_reference_weight"]) == weight:
            result[row["model_kind"]] = row
    return result


def _align_complete(left_bundle, right_bundle):
    left_info, left_latent, left_cells = left_bundle
    right_info, right_latent, right_cells = right_bundle
    left_rows, right_rows = matched_rows(left_cells, right_cells)
    if len(left_rows) != len(left_cells) or len(right_rows) != len(right_cells):
        raise ContractError("invariant comparison requires identical cell rosters")
    order = np.argsort(left_rows)
    left_rows, right_rows = left_rows[order], right_rows[order]
    if not np.array_equal(left_rows, np.arange(len(left_cells))):
        raise ContractError("invariant comparison failed to recover left cell order")
    return (
        left_info, np.asarray(left_latent), left_cells,
        right_info, np.asarray(right_latent[right_rows]), right_cells.iloc[right_rows].reset_index(drop=True),
    )


def evaluate_prediction_invariant(
    config: dict[str, Any], selection_lock: str | Path, output: str | Path,
) -> dict[str, Any]:
    selection_path = Path(selection_lock).resolve()
    selection = load_retargeted_bridge_selection_lock(selection_path, config)
    validate_program_firewall(config)
    selected = _selected_results(selection)
    expected = {"all_lineage", *config["lineages"]}
    if set(selected) != expected:
        raise ContractError("prediction invariant lacks the exact selected model roster")
    models, sources = {}, []
    for kind, result in selected.items():
        candidate_path = Path(result["sources"]["candidate_embedding"]).resolve()
        architecture_path = Path(result["sources"]["architecture_embedding"]).resolve()
        candidate = load_embedding(candidate_path)
        architecture = load_embedding(architecture_path)
        _, _, candidate_cells, _, _, architecture_cells = _align_complete(candidate, architecture)
        query = candidate_cells["analysis_eligible"].to_numpy(dtype=bool) & ~candidate_cells["strict_reference"].to_numpy(dtype=bool)
        prediction_equal = np.array_equal(
            candidate_cells.loc[query, "predicted_cell_type"].astype(str).to_numpy(),
            architecture_cells.loc[query, "predicted_cell_type"].astype(str).to_numpy(),
        )
        audit_equal = np.array_equal(
            candidate_cells.loc[query, "audit_cell_type"].astype(str).to_numpy(),
            architecture_cells.loc[query, "audit_cell_type"].astype(str).to_numpy(),
        )
        if not prediction_equal or not audit_equal:
            raise ContractError("V9 query labels differ from architecture surgery")
        evaluated_lineages = config["lineages"] if kind == "all_lineage" else [kind]
        candidate_scores = _label_scores((candidate[0], candidate[1], candidate_cells), evaluated_lineages)
        architecture_scores = _label_scores((architecture[0], architecture[1], architecture_cells), evaluated_lineages)
        overall_change = candidate_scores[0] - architecture_scores[0]
        lineage_changes = {
            lineage: candidate_scores[1][lineage] - architecture_scores[1][lineage]
            for lineage in evaluated_lineages
        }
        models[kind] = {
            "n_query_cells": int(query.sum()), "predictions_bitwise_equal": True,
            "audit_labels_equal": True, "query_macro_f1_change_vs_architecture_surgery": overall_change,
            "major_lineage_f1_changes_vs_architecture_surgery": lineage_changes,
            "architecture_noninferiority_pass": bool(
                overall_change >= -config["gates"]["query_macro_f1_margin"]
                and min(lineage_changes.values()) >= -config["gates"]["major_lineage_f1_margin"]
            ),
        }
        sources.append({
            "model_kind": kind,
            "candidate_embedding": str(candidate_path), "candidate_embedding_sha256": sha256_path(candidate_path),
            "architecture_embedding": str(architecture_path), "architecture_embedding_sha256": sha256_path(architecture_path),
        })
    result = {
        "schema_version": "masld-cl-orthogonal-bridge-prediction-invariant-v9",
        "config_sha256": config["_config_sha256"],
        "selection_lock": str(selection_path), "selection_lock_file_sha256": sha256_path(selection_path),
        "selection_lock_sha256": selection["lock_sha256"], "models": models, "sources": sources,
        "architecture_noninferiority_pass": bool(all(row["architecture_noninferiority_pass"] for row in models.values())),
        "best_non_cl_gate_status": "pending_de_novo_fine_tune_replay_only_and_ewc_only_comparators",
        "labels_are_audit_targets_and_were_not_overwritten": True,
    }
    write_json_exclusive(output, result)
    return result


def apply_offsets_in_order(
    parent: np.ndarray, datasets: np.ndarray, offsets: dict[str, np.ndarray], order: list[str],
) -> np.ndarray:
    result = np.asarray(parent).copy()
    seen = set()
    for dataset in order:
        if dataset not in offsets:
            continue
        if dataset in seen:
            raise ContractError("acquisition order contains a duplicate mapped dataset")
        seen.add(dataset)
        mask = datasets == dataset
        result[mask] = (np.asarray(parent[mask], dtype=np.float64) + offsets[dataset]).astype(result.dtype)
    if seen != set(offsets):
        raise ContractError("acquisition order omitted a mapped dataset")
    return result


def evaluate_order_invariant(
    config: dict[str, Any], selection_lock: str | Path,
    outcome_decision: str | Path, output: str | Path,
) -> dict[str, Any]:
    selection_path = Path(selection_lock).resolve()
    selection = load_retargeted_bridge_selection_lock(selection_path, config)
    validate_program_firewall(config)
    decision_path = Path(outcome_decision).resolve()
    with decision_path.open() as handle:
        decision = json.load(handle)
    if (
        decision.get("schema_version") != "masld-cl-orthogonal-bridge-outcome-decision-v9"
        or decision.get("selection_lock_sha256") != selection["lock_sha256"]
        or decision.get("disease_preservation_pass") is not True
    ):
        raise ContractError("order evaluation requires the passing locked V9 outcome decision")
    selected = _selected_results(selection)
    expected = {"all_lineage", *config["lineages"]}
    if set(selected) != expected:
        raise ContractError("order invariant lacks the exact selected model roster")
    models, sources = {}, []
    for kind, control in selected.items():
        candidate_path = Path(control["sources"]["candidate_embedding"]).resolve()
        bridge_path = Path(control["sources"]["bridge_manifest"]).resolve()
        with bridge_path.open() as handle:
            bridge = json.load(handle)
        parent_path = Path(bridge["parent_embedding_realpath"]).resolve()
        if sha256_path(parent_path) != bridge["parent_embedding_sha256"]:
            raise ContractError("order invariant parent embedding changed")
        candidate = load_embedding(candidate_path)
        parent = load_embedding(parent_path)
        _, candidate_latent, cells, _, parent_latent, parent_cells = _align_complete(candidate, parent)
        if not np.array_equal(cells["cell_id"].astype(str), parent_cells["cell_id"].astype(str)):
            raise ContractError("order invariant cell rosters differ")
        datasets = cells["dataset"].astype(str).to_numpy()
        offsets = {row["dataset"]: np.asarray(row["offset"], dtype=np.float64) for row in bridge["offsets"]}
        order_rows = []
        for order_name, order in config["acquisition_orders"].items():
            sequential = apply_offsets_in_order(parent_latent, datasets, offsets, order)
            maximum_error = float(np.max(np.abs(
                np.asarray(sequential, dtype=np.float64) - np.asarray(candidate_latent, dtype=np.float64)
            ), initial=0.0))
            order_rows.append({
                "order": order_name, "maximum_coordinate_error_vs_joint": maximum_error,
                "reference_and_disease_pass_inherited_from_identical_joint_coordinates": bool(maximum_error <= 1e-6),
                "control_alignment_pass_inherited_from_identical_joint_coordinates": bool(maximum_error <= 1e-6),
                "donor_centroid_distance_spearman_vs_joint": 1.0 if maximum_error <= 1e-6 else None,
            })
        models[kind] = {"orders": order_rows, "all_six_coordinate_equivalent": all(
            row["maximum_coordinate_error_vs_joint"] <= 1e-6 for row in order_rows
        )}
        sources.append({
            "model_kind": kind, "candidate_embedding": str(candidate_path),
            "candidate_embedding_sha256": sha256_path(candidate_path),
            "parent_embedding": str(parent_path), "parent_embedding_sha256": sha256_path(parent_path),
            "bridge_manifest": str(bridge_path), "bridge_manifest_sha256": sha256_path(bridge_path),
        })
    equivalent = bool(all(row["all_six_coordinate_equivalent"] for row in models.values()))
    result = {
        "schema_version": "masld-cl-orthogonal-bridge-order-invariant-v9",
        "config_sha256": config["_config_sha256"],
        "selection_lock": str(selection_path), "selection_lock_file_sha256": sha256_path(selection_path),
        "selection_lock_sha256": selection["lock_sha256"],
        "outcome_decision": str(decision_path), "outcome_decision_sha256": sha256_path(decision_path),
        "models": models, "sources": sources, "all_six_orders_coordinate_equivalent": equivalent,
        "passing_order_count": 6 if equivalent else 0,
        "worst_order_degradation": 0.0 if equivalent else None,
        "minimum_order_distance_spearman": 1.0 if equivalent else None,
        "order_variability": 0.0 if equivalent else None,
        "order_to_de_novo_seed_variability_ratio": 0.0 if equivalent else None,
        "order_gate_pass": equivalent,
    }
    write_json_exclusive(output, result)
    return result

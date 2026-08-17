"""Five-seed confirmation of the locked geometry-preserving V9 bridge."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .benchmark import _bundle_centroids, _paired_improvement
from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .control_adapter import _donor_balanced_centroid
from .embedding import load_embedding, matched_rows
from .firewall import load_retargeted_bridge_selection_lock, validate_program_firewall
from .orthogonal_bridge import (
    _balanced_fit_positions, apply_scaled_orthogonal_bridge, fit_scaled_orthogonal_bridge,
)
from .orthogonal_bridge_outcome import _scope_centroids, evaluate_matched_scope


def _selected(selection: dict[str, Any]) -> dict[str, dict[str, Any]]:
    weight = float(selection["selected_global_reference_weight"])
    result = {}
    for source in selection["control_results"]:
        with Path(source["path"]).open() as handle:
            row = json.load(handle)
        if float(row["global_reference_weight"]) == weight:
            result[row["model_kind"]] = row
    return result


def _aligned_raw(raw_bundle, cells):
    raw_info, raw_latent, raw_cells = raw_bundle
    left, right = matched_rows(raw_cells, cells)
    if len(right) != len(cells):
        raise ContractError("seed confirmation raw PCA lacks candidate cells")
    order = np.argsort(right)
    left, right = left[order], right[order]
    if not np.array_equal(right, np.arange(len(cells))):
        raise ContractError("seed confirmation failed to recover candidate order")
    return np.asarray(raw_latent[left], dtype=np.float32)


def _reference_positions(reference_cells, cells):
    positions = np.flatnonzero(cells["strict_reference"].to_numpy(dtype=bool))
    left, right = matched_rows(reference_cells, cells.iloc[positions].reset_index(drop=True))
    if len(right) != len(positions):
        raise ContractError("seed confirmation reference roster differs")
    order = np.argsort(right)
    return left[order], positions[right[order]]


def _one_seed_model(
    config: dict[str, Any], seed: int, model_kind: str, control: dict[str, Any],
    selected_weight: float,
) -> dict[str, Any]:
    candidate_path = Path(control["sources"]["candidate_embedding"]).resolve()
    selected_bundle = load_embedding(candidate_path)
    _, selected_latent, cells = selected_bundle
    bridge_path = Path(control["sources"]["bridge_manifest"]).resolve()
    with bridge_path.open() as handle:
        bridge = json.load(handle)
    parent_manifest_path = Path(bridge["parent_manifest_realpath"])
    if sha256_path(parent_manifest_path) != bridge["parent_manifest_sha256"]:
        raise ContractError("seed confirmation parent bridge changed")
    with parent_manifest_path.open() as handle:
        parent = json.load(handle)
    raw_path = Path(parent["raw_pca_embedding_realpath"])
    architecture_path = Path(parent["architecture_embedding_realpath"])
    reference_path = Path(parent["reference_embedding_realpath"])
    for path, expected in (
        (raw_path, parent["raw_pca_embedding_sha256"]),
        (architecture_path, parent["architecture_embedding_sha256"]),
        (reference_path, parent["reference_embedding_sha256"]),
    ):
        if sha256_path(path) != expected:
            raise ContractError("seed confirmation bridge source changed")
    architecture = load_embedding(architecture_path)
    reference = load_embedding(reference_path)
    raw_bundle = load_embedding(raw_path)
    if not np.array_equal(cells["cell_id"].astype(str), architecture[2]["cell_id"].astype(str)):
        raise ContractError("seed confirmation architecture roster differs")
    raw = _aligned_raw(raw_bundle, cells)
    ref_left, ref_positions = _reference_positions(reference[2], cells)
    fit_positions = _balanced_fit_positions(cells, ref_positions, model_kind, config, seed=seed)
    ref_index = {cell: index for index, cell in enumerate(reference[2]["cell_id"].astype(str))}
    fit_reference = np.asarray(
        [ref_index[cell] for cell in cells.iloc[fit_positions]["cell_id"].astype(str)], dtype=np.int64
    )
    source_mean, target_mean, rotation, scale = fit_scaled_orthogonal_bridge(
        raw[fit_positions], np.asarray(reference[1])[fit_reference]
    )
    adapted = apply_scaled_orthogonal_bridge(raw, source_mean, target_mean, rotation, scale)
    adapted[ref_positions] = np.asarray(reference[1])[ref_left].astype(np.float32)
    reference_mask = cells["strict_reference"].to_numpy(dtype=bool)
    primary = cells["primary_query"].to_numpy(dtype=bool)
    controls = primary & cells["query_control"].to_numpy(dtype=bool) & ~reference_mask
    preparations = cells["preparation"].astype(str).to_numpy()
    datasets = cells["dataset"].astype(str).to_numpy()
    global_centroid, n_reference = _donor_balanced_centroid(adapted, cells, reference_mask)
    compatible_mask = reference_mask & (preparations == "unsorted")
    compatible_centroid, n_compatible = _donor_balanced_centroid(adapted, cells, compatible_mask)
    target = (1.0 - selected_weight) * compatible_centroid + selected_weight * global_centroid
    centered_errors = []
    for dataset in config["roles"]["primary_query_datasets"]:
        control_centroid, _ = _donor_balanced_centroid(adapted, cells, controls & (datasets == dataset))
        offset = target - control_centroid
        mask = primary & (datasets == dataset)
        before = np.asarray(adapted[mask], dtype=np.float64)
        adapted[mask] = (before + offset).astype(np.float32)
        after = np.asarray(adapted[mask], dtype=np.float64)
        centered_errors.append(float(np.max(np.abs(
            (after - after.mean(axis=0)) - (before - before.mean(axis=0))
        ), initial=0.0)))
    if not np.array_equal(adapted[reference_mask], np.asarray(reference[1])[ref_left].astype(np.float32)):
        raise ContractError("seed confirmation changed frozen reference coordinates")
    candidate_bundle = (selected_bundle[0], adapted, cells)
    lineage = None if model_kind == "all_lineage" else model_kind
    harmony = load_embedding(control["sources"]["harmony_embedding"])
    harmony_result = _paired_improvement(
        _bundle_centroids(candidate_bundle, lineage), _bundle_centroids(harmony, lineage),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + seed + 7101,
    )
    architecture_result = _paired_improvement(
        _bundle_centroids(candidate_bundle, lineage), _bundle_centroids(architecture, lineage),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + seed + 7201,
    )
    protocol_harmony = _paired_improvement(
        _bundle_centroids(candidate_bundle, lineage, preparation="unsorted"),
        _bundle_centroids(harmony, lineage, preparation="unsorted"),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + seed + 7301,
    )
    protocol_architecture = _paired_improvement(
        _bundle_centroids(candidate_bundle, lineage, preparation="unsorted"),
        _bundle_centroids(architecture, lineage, preparation="unsorted"),
        config["bootstrap"]["replicates"], config["bootstrap"]["seed"] + seed + 7401,
    )
    disease_values = []
    evaluated_lineages = config["lineages"] if model_kind == "all_lineage" else [model_kind]
    for evaluated_lineage in evaluated_lineages:
        for study in [None, *config["evaluation"]["powered_query_studies"]]:
            candidate_centroids, raw_centroids, _, is_control, _ = _scope_centroids(
                adapted, np.asarray(raw, dtype=np.float64), cells, evaluated_lineage, study
            )
            result = evaluate_matched_scope(candidate_centroids, raw_centroids, is_control)
            if result is not None:
                disease_values.extend([result["retention"], result["distance_spearman"]])
    if not disease_values:
        raise ContractError("seed confirmation has no testable disease scope")
    threshold = config["gates"]["minimum_control_alignment_improvement"]
    gates = {
        "harmony_improvement": harmony_result["improvement"] >= threshold and harmony_result["ci_low"] > 0,
        "architecture_improvement": architecture_result["improvement"] >= threshold and architecture_result["ci_low"] > 0,
        "protocol_harmony_no_material_worsening": protocol_harmony["improvement"] >= -config["gates"]["maximum_stratum_worsening"],
        "protocol_architecture_no_material_worsening": protocol_architecture["improvement"] >= -config["gates"]["maximum_stratum_worsening"],
        "reference_coordinates_frozen": True,
        "uniform_translation": max(centered_errors) <= 1e-5,
        "disease_geometry": min(disease_values) >= min(
            config["gates"]["per_study_disease_retention"],
            config["gates"]["within_study_distance_spearman"],
        ),
    }
    seed17_error = None
    if seed == config["screen"]["seed"]:
        seed17_error = float(np.max(np.abs(
            np.asarray(adapted, dtype=np.float64) - np.asarray(selected_latent, dtype=np.float64)
        ), initial=0.0))
        gates["reproduces_locked_seed17"] = seed17_error <= 1e-6
    return {
        "model_kind": model_kind, "seed": seed, "scale": scale,
        "n_fit_cells": int(len(fit_positions)), "n_reference_donors": n_reference,
        "n_compatible_reference_donors": n_compatible,
        "orthonormality_max_abs_error": float(np.max(np.abs(rotation.T @ rotation - np.eye(len(rotation))))),
        "maximum_centered_coordinate_error": max(centered_errors),
        "seed17_maximum_coordinate_error_vs_locked": seed17_error,
        "alignment": {
            "harmony": harmony_result, "architecture_surgery": architecture_result,
            "protocol_harmony": protocol_harmony, "protocol_architecture_surgery": protocol_architecture,
        },
        "minimum_disease_geometry_metric": min(disease_values), "gates": gates,
        "pass": bool(all(gates.values())),
        "sources": {
            "candidate_embedding": str(candidate_path), "bridge_manifest": str(bridge_path),
            "raw_pca_embedding": str(raw_path), "architecture_embedding": str(architecture_path),
            "reference_embedding": str(reference_path),
        },
    }


def confirm_orthogonal_bridge_seeds(
    config: dict[str, Any], selection_lock: str | Path, output: str | Path,
) -> dict[str, Any]:
    selection_path = Path(selection_lock).resolve()
    selection = load_retargeted_bridge_selection_lock(selection_path, config)
    validate_program_firewall(config)
    selected = _selected(selection)
    expected = {"all_lineage", *config["lineages"]}
    if set(selected) != expected:
        raise ContractError("seed confirmation lacks the exact six-model roster")
    results, seed_summaries = [], {}
    for seed in config["screen"]["confirmation_seeds"]:
        rows = [
            _one_seed_model(
                config, int(seed), kind, selected[kind],
                float(selection["selected_global_reference_weight"]),
            )
            for kind in ("all_lineage", *config["lineages"])
        ]
        results.extend(rows)
        all_lineage = next(row for row in rows if row["model_kind"] == "all_lineage")
        lineage_rows = [row for row in rows if row["model_kind"] in config["lineages"]]
        improved_lineages = sum(
            row["gates"]["harmony_improvement"] and row["gates"]["architecture_improvement"]
            for row in lineage_rows
        )
        preservation = all(
            row["gates"]["reference_coordinates_frozen"]
            and row["gates"]["uniform_translation"]
            and row["gates"]["disease_geometry"] for row in rows
        )
        protocol = all(
            row["gates"]["protocol_harmony_no_material_worsening"]
            and row["gates"]["protocol_architecture_no_material_worsening"] for row in rows
        )
        passed = bool(
            preservation and protocol
            and all_lineage["gates"]["harmony_improvement"]
            and all_lineage["gates"]["architecture_improvement"]
            and improved_lineages >= config["gates"]["minimum_improved_lineages"]
        )
        seed_summaries[str(seed)] = {
            "pass": passed, "preservation": preservation, "protocol": protocol,
            "improved_lineages": improved_lineages,
        }
    result = {
        "schema_version": "masld-cl-orthogonal-bridge-seed-confirmation-v9",
        "config_sha256": config["_config_sha256"],
        "selection_lock": str(selection_path), "selection_lock_file_sha256": sha256_path(selection_path),
        "selection_lock_sha256": selection["lock_sha256"],
        "selected_global_reference_weight": selection["selected_global_reference_weight"],
        "seeds": list(map(int, config["screen"]["confirmation_seeds"])),
        "seed_summaries": seed_summaries, "models": results,
        "all_five_seeds_pass": bool(all(row["pass"] for row in seed_summaries.values())),
        "gpu_training_seed_gate": "not_applicable_to_closed_form_bridge; seeds_resample_the_donor_balanced_reference_fit",
    }
    write_json_exclusive(output, result)
    return result

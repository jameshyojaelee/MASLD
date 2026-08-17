#!/usr/bin/env python3
"""Select a dual-head blend on strict-reference donors, then infer query labels."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import torch
from scvi.model import SCANVI

from masld_cl.config import load_config, write_json_exclusive
from masld_cl.contracts import ContractError, DeterministicGzipTextWriter, sha256_path
from masld_cl.embedding import load_embedding
from masld_cl.metrics import macro_f1
from masld_cl.training import _set_model_labels, set_all_seeds


def _digest(state: dict[str, torch.Tensor], keys: list[str]) -> str:
    result = hashlib.sha256()
    for key in keys:
        value = state[key].detach().cpu().contiguous()
        result.update(key.encode())
        result.update(str(value.dtype).encode())
        result.update(np.asarray(value.shape, dtype=np.int64).tobytes())
        result.update(value.numpy().tobytes())
    return result.hexdigest()


def _donor_scores(cells, predictions: np.ndarray) -> dict[str, float]:
    reference = cells.loc[cells["strict_reference"]].copy()
    reference["candidate_prediction"] = predictions[
        cells["strict_reference"].to_numpy(dtype=bool)
    ]
    scores = {
        str(donor): macro_f1(
            donor_cells["audit_cell_type"].astype(str).to_numpy(),
            donor_cells["candidate_prediction"].astype(str).to_numpy(),
        )
        for donor, donor_cells in reference.groupby("donor_id", sort=True)
    }
    if len(scores) != 7:
        raise ContractError(f"V17 requires exactly seven strict-reference donors, observed {len(scores)}")
    return scores


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--prepared", required=True)
    parser.add_argument("--reference-model", required=True)
    parser.add_argument("--adapted-model", required=True)
    parser.add_argument("--v14-embedding", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    policy_path = Path(args.policy).resolve()
    with policy_path.open() as handle:
        policy = json.load(handle)
    rule = policy.get("selection_rule", {})
    weights = rule.get("adapted_head_probability_weights")
    if (
        policy.get("schema_version")
        != "masld-cl-reference-only-dual-head-selection-policy-v17"
        or policy.get("config_sha256") != config["_config_sha256"]
        or weights != [0.0, 0.25, 0.5, 0.75, 1.0]
        or rule.get("query_labels_used") is not False
        or rule.get("query_predictions_generated_only_after_selection_lock") is not True
    ):
        raise ContractError("V17 policy identity or selection rule differs")

    v14_path = Path(args.v14_embedding).resolve()
    v14_info, _, v14_cells = load_embedding(v14_path)
    if v14_info.get("method") != "geometry_preserving_continual_adapter":
        raise ContractError("V17 V14 source differs")
    set_all_seeds(17)
    reference = SCANVI.load(Path(args.reference_model).resolve())
    adapted = SCANVI.load(Path(args.adapted_model).resolve())
    reference_state = reference.module.state_dict()
    adapted_state = adapted.module.state_dict()
    classifier = sorted(key for key in adapted_state if key.startswith("classifier."))
    reference_classifier = sorted(
        key for key in reference_state if key.startswith("classifier.")
    )
    if (
        not classifier
        or classifier != reference_classifier
        or any(adapted_state[key].shape != reference_state[key].shape for key in classifier)
    ):
        raise ContractError("V17 classifier state rosters differ")
    nonclassifier = sorted(set(adapted_state) - set(classifier))
    nonclassifier_before = _digest(adapted_state, nonclassifier)
    adapted_classifier_digest = _digest(adapted_state, classifier)
    reference_classifier_digest = _digest(reference_state, classifier)

    full = ad.read_h5ad(Path(args.prepared).resolve())
    cell_ids = v14_cells["cell_id"].astype(str).tolist()
    if (
        len(cell_ids) != len(set(cell_ids))
        or not set(cell_ids).issubset(set(full.obs_names.astype(str)))
    ):
        raise ContractError("V17 embedding and prepared cell rosters differ")
    mapped = full[cell_ids].copy()
    if not np.array_equal(mapped.obs_names.astype(str).to_numpy(), np.asarray(cell_ids)):
        raise ContractError("V17 prepared cell order differs")
    _set_model_labels(
        mapped,
        query_unknown=True,
        unlabeled=config["features"]["unlabeled_category"],
    )
    adapted_probabilities = adapted.predict(
        mapped, soft=True, batch_size=config["architecture"]["batch_size"]
    )
    with torch.no_grad():
        for key in classifier:
            adapted_state[key].copy_(reference_state[key])
    restored_state = adapted.module.state_dict()
    if (
        _digest(restored_state, nonclassifier) != nonclassifier_before
        or _digest(restored_state, classifier) != reference_classifier_digest
    ):
        raise ContractError("V17 restore changed nonclassifier state or failed")
    reference_probabilities = adapted.predict(
        mapped, soft=True, batch_size=config["architecture"]["batch_size"]
    )
    if (
        list(adapted_probabilities.columns) != list(reference_probabilities.columns)
        or adapted_probabilities.shape != reference_probabilities.shape
        or adapted_probabilities.shape[0] != len(v14_cells)
    ):
        raise ContractError("V17 probability outputs differ in class or cell roster")
    adapted_values = adapted_probabilities.to_numpy(dtype=np.float64)
    reference_values = reference_probabilities.to_numpy(dtype=np.float64)
    if not np.isfinite(adapted_values).all() or not np.isfinite(reference_values).all():
        raise ContractError("V17 probabilities contain non-finite values")
    classes = np.asarray(adapted_probabilities.columns.astype(str))

    weight_metrics = []
    for weight in weights:
        probabilities = weight * adapted_values + (1.0 - weight) * reference_values
        predictions = classes[np.argmax(probabilities, axis=1)]
        donor_scores = _donor_scores(v14_cells, predictions)
        donor_values = np.asarray(list(donor_scores.values()), dtype=float)
        weight_metrics.append(
            {
                "adapted_weight": weight,
                "reference_weight": 1.0 - weight,
                "donor_macro_f1": donor_scores,
                "mean_donor_macro_f1": float(donor_values.mean()),
                "standard_error": float(donor_values.std(ddof=1) / math.sqrt(len(donor_values))),
            }
        )
    best = max(weight_metrics, key=lambda item: item["mean_donor_macro_f1"])
    cutoff = best["mean_donor_macro_f1"] - best["standard_error"]
    eligible = [
        item for item in weight_metrics if item["mean_donor_macro_f1"] >= cutoff
    ]
    selected = min(eligible, key=lambda item: item["adapted_weight"])

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    lock = {
        "schema_version": "masld-cl-reference-only-head-selection-lock-v17",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "selection_cells": "strict_reference_only",
        "selection_unit": "biological_donor",
        "strict_reference_donors": sorted(
            v14_cells.loc[v14_cells["strict_reference"], "donor_id"].astype(str).unique()
        ),
        "n_strict_reference_donors": 7,
        "query_audit_labels_read": False,
        "case_stage_program_hero_gene_umap_cas13_read": False,
        "weight_metrics": weight_metrics,
        "best_mean_weight": best["adapted_weight"],
        "best_mean": best["mean_donor_macro_f1"],
        "best_standard_error": best["standard_error"],
        "one_standard_error_cutoff": cutoff,
        "eligible_adapted_weights": [item["adapted_weight"] for item in eligible],
        "selected_adapted_weight": selected["adapted_weight"],
        "selected_reference_weight": selected["reference_weight"],
        "tie_break": "smallest_adapted_weight",
        "query_predictions_generated_after_this_lock": True,
    }
    lock_path = output / "reference_selection_lock.json"
    write_json_exclusive(lock_path, lock)

    selected_weight = selected["adapted_weight"]
    probabilities = selected_weight * adapted_values + (1.0 - selected_weight) * reference_values
    predictions = classes[np.argmax(probabilities, axis=1)]
    reference_mask = v14_cells["strict_reference"].to_numpy(dtype=bool)
    predictions[reference_mask] = v14_cells.loc[
        reference_mask, "predicted_cell_type"
    ].astype(str).to_numpy()

    source_latent = v14_path.parent / v14_info["latent_file"]
    latent_path = output / v14_info["latent_file"]
    shutil.copy2(source_latent, latent_path)
    if sha256_path(latent_path) != sha256_path(source_latent):
        raise ContractError("V17 latent copy differs")
    output_cells = v14_cells.copy()
    output_cells["predicted_cell_type"] = predictions
    cells_path = output / v14_info["cells_file"]
    with DeterministicGzipTextWriter(cells_path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(output_cells.columns)
        writer.writerows(output_cells.itertuples(index=False, name=None))
    embedding = dict(v14_info)
    embedding.update(
        {
            "method": "geometry_preserving_continual_adapter_with_reference_selected_dual_head",
            "latent_sha256": sha256_path(latent_path),
            "cells_sha256": sha256_path(cells_path),
        }
    )
    write_json_exclusive(output / "embedding_manifest.json", embedding)
    manifest = {
        "schema_version": "masld-cl-reference-selected-dual-head-v17",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "selection_lock": {"path": str(lock_path.resolve()), "sha256": sha256_path(lock_path)},
        "query_labels_used_for_selection": False,
        "weights": {"adapted": selected_weight, "reference": 1.0 - selected_weight},
        "class_order": classes.tolist(),
        "classifier_keys": classifier,
        "adapted_classifier_digest": adapted_classifier_digest,
        "reference_classifier_digest": reference_classifier_digest,
        "nonclassifier_digest_before": nonclassifier_before,
        "nonclassifier_digest_after": _digest(restored_state, nonclassifier),
        "latent_bitwise_identical_to_v14": True,
        "v14_embedding": {"path": str(v14_path), "sha256": sha256_path(v14_path)},
        "prepared": {
            "path": str(Path(args.prepared).resolve()),
            "sha256": sha256_path(args.prepared),
        },
        "embedding": embedding,
        "development_query_label_results_previously_read": True,
        "promotion_requires_independent_held_study_confirmation": True,
    }
    write_json_exclusive(output / "reference_selected_dual_head_manifest.json", manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

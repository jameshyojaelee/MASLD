#!/usr/bin/env python3
"""Calibrate two scANVI heads using strict-reference donors only."""

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
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from masld_cl.config import load_config, write_json_exclusive
from masld_cl.contracts import ContractError, DeterministicGzipTextWriter, sha256_path
from masld_cl.embedding import load_embedding
from masld_cl.metrics import macro_f1, positive_class_f1
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


def _balanced_reference_positions(cells, *, cap: int, seed: int) -> np.ndarray:
    reference_positions = np.flatnonzero(cells["strict_reference"].to_numpy(dtype=bool))
    reference = cells.iloc[reference_positions].copy()
    rng = np.random.default_rng(seed)
    selected = []
    for _, group in reference.groupby(["donor_id", "audit_cell_type"], sort=True):
        positions = group.index.to_numpy(dtype=np.int64)
        if len(positions) > cap:
            positions = np.sort(rng.choice(positions, cap, replace=False))
        selected.extend(positions.tolist())
    selected = np.asarray(sorted(selected), dtype=np.int64)
    if not len(selected) or not np.all(cells.iloc[selected]["strict_reference"]):
        raise ContractError("V18 balanced sample is empty or contains query cells")
    return selected


def _fit(
    features, labels, *, c_value: float, solver: str, tolerance: float, max_iter: int,
):
    scaler = StandardScaler().fit(features)
    transformed = scaler.transform(features)
    model = LogisticRegression(
        penalty="l2",
        C=c_value,
        solver=solver,
        tol=tolerance,
        max_iter=max_iter,
        random_state=17,
    ).fit(transformed, labels)
    if int(np.max(model.n_iter_)) >= max_iter:
        raise ContractError(f"V18 calibrator did not converge for C={c_value}")
    return scaler, model


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
    rule = policy.get("calibration_rule", {})
    c_grid = rule.get("inverse_regularization_grid")
    policy_schema = policy.get("schema_version")
    versions = {
        "masld-cl-reference-only-stacked-calibration-policy-v18": "v18",
        "masld-cl-reference-only-stacked-calibration-policy-v19": "v19",
        "masld-cl-reference-only-maximin-ensemble-policy-v20": "v20",
    }
    run_version = versions.get(policy_schema)
    if (
        run_version is None
        or policy.get("config_sha256") != config["_config_sha256"]
        or c_grid != [0.01, 0.1, 1.0, 10.0]
        or rule.get("query_labels_used") is not False
        or rule.get("query_predictions_generated_only_after_selection_lock") is not True
    ):
        raise ContractError("V18 policy identity or calibration rule differs")
    clip = float(rule["probability_clip"])
    solver = str(rule["solver"])
    tolerance = float(rule["tolerance"])
    max_iter = int(rule["maximum_iterations"])
    expected_optimizer = {
        "v18": ("lbfgs", 500),
        "v19": ("newton-cholesky", 100),
        "v20": ("newton-cholesky", 100),
    }
    if (solver, max_iter) != expected_optimizer[run_version]:
        raise ContractError("calibration optimizer differs from the versioned policy")

    v14_path = Path(args.v14_embedding).resolve()
    v14_info, _, v14_cells = load_embedding(v14_path)
    if v14_info.get("method") != "geometry_preserving_continual_adapter":
        raise ContractError("V18 V14 source differs")
    set_all_seeds(17)
    reference = SCANVI.load(Path(args.reference_model).resolve())
    adapted = SCANVI.load(Path(args.adapted_model).resolve())
    reference_state = reference.module.state_dict()
    adapted_state = adapted.module.state_dict()
    classifier = sorted(key for key in adapted_state if key.startswith("classifier."))
    reference_classifier = sorted(key for key in reference_state if key.startswith("classifier."))
    if (
        not classifier
        or classifier != reference_classifier
        or any(adapted_state[key].shape != reference_state[key].shape for key in classifier)
    ):
        raise ContractError("V18 classifier state rosters differ")
    nonclassifier = sorted(set(adapted_state) - set(classifier))
    nonclassifier_before = _digest(adapted_state, nonclassifier)
    adapted_classifier_digest = _digest(adapted_state, classifier)
    reference_classifier_digest = _digest(reference_state, classifier)

    full = ad.read_h5ad(Path(args.prepared).resolve())
    cell_ids = v14_cells["cell_id"].astype(str).tolist()
    if len(cell_ids) != len(set(cell_ids)) or not set(cell_ids).issubset(set(full.obs_names.astype(str))):
        raise ContractError("V18 embedding and prepared cell rosters differ")
    mapped = full[cell_ids].copy()
    if not np.array_equal(mapped.obs_names.astype(str).to_numpy(), np.asarray(cell_ids)):
        raise ContractError("V18 prepared cell order differs")
    _set_model_labels(mapped, query_unknown=True, unlabeled=config["features"]["unlabeled_category"])
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
        raise ContractError("V18 restore changed nonclassifier state or failed")
    reference_probabilities = adapted.predict(
        mapped, soft=True, batch_size=config["architecture"]["batch_size"]
    )
    if (
        list(adapted_probabilities.columns) != list(reference_probabilities.columns)
        or adapted_probabilities.shape != reference_probabilities.shape
        or adapted_probabilities.shape[0] != len(v14_cells)
    ):
        raise ContractError("V18 probability outputs differ in class or cell roster")
    class_order = np.asarray(adapted_probabilities.columns.astype(str))
    adapted_values = adapted_probabilities.to_numpy(dtype=np.float64)
    reference_values = reference_probabilities.to_numpy(dtype=np.float64)
    features = np.concatenate(
        [np.log(np.clip(adapted_values, clip, 1.0)), np.log(np.clip(reference_values, clip, 1.0))],
        axis=1,
    )
    if not np.isfinite(features).all():
        raise ContractError("V18 stacked features contain non-finite values")

    sample_positions = _balanced_reference_positions(v14_cells, cap=2000, seed=17)
    sample_cells = v14_cells.iloc[sample_positions].reset_index(drop=True)
    sample_features = features[sample_positions]
    sample_labels = sample_cells["audit_cell_type"].astype(str).to_numpy()
    sample_donors = sample_cells["donor_id"].astype(str).to_numpy()
    donors = sorted(set(sample_donors))
    if len(donors) != 7 or set(sample_labels) != set(class_order):
        raise ContractError("V18 reference sample has the wrong donor or class roster")

    cv_metrics = []
    models_by_c = {}
    for c_value in c_grid:
        donor_scores = {}
        iterations = {}
        lineage_donor_scores = {lineage: {} for lineage in config["lineages"]}
        fold_models = {}
        for donor in donors:
            validation = sample_donors == donor
            training = ~validation
            scaler, model = _fit(
                sample_features[training], sample_labels[training],
                c_value=c_value, solver=solver, tolerance=tolerance, max_iter=max_iter,
            )
            predictions = model.predict(scaler.transform(sample_features[validation]))
            donor_scores[donor] = macro_f1(sample_labels[validation], predictions)
            iterations[donor] = int(np.max(model.n_iter_))
            fold_models[donor] = (scaler, model)
            for lineage in config["lineages"]:
                truth = sample_labels[validation] == lineage
                if truth.any():
                    lineage_donor_scores[lineage][donor] = positive_class_f1(
                        truth, predictions == lineage
                    )
        values = np.asarray(list(donor_scores.values()), dtype=float)
        metric = {
            "C": c_value,
            "donor_macro_f1": donor_scores,
            "iterations": iterations,
            "mean_donor_macro_f1": float(values.mean()),
            "standard_error": float(values.std(ddof=1) / math.sqrt(len(values))),
        }
        if run_version == "v20":
            lineage_summary = {}
            for lineage, scores in lineage_donor_scores.items():
                if not scores:
                    raise ContractError(f"V20 lineage is not reference-evaluable: {lineage}")
                lineage_values = np.asarray(list(scores.values()), dtype=float)
                lineage_summary[lineage] = {
                    "donor_positive_f1": scores,
                    "mean_donor_positive_f1": float(lineage_values.mean()),
                    "standard_error": float(
                        0.0 if len(lineage_values) == 1
                        else lineage_values.std(ddof=1) / math.sqrt(len(lineage_values))
                    ),
                }
            metric["major_lineage_scores"] = lineage_summary
            metric["maximin_major_lineage_f1"] = min(
                item["mean_donor_positive_f1"] for item in lineage_summary.values()
            )
        cv_metrics.append(metric)
        models_by_c[c_value] = fold_models
    if run_version == "v20":
        best = max(cv_metrics, key=lambda item: item["maximin_major_lineage_f1"])
        best_value = best["maximin_major_lineage_f1"]
        eligible = [item for item in cv_metrics if item["maximin_major_lineage_f1"] == best_value]
        selected = min(eligible, key=lambda item: item["C"])
        cutoff = None
        final_scaler = final_model = None
    else:
        best = max(cv_metrics, key=lambda item: item["mean_donor_macro_f1"])
        cutoff = best["mean_donor_macro_f1"] - best["standard_error"]
        eligible = [item for item in cv_metrics if item["mean_donor_macro_f1"] >= cutoff]
        selected = min(eligible, key=lambda item: item["C"])
        final_scaler, final_model = _fit(
            sample_features, sample_labels, c_value=selected["C"],
            solver=solver, tolerance=tolerance, max_iter=max_iter,
        )

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    calibrator_path = output / "calibrator_parameters.npz"
    if run_version == "v20":
        selected_models = models_by_c[selected["C"]]
        model_classes = [selected_models[donor][1].classes_ for donor in donors]
        if any(not np.array_equal(classes, model_classes[0]) for classes in model_classes):
            raise ContractError("V20 fold-model class rosters differ")
        np.savez_compressed(
            calibrator_path,
            donor_ids=np.asarray(donors),
            scaler_mean=np.stack([selected_models[donor][0].mean_ for donor in donors]),
            scaler_scale=np.stack([selected_models[donor][0].scale_ for donor in donors]),
            classes=model_classes[0],
            coefficients=np.stack([selected_models[donor][1].coef_ for donor in donors]),
            intercept=np.stack([selected_models[donor][1].intercept_ for donor in donors]),
        )
        final_iterations = {
            donor: int(np.max(selected_models[donor][1].n_iter_)) for donor in donors
        }
    else:
        np.savez_compressed(
            calibrator_path,
            scaler_mean=final_scaler.mean_,
            scaler_scale=final_scaler.scale_,
            classes=final_model.classes_,
            coefficients=final_model.coef_,
            intercept=final_model.intercept_,
        )
        final_iterations = int(np.max(final_model.n_iter_))
    lock = {
        "schema_version": f"masld-cl-reference-only-stacked-calibration-lock-{run_version}",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "selection_cells": "strict_reference_only",
        "selection_unit": "biological_donor",
        "strict_reference_donors": donors,
        "n_strict_reference_donors": 7,
        "balanced_sample_cells": int(len(sample_positions)),
        "query_audit_labels_read": False,
        "case_stage_program_hero_gene_umap_cas13_read": False,
        "cv_metrics": cv_metrics,
        "best_mean_C": best["C"],
        "selection_metric": (
            "maximin_major_lineage_f1" if run_version == "v20" else "mean_donor_macro_f1"
        ),
        "best_mean": (
            best["maximin_major_lineage_f1"]
            if run_version == "v20" else best["mean_donor_macro_f1"]
        ),
        "best_standard_error": None if run_version == "v20" else best["standard_error"],
        "one_standard_error_cutoff": cutoff,
        "eligible_C": [item["C"] for item in eligible],
        "selected_C": selected["C"],
        "tie_break": "smallest_C",
        "final_iterations": final_iterations,
        "calibrator_parameters": {
            "path": str(calibrator_path.resolve()), "sha256": sha256_path(calibrator_path)
        },
        "query_predictions_generated_after_this_lock": True,
    }
    lock_path = output / "reference_calibration_lock.json"
    write_json_exclusive(lock_path, lock)

    if run_version == "v20":
        averaged = np.zeros((len(features), len(model_classes[0])), dtype=np.float64)
        for donor in donors:
            scaler, model = selected_models[donor]
            averaged += model.predict_proba(scaler.transform(features))
        averaged /= len(donors)
        predictions = model_classes[0][np.argmax(averaged, axis=1)].astype(object)
    else:
        predictions = final_model.predict(final_scaler.transform(features)).astype(object)
    reference_mask = v14_cells["strict_reference"].to_numpy(dtype=bool)
    predictions[reference_mask] = v14_cells.loc[
        reference_mask, "predicted_cell_type"
    ].astype(str).to_numpy()
    source_latent = v14_path.parent / v14_info["latent_file"]
    latent_path = output / v14_info["latent_file"]
    shutil.copy2(source_latent, latent_path)
    if sha256_path(latent_path) != sha256_path(source_latent):
        raise ContractError("V18 latent copy differs")
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
            "method": (
                "geometry_preserving_continual_adapter_with_reference_maximin_head_ensemble"
                if run_version == "v20"
                else "geometry_preserving_continual_adapter_with_reference_calibrated_heads"
            ),
            "latent_sha256": sha256_path(latent_path),
            "cells_sha256": sha256_path(cells_path),
        }
    )
    write_json_exclusive(output / "embedding_manifest.json", embedding)
    manifest = {
        "schema_version": f"masld-cl-reference-calibrated-heads-{run_version}",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "calibration_lock": {"path": str(lock_path.resolve()), "sha256": sha256_path(lock_path)},
        "query_labels_used_for_selection_or_fit": False,
        "class_order": class_order.tolist(),
        "classifier_keys": classifier,
        "adapted_classifier_digest": adapted_classifier_digest,
        "reference_classifier_digest": reference_classifier_digest,
        "nonclassifier_digest_before": nonclassifier_before,
        "nonclassifier_digest_after": _digest(restored_state, nonclassifier),
        "latent_bitwise_identical_to_v14": True,
        "v14_embedding": {"path": str(v14_path), "sha256": sha256_path(v14_path)},
        "prepared": {"path": str(Path(args.prepared).resolve()), "sha256": sha256_path(args.prepared)},
        "embedding": embedding,
        "development_query_label_results_previously_read": True,
        "promotion_requires_independent_held_study_confirmation": True,
    }
    write_json_exclusive(output / "reference_calibrated_heads_manifest.json", manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run one donor-by-genomic-fold observed-multiome surface without ranking."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import time
from typing import Any, Mapping, Sequence

import h5py
import numpy as np

from masld_bench.adapters.observed_multiome_factorized import MODEL_IDS, export, fit, load, predict
from masld_bench.artifacts import freeze_tree, reject_symlink_components, verify_frozen_tree, write_json_exclusive
from masld_bench.observed_multiome_surface import (
    SurfacePreprocessState,
    apply_nonnegative_reduction,
    apply_sparse_reduction,
    apply_target_reduction,
    decode,
    export_preprocess,
    fit_nonnegative_reduction,
    fit_sparse_reduction,
    fit_target_reduction,
    load_preprocess,
    log_cpm,
    profile_deviance_skill,
    read_csr,
)


SCHEMA = "masld-bench-observed-multiome-factorized-surface-smoke-v1"
CAMPAIGN_ID = "gse296875_observed_multiome_factorized_surface_smoke_20260825"


class SurfaceSmokeError(ValueError):
    """Raised when the prespecified surface or its process separation differs."""


def _digest(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SurfaceSmokeError("JSON object required")
    return value


def _tree(root: Path, record: Mapping[str, Any], label: str) -> Path:
    lexical = reject_symlink_components(root / str(record.get("path", "")), label=label)
    path = lexical.resolve(strict=True)
    path.relative_to(root)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise SurfaceSmokeError(f"{label} drifted")
    return path


def validate_config(root: Path, config: Mapping[str, Any]) -> tuple[dict[str, Path], dict[str, Path]]:
    if config.get("schema_version") != SCHEMA or config.get("campaign_id") != CAMPAIGN_ID or config.get("dataset_id") != "gse296875" or config.get("stage") != "smoke":
        raise SurfaceSmokeError("surface smoke identity differs")
    if config.get("surface") != {"held_genomic_fold": 0, "held_donor_fold": 0, "seed": 20260825, "outer_training_donors": 28, "outer_training_donor_lineage_units": 140, "held_donors": 11, "held_donor_lineage_units": 55, "training_targets": 1000, "held_targets": 1000}:
        raise SurfaceSmokeError("surface census differs")
    if tuple(config.get("models", ())) != MODEL_IDS:
        raise SurfaceSmokeError("surface model roster differs")
    preprocessing = config.get("preprocessing")
    if not isinstance(preprocessing, dict) or preprocessing.get("components") != 4 or preprocessing.get("fit_on_outer_training_rows_only") is not True or preprocessing.get("fit_on_training_targets_only") is not True:
        raise SurfaceSmokeError("surface preprocessing differs")
    evaluation = config.get("evaluation")
    if not isinstance(evaluation, dict) or evaluation.get("ranking_authorized") is not False or evaluation.get("bootstrap_replicates") != 0 or evaluation.get("p_values") is not False:
        raise SurfaceSmokeError("surface evaluation differs")
    firewall = config.get("firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise SurfaceSmokeError("surface firewall is open")
    parents = config.get("parents")
    expected_parents = {"prediction_broker", "model_inputs", "training_labels", "target_features", "evaluator_outcomes", "verified_model_inputs", "verified_training_labels", "verified_target_features", "verified_evaluator_outcomes"}
    if not isinstance(parents, dict) or set(parents) != expected_parents:
        raise SurfaceSmokeError("surface parent roster differs")
    resolved: dict[str, Path] = {}
    for key, record in parents.items():
        if key == "evaluator_outcomes":
            path = reject_symlink_components(root / str(record.get("path", "")), label=key).resolve(strict=True)
            path.relative_to(root)
            resolved[key] = path
        else:
            resolved[key] = _tree(root, record, key)
    if _json(resolved["prediction_broker"] / "receipt.json").get("single_surface_single_seed_biological_smoke_authorized") is not True:
        raise SurfaceSmokeError("prediction broker did not authorize this gate")
    for key in ("verified_model_inputs", "verified_training_labels", "verified_target_features", "verified_evaluator_outcomes"):
        if _json(resolved[key] / "receipt.json").get("promotion_gate_passed") is not True:
            raise SurfaceSmokeError(f"{key} is not promoted")
    children = config.get("children")
    if not isinstance(children, dict) or set(children) != {"model_input", "training_label", "training_target_features", "held_target_features", "evaluator_outcome"}:
        raise SurfaceSmokeError("surface child roster differs")
    child_paths: dict[str, Path] = {}
    for key, record in children.items():
        if not isinstance(record, dict) or record.get("parent") not in resolved:
            raise SurfaceSmokeError("surface child parent differs")
        path = (resolved[record["parent"]] / str(record.get("path", ""))).resolve(strict=True)
        path.relative_to(resolved[record["parent"]])
        if key != "evaluator_outcome":
            verify_frozen_tree(path)
            if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
                raise SurfaceSmokeError(f"{key} child drifted")
        child_paths[key] = path
    runtime = config.get("runtime")
    if not isinstance(runtime, dict):
        raise SurfaceSmokeError("surface runtime differs")
    runtime_path = (root / str(runtime.get("path", ""))).resolve(strict=True)
    runtime_path.relative_to(root)
    if _digest(runtime_path / "ARTIFACTS.json") != runtime.get("artifacts_sha256") or _digest(runtime_path / "pip-freeze.txt") != runtime.get("pip_freeze_sha256"):
        raise SurfaceSmokeError("surface runtime drifted")
    return resolved, child_paths


def _axis_hash(values: Sequence[str]) -> str:
    return sha256("\n".join(values).encode("utf-8")).hexdigest()


def _surface_id(surface: Mapping[str, Any]) -> str:
    return f"genomic_{int(surface['held_genomic_fold'])}_donor_{int(surface['held_donor_fold'])}_seed_{int(surface['seed'])}"


def _verify_bound_tree(record: Mapping[str, Any], label: str) -> Path:
    path = Path(str(record.get("path", ""))).resolve(strict=True)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise SurfaceSmokeError(f"bound {label} drifted")
    return path


def _fit_stage(args: argparse.Namespace) -> None:
    binding = _json(args.binding)
    if set(binding) != {"model_input", "training_label", "training_target_features", "surface", "preprocessing", "fit"}:
        raise SurfaceSmokeError("fit binding roster differs")
    model_path = _verify_bound_tree(binding["model_input"], "model input")
    label_path = _verify_bound_tree(binding["training_label"], "training label")
    target_path = _verify_bound_tree(binding["training_target_features"], "training target features")
    surface = binding["surface"]
    held_donor_fold = int(surface["held_donor_fold"])
    with h5py.File(model_path / "model_input.h5", "r") as model, h5py.File(label_path / "training_labels.h5", "r") as labels, h5py.File(target_path / "target_features.h5", "r") as targets:
        model_rows = decode(model["rows/row_hash"][:], unique=True)
        model_lineages = decode(model["rows/lineage"][:], unique=False)
        donor_folds = np.asarray(model["rows/donor_fold"][:], dtype=int)
        train_indices = np.flatnonzero(donor_folds != held_donor_fold)
        train_rows = [model_rows[index] for index in train_indices]
        train_lineages = [model_lineages[index] for index in train_indices]
        label_rows = decode(labels["rows/row_hash"][:], unique=True)
        label_lineages = decode(labels["rows/lineage"][:], unique=False)
        if train_rows != label_rows or train_lineages != label_lineages or np.any(np.asarray(labels["rows/donor_fold"][:], dtype=int) == held_donor_fold):
            raise SurfaceSmokeError("outer-training row join differs")
        training_target_ids = decode(targets["training_target_hash"][:], unique=True)
        if training_target_ids != decode(labels["training_target_atac/training_target_hash"][:], unique=True):
            raise SurfaceSmokeError("training target join differs")
        rna_counts = read_csr(model["rna/counts_csr"])[train_indices]
        atac_counts = read_csr(model["observed_atac_input/counts_csr"])[train_indices]
        target_counts = read_csr(labels["training_target_atac/counts_csr"]).toarray()
        target_values = np.asarray(targets["features"][:], dtype=np.float64)
        gene_ids = decode(model["rna/ensembl_id"][:], unique=True)
        input_ids = decode(model["observed_atac_input/input_hash"][:], unique=True)
    expected = surface
    if len(train_rows) != expected["outer_training_donor_lineage_units"] or len(training_target_ids) != expected["training_targets"] or target_counts.shape != (len(train_rows), len(training_target_ids)):
        raise SurfaceSmokeError("outer-training surface shape differs")
    rna_normalized, rna_zero = log_cpm(rna_counts)
    atac_normalized, atac_zero = log_cpm(atac_counts)
    dimensions = int(binding["preprocessing"]["components"])
    seed = int(surface["seed"])
    iterations = int(binding["preprocessing"]["svd_iterations"])
    rna_embedding, rna_basis = fit_sparse_reduction(rna_normalized, components=dimensions, seed=seed, iterations=iterations)
    atac_embedding, atac_basis = fit_nonnegative_reduction(atac_normalized, components=dimensions, seed=seed, max_iter=int(binding["preprocessing"]["nmf_max_iter"]))
    target_embedding, target_mean, target_scale, target_basis = fit_target_reduction(target_values, components=dimensions)
    args.output.mkdir(parents=True, exist_ok=False)
    preprocess_state = SurfacePreprocessState(rna_basis, atac_basis, target_mean, target_scale, target_basis, dimensions, seed)
    export_preprocess(preprocess_state, args.output / "preprocess", {"training_rows": len(train_rows), "training_targets": len(training_target_ids), "gene_axis_sha256": _axis_hash(gene_ids), "observed_atac_axis_sha256": _axis_hash(input_ids), "training_target_axis_sha256": _axis_hash(training_target_ids), "rna_zero_library_rows": rna_zero, "observed_atac_zero_library_rows": atac_zero, "fit_on_outer_training_rows_only": True, "fit_on_training_targets_only": True})
    fit_config = binding["fit"]
    for model_id in MODEL_IDS:
        state = fit(model_id, row_ids=train_rows, strata=train_lineages, target_ids=training_target_ids, rna=rna_embedding, observed_atac=atac_embedding, sequence_features=target_embedding, target_covariates=target_embedding, target_counts=target_counts, seed=seed, ridge_alpha=float(fit_config["ridge_alpha"]), poisson_alpha=float(fit_config["poisson_alpha"]), poisson_max_iter=int(fit_config["poisson_max_iter"]), poisson_tolerance=float(fit_config["poisson_tolerance"]), maximum_expanded_pairs=int(fit_config["maximum_expanded_pairs"]))
        export(state, args.output / "models" / model_id)
    write_json_exclusive(args.output / "receipt.json", {"schema_version": "masld-bench-observed-multiome-factorized-surface-fit-receipt-v1", "held_genomic_fold": int(surface["held_genomic_fold"]), "held_donor_fold": held_donor_fold, "seed": seed, "models": list(MODEL_IDS), "training_rows": len(train_rows), "training_targets": len(training_target_ids), "rna_zero_library_rows": rna_zero, "observed_atac_zero_library_rows": atac_zero, "evaluator_artifact_bound": False, "held_target_features_bound": False, "model_fit": True, "benchmark_metric_calculated": False, "sealed_outcomes_read": False})
    digest = freeze_tree(args.output, metadata={"artifact_class": "gse296875_observed_multiome_factorized_surface_fit", "surface": _surface_id(surface), "models": len(MODEL_IDS), "evaluator_artifact_bound": False})
    print(json.dumps({"artifacts_sha256": digest}, sort_keys=True))


def _predict_stage(args: argparse.Namespace) -> None:
    binding = _json(args.binding)
    if set(binding) != {"model_input", "frozen_fit_state", "held_target_features", "surface", "fit"}:
        raise SurfaceSmokeError("predict binding roster differs")
    model_path = _verify_bound_tree(binding["model_input"], "model input")
    state_path = _verify_bound_tree(binding["frozen_fit_state"], "fit state")
    target_path = _verify_bound_tree(binding["held_target_features"], "held target features")
    held_donor_fold = int(binding["surface"]["held_donor_fold"])
    preprocess, metadata = load_preprocess(state_path / "preprocess")
    with h5py.File(model_path / "model_input.h5", "r") as model, h5py.File(target_path / "target_features.h5", "r") as targets:
        all_rows = decode(model["rows/row_hash"][:], unique=True)
        all_units = decode(model["rows/unit_hash"][:], unique=False)
        all_lineages = decode(model["rows/lineage"][:], unique=False)
        donor_folds = np.asarray(model["rows/donor_fold"][:], dtype=int)
        query_indices = np.flatnonzero(donor_folds == held_donor_fold)
        query_rows = [all_rows[index] for index in query_indices]
        query_units = [all_units[index] for index in query_indices]
        query_lineages = [all_lineages[index] for index in query_indices]
        gene_ids = decode(model["rna/ensembl_id"][:], unique=True)
        input_ids = decode(model["observed_atac_input/input_hash"][:], unique=True)
        if _axis_hash(gene_ids) != metadata["gene_axis_sha256"] or _axis_hash(input_ids) != metadata["observed_atac_axis_sha256"]:
            raise SurfaceSmokeError("query feature axis drifted")
        rna, rna_zero = log_cpm(read_csr(model["rna/counts_csr"])[query_indices])
        atac, atac_zero = log_cpm(read_csr(model["observed_atac_input/counts_csr"])[query_indices])
        held_target_ids = decode(targets["target_hash"][:], unique=True)
        held_values = np.asarray(targets["features"][:], dtype=np.float64)
    if len(query_rows) != binding["surface"]["held_donor_lineage_units"] or len(set(query_units)) != binding["surface"]["held_donors"] or len(held_target_ids) != binding["surface"]["held_targets"]:
        raise SurfaceSmokeError("query surface census differs")
    rna_embedding = apply_sparse_reduction(rna, preprocess.rna_components)
    atac_embedding = apply_nonnegative_reduction(atac, preprocess.atac_components)
    target_embedding = apply_target_reduction(held_values, preprocess.target_mean, preprocess.target_scale, preprocess.target_components)
    args.output.mkdir(parents=True, exist_ok=False)
    for model_id in MODEL_IDS:
        state = load(state_path / "models" / model_id)
        values = predict(state, row_ids=query_rows, strata=query_lineages, target_ids=held_target_ids, rna=rna_embedding, observed_atac=atac_embedding, sequence_features=target_embedding, target_covariates=target_embedding, maximum_expanded_pairs=int(binding["fit"]["maximum_expanded_pairs"]))
        with (args.output / f"{model_id}.npy").open("xb") as handle:
            np.save(handle, values, allow_pickle=False)
    write_json_exclusive(args.output / "identifiers.json", {"row_hash": query_rows, "unit_hash": query_units, "lineage": query_lineages, "target_hash": held_target_ids})
    write_json_exclusive(args.output / "receipt.json", {"schema_version": "masld-bench-observed-multiome-factorized-surface-prediction-receipt-v1", "held_genomic_fold": int(binding["surface"]["held_genomic_fold"]), "held_donor_fold": held_donor_fold, "models": list(MODEL_IDS), "query_rows": len(query_rows), "query_donors": len(set(query_units)), "held_targets": len(held_target_ids), "rna_zero_library_rows": rna_zero, "observed_atac_zero_library_rows": atac_zero, "training_label_artifact_bound": False, "evaluator_artifact_bound": False, "benchmark_metric_calculated": False, "sealed_outcomes_read": False})
    digest = freeze_tree(args.output, metadata={"artifact_class": "gse296875_observed_multiome_factorized_surface_predictions", "surface": _surface_id(binding["surface"]), "models": len(MODEL_IDS), "evaluator_artifact_bound": False})
    print(json.dumps({"artifacts_sha256": digest}, sort_keys=True))


def _commit_stage(args: argparse.Namespace) -> None:
    binding = _json(args.binding)
    if set(binding) != {"frozen_prediction_bundle"}:
        raise SurfaceSmokeError("commit binding roster differs")
    predictions = _verify_bound_tree(binding["frozen_prediction_bundle"], "prediction bundle")
    args.output.mkdir(parents=True, exist_ok=False)
    bundle_hash = _digest(predictions / "ARTIFACTS.json")
    write_json_exclusive(args.output / "commit.json", {"schema_version": "masld-bench-observed-multiome-factorized-surface-prediction-commit-v1", "prediction_bundle_path": str(predictions), "prediction_bundle_artifacts_sha256": bundle_hash, "committed_unix_time_ns": time.time_ns(), "evaluator_artifact_bound": False, "sealed_outcomes_read": False})
    digest = freeze_tree(args.output, metadata={"artifact_class": "gse296875_observed_multiome_factorized_surface_prediction_commit", "prediction_bundle_artifacts_sha256": bundle_hash, "evaluator_artifact_bound": False})
    print(json.dumps({"artifacts_sha256": digest}, sort_keys=True))


def _evaluate_stage(args: argparse.Namespace) -> None:
    binding = _json(args.binding)
    if set(binding) != {"prediction_commit", "evaluator_outcome", "surface", "evaluation"}:
        raise SurfaceSmokeError("evaluator binding roster differs")
    commit_path = _verify_bound_tree(binding["prediction_commit"], "prediction commit")
    evaluator_path = _verify_bound_tree(binding["evaluator_outcome"], "evaluator outcome")
    commit = _json(commit_path / "commit.json")
    predictions = Path(str(commit["prediction_bundle_path"])).resolve(strict=True)
    verify_frozen_tree(predictions)
    if _digest(predictions / "ARTIFACTS.json") != commit["prediction_bundle_artifacts_sha256"]:
        raise SurfaceSmokeError("committed prediction bundle drifted")
    identifiers = _json(predictions / "identifiers.json")
    held_donor_fold = int(binding["surface"]["held_donor_fold"])
    with h5py.File(evaluator_path / "evaluator_outcomes.h5", "r") as evaluator:
        rows = decode(evaluator["rows/row_hash"][:], unique=True)
        units = decode(evaluator["rows/unit_hash"][:], unique=False)
        lineages = decode(evaluator["rows/lineage"][:], unique=False)
        folds = np.asarray(evaluator["rows/donor_fold"][:], dtype=int)
        query = np.flatnonzero(folds == held_donor_fold)
        target_ids = decode(evaluator["target_atac/target_hash"][:], unique=True)
        counts = read_csr(evaluator["target_atac/counts_csr"])[query].toarray()
    query_rows = [rows[index] for index in query]
    query_units = [units[index] for index in query]
    query_lineages = [lineages[index] for index in query]
    if identifiers != {"row_hash": query_rows, "unit_hash": query_units, "lineage": query_lineages, "target_hash": target_ids}:
        raise SurfaceSmokeError("evaluator join differs")
    metrics = {}
    unit_skills = []
    unit_model_deviances = []
    unit_null_deviance = None
    pseudocount = float(binding["evaluation"]["pseudocount"])
    for model_id in MODEL_IDS:
        with (predictions / f"{model_id}.npy").open("rb") as handle:
            estimated = np.load(handle, allow_pickle=False)
        skill, model_deviance, null_deviance = profile_deviance_skill(counts, estimated, pseudocount=pseudocount)
        unit_skills.append(skill)
        unit_model_deviances.append(model_deviance)
        if unit_null_deviance is None:
            unit_null_deviance = null_deviance
        elif not np.array_equal(unit_null_deviance, null_deviance):
            raise SurfaceSmokeError("model-specific null deviance differs")
        valid = np.isfinite(skill)
        donor_values = [float(np.mean(skill[[index for index, value in enumerate(query_units) if value == unit and valid[index]]])) for unit in sorted(set(query_units)) if any(value == unit and valid[index] for index, value in enumerate(query_units))]
        lineage_values = {lineage: float(np.mean(skill[[index for index, value in enumerate(query_lineages) if value == lineage and valid[index]]])) for lineage in sorted(set(query_lineages)) if any(value == lineage and valid[index] for index, value in enumerate(query_lineages))}
        metrics[model_id] = {"donor_macro_profile_deviance_skill": float(np.mean(donor_values)), "donors_evaluable": len(donor_values), "donor_lineage_units_evaluable": int(valid.sum()), "lineage_profile_deviance_skill": lineage_values, "mean_model_deviance": float(np.mean(model_deviance[valid])), "mean_uniform_null_deviance": float(np.mean(null_deviance[valid]))}
    args.output.mkdir(parents=True, exist_ok=False)
    with (args.output / "unit_metrics.npz").open("xb") as handle:
        np.savez_compressed(
            handle,
            row_hash=np.asarray(query_rows),
            unit_hash=np.asarray(query_units),
            lineage=np.asarray(query_lineages),
            model_id=np.asarray(MODEL_IDS),
            profile_deviance_skill=np.column_stack(unit_skills),
            model_deviance=np.column_stack(unit_model_deviances),
            uniform_null_deviance=np.asarray(unit_null_deviance),
        )
    write_json_exclusive(args.output / "receipt.json", {"schema_version": "masld-bench-observed-multiome-factorized-surface-evaluator-receipt-v1", "held_genomic_fold": int(binding["surface"]["held_genomic_fold"]), "held_donor_fold": held_donor_fold, "biological_donors": len(set(query_units)), "donor_lineage_units": len(query_rows), "held_targets": len(target_ids), "endpoint": binding["evaluation"]["endpoint"], "metrics_unranked": metrics, "unit_metrics_written": True, "unit_metrics_shape": [len(query_rows), len(MODEL_IDS)], "prediction_commit_verified": True, "fit_state_bound": False, "model_input_bound": False, "training_label_bound": False, "bootstrap_replicates": 0, "p_values_calculated": False, "partial_ranking_authorized": False, "promotion_authorized": False, "sealed_outcomes_read": False})
    digest = freeze_tree(args.output, metadata={"artifact_class": "gse296875_observed_multiome_factorized_surface_evaluator", "surface": _surface_id(binding["surface"]), "biological_donors": len(set(query_units)), "partial_ranking_authorized": False})
    print(json.dumps({"artifacts_sha256": digest}, sort_keys=True))


def _invoke(output: Path, label: str, stage: str, binding: Path, stage_output: Path) -> None:
    completed = subprocess.run([sys.executable, str(Path(__file__).resolve(strict=True)), stage, "--binding", str(binding), "--output", str(stage_output)], check=False, capture_output=True, text=True)
    (output / f"{label}.stdout.log").write_text(completed.stdout, encoding="utf-8")
    (output / f"{label}.stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise SurfaceSmokeError(f"{stage} subprocess failed with exit code {completed.returncode}")


def _binding_record(path: Path) -> dict[str, str]:
    return {"path": str(path), "artifacts_sha256": _digest(path / "ARTIFACTS.json")}


def _orchestrate(args: argparse.Namespace) -> None:
    root = args.root.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    config_path.relative_to(root)
    config = _json(config_path)
    _, children = validate_config(root, config)
    output = reject_symlink_components(args.output, label="factorized surface smoke output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    bindings = output / "bindings"
    bindings.mkdir()
    fit_binding = bindings / "fit.json"
    write_json_exclusive(fit_binding, {"model_input": _binding_record(children["model_input"]), "training_label": _binding_record(children["training_label"]), "training_target_features": _binding_record(children["training_target_features"]), "surface": config["surface"], "preprocessing": config["preprocessing"], "fit": config["fit"]})
    state = output / "fit_state"
    _invoke(output, "fit", "fit", fit_binding, state)
    predict_binding = bindings / "predict.json"
    write_json_exclusive(predict_binding, {"model_input": _binding_record(children["model_input"]), "frozen_fit_state": _binding_record(state), "held_target_features": _binding_record(children["held_target_features"]), "surface": config["surface"], "fit": config["fit"]})
    predictions = output / "predictions"
    _invoke(output, "predict", "predict", predict_binding, predictions)
    commit_binding = bindings / "commit.json"
    write_json_exclusive(commit_binding, {"frozen_prediction_bundle": _binding_record(predictions)})
    commit = output / "prediction_commit"
    _invoke(output, "commit", "commit", commit_binding, commit)
    committed_ns = _json(commit / "commit.json")["committed_unix_time_ns"]
    evaluator_launch_ns = time.time_ns()
    if evaluator_launch_ns <= committed_ns:
        raise SurfaceSmokeError("evaluator launch preceded prediction commit")
    evaluator_binding = bindings / "evaluate.json"
    evaluator_record = config["children"]["evaluator_outcome"]
    write_json_exclusive(evaluator_binding, {"prediction_commit": _binding_record(commit), "evaluator_outcome": {"path": str(children["evaluator_outcome"]), "artifacts_sha256": evaluator_record["artifacts_sha256"]}, "surface": config["surface"], "evaluation": config["evaluation"]})
    evaluator = output / "evaluator"
    _invoke(output, "evaluate", "evaluate", evaluator_binding, evaluator)
    receipt = {"schema_version": "masld-bench-observed-multiome-factorized-surface-smoke-receipt-v1", "campaign_id": CAMPAIGN_ID, "dataset_id": "gse296875", "stage": "smoke", "surface": config["surface"], "models": list(MODEL_IDS), "separate_fit_process": True, "separate_predict_process": True, "prediction_committed_before_evaluator_launch": True, "evaluator_artifact_read_before_prediction_commit": False, "separate_evaluator_process": True, "biological_smoke_completed": True, "partial_ranking_authorized": False, "promotion_authorized": False, "sealed_outcomes_read": False, "next_gate": "independent_surface_smoke_verification_then_prespecified_full_rectangle_authorization"}
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [config_path, Path(__file__).resolve(strict=True), root / "src/masld_bench/observed_multiome_surface.py", root / "src/masld_bench/adapters/observed_multiome_factorized.py", root / "tests/unit/test_observed_multiome_surface.py", root / "tests/unit/test_observed_multiome_factorized_surface_smoke.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_factorized_surface_smoke", "campaign_id": CAMPAIGN_ID, "surface": "genomic_0_donor_0_seed_20260825", "partial_ranking_authorized": False, "sealed_outcomes_accessed": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    stages = parser.add_subparsers(dest="stage", required=True)
    orchestrate = stages.add_parser("orchestrate")
    orchestrate.add_argument("--root", type=Path, required=True)
    orchestrate.add_argument("--config", type=Path, required=True)
    orchestrate.add_argument("--output", type=Path, required=True)
    for stage in ("fit", "predict", "commit", "evaluate"):
        subparser = stages.add_parser(stage)
        subparser.add_argument("--binding", type=Path, required=True)
        subparser.add_argument("--output", type=Path, required=True)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    {"orchestrate": _orchestrate, "fit": _fit_stage, "predict": _predict_stage, "commit": _commit_stage, "evaluate": _evaluate_stage}[args.stage](args)


if __name__ == "__main__":
    main()

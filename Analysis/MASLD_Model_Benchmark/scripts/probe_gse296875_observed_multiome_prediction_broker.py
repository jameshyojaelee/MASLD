#!/usr/bin/env python3
"""Exercise fit, predict, commit, and evaluator process boundaries synthetically."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import time
from typing import Any, Mapping, Sequence

import numpy as np

from masld_bench.adapters.observed_multiome_factorized import MODEL_IDS, export, fit, load, predict
from masld_bench.artifacts import freeze_tree, reject_symlink_components, verify_frozen_tree, write_json_exclusive


SCHEMA = "masld-bench-observed-multiome-prediction-broker-fixture-v1"
CAMPAIGN_ID = "gse296875_observed_multiome_prediction_broker_fixture_20260825"


class PredictionBrokerError(ValueError):
    """Raised when process binding, state, prediction, or commit differs."""


def _digest(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise PredictionBrokerError("JSON object required")
    return value


def _tree(root: Path, record: Mapping[str, Any], label: str) -> Path:
    lexical = reject_symlink_components(root / str(record.get("path", "")), label=label)
    path = lexical.resolve(strict=True)
    path.relative_to(root)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise PredictionBrokerError(f"{label} drifted")
    return path


def validate_config(root: Path, config: Mapping[str, Any]) -> dict[str, Path]:
    if config.get("schema_version") != SCHEMA or config.get("campaign_id") != CAMPAIGN_ID or config.get("dataset_id") != "gse296875" or config.get("stage") != "smoke":
        raise PredictionBrokerError("broker fixture identity differs")
    if tuple(config.get("model_ids", ())) != MODEL_IDS or config.get("seed") != 20260825:
        raise PredictionBrokerError("broker model or seed roster differs")
    if config.get("synthetic_dimensions") != {"training_rows": 12, "query_rows": 6, "training_targets": 7, "held_targets": 5, "rna_features": 3, "observed_atac_features": 4, "target_features": 3}:
        raise PredictionBrokerError("broker synthetic dimensions differ")
    binding = config.get("process_bindings")
    if binding != {
        "fit": ["model_input", "training_labels", "training_target_features"],
        "predict": ["model_input", "frozen_fit_state", "held_target_features"],
        "commit": ["frozen_prediction_bundle"],
        "evaluate": ["prediction_commit", "synthetic_evaluator_outcomes"],
    }:
        raise PredictionBrokerError("broker process binding differs")
    firewall = config.get("firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise PredictionBrokerError("broker firewall is open")
    parents = config.get("parents")
    expected = {"supervised_contract", "verified_model_inputs", "verified_training_labels", "verified_target_features", "verified_evaluator_outcomes"}
    if not isinstance(parents, dict) or set(parents) != expected:
        raise PredictionBrokerError("broker parent roster differs")
    resolved = {key: _tree(root, parents[key], key) for key in parents}
    checks = {
        "verified_model_inputs": "promotion_gate_passed",
        "verified_training_labels": "promotion_gate_passed",
        "verified_target_features": "promotion_gate_passed",
        "verified_evaluator_outcomes": "promotion_gate_passed",
    }
    for key, field in checks.items():
        if _json(resolved[key] / "receipt.json").get(field) is not True:
            raise PredictionBrokerError(f"{key} is not promoted")
    contract = _json(resolved["supervised_contract"] / "receipt.json")
    if contract.get("contract_id") != "gse296875_observed_multiome_supervised_training_20260825" or contract.get("model_fit_may_bind_evaluator_outcomes") is not False:
        raise PredictionBrokerError("supervised contract differs")
    return resolved


def stage_binding_keys(stage: str) -> tuple[str, ...]:
    values = {
        "fit": ("model_input", "training_labels", "training_target_features"),
        "predict": ("model_input", "frozen_fit_state", "held_target_features"),
        "commit": ("frozen_prediction_bundle",),
        "evaluate": ("prediction_commit", "synthetic_evaluator_outcomes"),
    }
    try:
        return values[stage]
    except KeyError as error:
        raise PredictionBrokerError("unsupported broker stage") from error


def _write_npz(path: Path, **values: np.ndarray) -> None:
    with path.open("xb") as handle:
        np.savez_compressed(handle, **values)


def _load_npz(path: Path, expected: set[str]) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as archive:
        if set(archive.files) != expected:
            raise PredictionBrokerError(f"{path.name} array roster differs")
        return {name: np.asarray(archive[name]).copy() for name in archive.files}


def _strings(values: np.ndarray) -> list[str]:
    if values.ndim != 1:
        raise PredictionBrokerError("identifier axis differs")
    output = [str(value) for value in values.tolist()]
    if not output or any(not value for value in output) or len(set(output)) != len(output):
        raise PredictionBrokerError("identifier values differ")
    return output


def _labels(values: np.ndarray) -> list[str]:
    if values.ndim != 1:
        raise PredictionBrokerError("label axis differs")
    output = [str(value) for value in values.tolist()]
    if not output or any(not value for value in output):
        raise PredictionBrokerError("label values differ")
    return output


def _make_synthetic(config: Mapping[str, Any], output: Path) -> dict[str, Path]:
    dimensions = config["synthetic_dimensions"]
    generator = np.random.default_rng(int(config["seed"]))
    n_train = int(dimensions["training_rows"])
    n_query = int(dimensions["query_rows"])
    n_train_targets = int(dimensions["training_targets"])
    n_held_targets = int(dimensions["held_targets"])
    rna = generator.normal(size=(n_train + n_query, int(dimensions["rna_features"])))
    atac = generator.poisson(3.0, size=(n_train + n_query, int(dimensions["observed_atac_features"]))).astype(float)
    target = generator.normal(size=(n_train_targets + n_held_targets, int(dimensions["target_features"])))
    train_rate = np.exp(np.clip(0.3 + rna[:n_train] @ target[:n_train_targets].T * 0.1, -2.0, 2.0))
    held_rate = np.exp(np.clip(0.3 + rna[n_train:] @ target[n_train_targets:].T * 0.1, -2.0, 2.0))
    source = output / "synthetic_sources"
    source.mkdir()
    paths = {
        "model_input": source / "model_input.npz",
        "training_labels": source / "training_labels.npz",
        "training_target_features": source / "training_target_features.npz",
        "held_target_features": source / "held_target_features.npz",
        "synthetic_evaluator_outcomes": source / "synthetic_evaluator_outcomes.npz",
    }
    _write_npz(
        paths["model_input"],
        training_row_ids=np.asarray([f"train-row-{index:03d}" for index in range(n_train)]),
        query_row_ids=np.asarray([f"query-row-{index:03d}" for index in range(n_query)]),
        training_strata=np.asarray([f"lineage-{index % 2}" for index in range(n_train)]),
        query_strata=np.asarray([f"lineage-{index % 2}" for index in range(n_query)]),
        training_rna=rna[:n_train],
        query_rna=rna[n_train:],
        training_observed_atac=atac[:n_train],
        query_observed_atac=atac[n_train:],
    )
    _write_npz(
        paths["training_labels"],
        training_row_ids=np.asarray([f"train-row-{index:03d}" for index in range(n_train)]),
        training_target_ids=np.asarray([f"train-target-{index:03d}" for index in range(n_train_targets)]),
        training_target_counts=generator.poisson(train_rate).astype(float),
    )
    _write_npz(
        paths["training_target_features"],
        training_target_ids=np.asarray([f"train-target-{index:03d}" for index in range(n_train_targets)]),
        sequence_features=target[:n_train_targets],
        target_covariates=target[:n_train_targets] * 0.5,
    )
    _write_npz(
        paths["held_target_features"],
        held_target_ids=np.asarray([f"held-target-{index:03d}" for index in range(n_held_targets)]),
        sequence_features=target[n_train_targets:],
        target_covariates=target[n_train_targets:] * 0.5,
    )
    _write_npz(
        paths["synthetic_evaluator_outcomes"],
        query_row_ids=np.asarray([f"query-row-{index:03d}" for index in range(n_query)]),
        held_target_ids=np.asarray([f"held-target-{index:03d}" for index in range(n_held_targets)]),
        held_target_counts=generator.poisson(held_rate).astype(float),
    )
    return paths


def _fit_stage(args: argparse.Namespace) -> None:
    model = _load_npz(args.model_input, {"training_row_ids", "query_row_ids", "training_strata", "query_strata", "training_rna", "query_rna", "training_observed_atac", "query_observed_atac"})
    labels = _load_npz(args.training_labels, {"training_row_ids", "training_target_ids", "training_target_counts"})
    targets = _load_npz(args.training_target_features, {"training_target_ids", "sequence_features", "target_covariates"})
    rows = _strings(model["training_row_ids"])
    target_ids = _strings(targets["training_target_ids"])
    if rows != _strings(labels["training_row_ids"]) or target_ids != _strings(labels["training_target_ids"]):
        raise PredictionBrokerError("fit join differs")
    args.output.mkdir(parents=True, exist_ok=False)
    for model_id in MODEL_IDS:
        state = fit(
            model_id,
            row_ids=rows,
            strata=_labels(model["training_strata"]),
            target_ids=target_ids,
            rna=model["training_rna"],
            observed_atac=model["training_observed_atac"],
            sequence_features=targets["sequence_features"],
            target_covariates=targets["target_covariates"],
            target_counts=labels["training_target_counts"],
            seed=args.seed,
        )
        export(state, args.output / model_id)
    write_json_exclusive(args.output / "receipt.json", {"schema_version": "masld-bench-observed-multiome-broker-fit-receipt-v1", "stage": "fit", "models": list(MODEL_IDS), "binding_keys": list(stage_binding_keys("fit")), "training_rows": len(rows), "training_targets": len(target_ids), "evaluator_outcomes_bound": False, "held_target_features_bound": False, "biological_data_read": False, "sealed_outcomes_read": False})
    digest = freeze_tree(args.output, metadata={"artifact_class": "observed_multiome_synthetic_broker_fit_state", "models": len(MODEL_IDS), "evaluator_outcomes_bound": False})
    print(json.dumps({"artifacts_sha256": digest}, sort_keys=True))


def _predict_stage(args: argparse.Namespace) -> None:
    model = _load_npz(args.model_input, {"training_row_ids", "query_row_ids", "training_strata", "query_strata", "training_rna", "query_rna", "training_observed_atac", "query_observed_atac"})
    held = _load_npz(args.held_target_features, {"held_target_ids", "sequence_features", "target_covariates"})
    verify_frozen_tree(args.frozen_fit_state)
    rows = _strings(model["query_row_ids"])
    targets = _strings(held["held_target_ids"])
    args.output.mkdir(parents=True, exist_ok=False)
    for model_id in MODEL_IDS:
        state = load(args.frozen_fit_state / model_id)
        values = predict(
            state,
            row_ids=rows,
            strata=_labels(model["query_strata"]),
            target_ids=targets,
            rna=model["query_rna"],
            observed_atac=model["query_observed_atac"],
            sequence_features=held["sequence_features"],
            target_covariates=held["target_covariates"],
        )
        with (args.output / f"{model_id}.npy").open("xb") as handle:
            np.save(handle, values, allow_pickle=False)
    write_json_exclusive(args.output / "identifiers.json", {"query_row_ids": rows, "held_target_ids": targets})
    write_json_exclusive(args.output / "receipt.json", {"schema_version": "masld-bench-observed-multiome-broker-predict-receipt-v1", "stage": "predict", "models": list(MODEL_IDS), "binding_keys": list(stage_binding_keys("predict")), "query_rows": len(rows), "held_targets": len(targets), "training_labels_bound": False, "evaluator_outcomes_bound": False, "biological_data_read": False, "benchmark_metric_calculated": False, "sealed_outcomes_read": False})
    digest = freeze_tree(args.output, metadata={"artifact_class": "observed_multiome_synthetic_broker_prediction_bundle", "models": len(MODEL_IDS), "training_labels_bound": False, "evaluator_outcomes_bound": False})
    print(json.dumps({"artifacts_sha256": digest}, sort_keys=True))


def _commit_stage(args: argparse.Namespace) -> None:
    verify_frozen_tree(args.frozen_prediction_bundle)
    args.output.mkdir(parents=True, exist_ok=False)
    bundle_hash = _digest(args.frozen_prediction_bundle / "ARTIFACTS.json")
    write_json_exclusive(args.output / "commit.json", {"schema_version": "masld-bench-observed-multiome-broker-prediction-commit-v1", "stage": "commit", "binding_keys": list(stage_binding_keys("commit")), "prediction_bundle_path": str(args.frozen_prediction_bundle.resolve(strict=True)), "prediction_bundle_artifacts_sha256": bundle_hash, "committed_unix_time_ns": time.time_ns(), "evaluator_outcomes_bound": False, "sealed_outcomes_read": False})
    digest = freeze_tree(args.output, metadata={"artifact_class": "observed_multiome_synthetic_broker_prediction_commit", "prediction_bundle_artifacts_sha256": bundle_hash, "evaluator_outcomes_bound": False})
    print(json.dumps({"artifacts_sha256": digest}, sort_keys=True))


def _evaluate_stage(args: argparse.Namespace) -> None:
    verify_frozen_tree(args.prediction_commit)
    commit = _json(args.prediction_commit / "commit.json")
    bundle = Path(str(commit["prediction_bundle_path"])).resolve(strict=True)
    verify_frozen_tree(bundle)
    if _digest(bundle / "ARTIFACTS.json") != commit.get("prediction_bundle_artifacts_sha256"):
        raise PredictionBrokerError("committed prediction bundle drifted")
    outcomes = _load_npz(args.synthetic_evaluator_outcomes, {"query_row_ids", "held_target_ids", "held_target_counts"})
    identifiers = _json(bundle / "identifiers.json")
    if identifiers["query_row_ids"] != _strings(outcomes["query_row_ids"]) or identifiers["held_target_ids"] != _strings(outcomes["held_target_ids"]):
        raise PredictionBrokerError("evaluator join differs")
    scores = {}
    for model_id in MODEL_IDS:
        with (bundle / f"{model_id}.npy").open("rb") as handle:
            prediction = np.load(handle, allow_pickle=False)
        if prediction.shape != outcomes["held_target_counts"].shape or not np.isfinite(prediction).all() or np.any(prediction < 0):
            raise PredictionBrokerError("committed prediction differs")
        scores[model_id] = float(np.sqrt(np.mean((prediction - outcomes["held_target_counts"]) ** 2)))
    args.output.mkdir(parents=True, exist_ok=False)
    write_json_exclusive(args.output / "receipt.json", {"schema_version": "masld-bench-observed-multiome-broker-evaluator-receipt-v1", "stage": "evaluate", "binding_keys": list(stage_binding_keys("evaluate")), "models": list(MODEL_IDS), "synthetic_rmse": scores, "prediction_commit_verified": True, "fit_state_bound": False, "model_input_bound": False, "training_labels_bound": False, "synthetic_outcomes_only": True, "biological_data_read": False, "sealed_outcomes_read": False})
    digest = freeze_tree(args.output, metadata={"artifact_class": "observed_multiome_synthetic_broker_evaluator", "prediction_commit_verified": True, "biological_data_read": False})
    print(json.dumps({"artifacts_sha256": digest}, sort_keys=True))


def _invoke(output: Path, label: str, arguments: Sequence[str]) -> None:
    completed = subprocess.run([sys.executable, str(Path(__file__).resolve(strict=True)), *arguments], check=True, capture_output=True, text=True)
    (output / f"{label}.stdout.log").write_text(completed.stdout, encoding="utf-8")
    (output / f"{label}.stderr.log").write_text(completed.stderr, encoding="utf-8")


def _orchestrate(args: argparse.Namespace) -> None:
    root = args.root.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    config_path.relative_to(root)
    config = _json(config_path)
    validate_config(root, config)
    output = reject_symlink_components(args.output, label="prediction broker fixture output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    paths = _make_synthetic(config, output)
    state = output / "fit_state"
    predictions = output / "predictions"
    commit = output / "prediction_commit"
    evaluator = output / "evaluator"
    _invoke(output, "fit", ["fit", "--model-input", str(paths["model_input"]), "--training-labels", str(paths["training_labels"]), "--training-target-features", str(paths["training_target_features"]), "--seed", str(config["seed"]), "--output", str(state)])
    _invoke(output, "predict", ["predict", "--model-input", str(paths["model_input"]), "--frozen-fit-state", str(state), "--held-target-features", str(paths["held_target_features"]), "--output", str(predictions)])
    _invoke(output, "commit", ["commit", "--frozen-prediction-bundle", str(predictions), "--output", str(commit)])
    committed_ns = _json(commit / "commit.json")["committed_unix_time_ns"]
    evaluator_launch_ns = time.time_ns()
    if evaluator_launch_ns <= committed_ns:
        raise PredictionBrokerError("evaluator launch did not follow prediction commit")
    _invoke(output, "evaluate", ["evaluate", "--prediction-commit", str(commit), "--synthetic-evaluator-outcomes", str(paths["synthetic_evaluator_outcomes"]), "--output", str(evaluator)])
    receipt = {"schema_version": "masld-bench-observed-multiome-prediction-broker-fixture-receipt-v1", "campaign_id": CAMPAIGN_ID, "dataset_id": "gse296875", "stage": "smoke", "models": list(MODEL_IDS), "separate_fit_process": True, "separate_predict_process": True, "prediction_committed_before_evaluator_launch": True, "separate_evaluator_process": True, "fit_bound_evaluator_outcomes": False, "predict_bound_training_labels": False, "predict_bound_evaluator_outcomes": False, "evaluator_bound_fit_state": False, "evaluator_bound_model_input": False, "biological_data_read": False, "sealed_outcomes_read": False, "broker_fixture_passed": True, "single_surface_single_seed_biological_smoke_authorized": True, "next_gate": "single_surface_single_seed_biological_smoke"}
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [config_path, Path(__file__).resolve(strict=True), root / "src/masld_bench/adapters/observed_multiome_factorized.py", root / "tests/unit/test_observed_multiome_prediction_broker.py", root / "tests/unit/test_observed_multiome_factorized.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_prediction_broker_fixture", "campaign_id": CAMPAIGN_ID, "broker_fixture_passed": True, "biological_data_read": False, "sealed_outcomes_accessed": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    stages = parser.add_subparsers(dest="stage", required=True)
    orchestrate = stages.add_parser("orchestrate")
    orchestrate.add_argument("--root", type=Path, required=True)
    orchestrate.add_argument("--config", type=Path, required=True)
    orchestrate.add_argument("--output", type=Path, required=True)
    fit_parser = stages.add_parser("fit")
    fit_parser.add_argument("--model-input", type=Path, required=True)
    fit_parser.add_argument("--training-labels", type=Path, required=True)
    fit_parser.add_argument("--training-target-features", type=Path, required=True)
    fit_parser.add_argument("--seed", type=int, required=True)
    fit_parser.add_argument("--output", type=Path, required=True)
    predict_parser = stages.add_parser("predict")
    predict_parser.add_argument("--model-input", type=Path, required=True)
    predict_parser.add_argument("--frozen-fit-state", type=Path, required=True)
    predict_parser.add_argument("--held-target-features", type=Path, required=True)
    predict_parser.add_argument("--output", type=Path, required=True)
    commit_parser = stages.add_parser("commit")
    commit_parser.add_argument("--frozen-prediction-bundle", type=Path, required=True)
    commit_parser.add_argument("--output", type=Path, required=True)
    evaluate_parser = stages.add_parser("evaluate")
    evaluate_parser.add_argument("--prediction-commit", type=Path, required=True)
    evaluate_parser.add_argument("--synthetic-evaluator-outcomes", type=Path, required=True)
    evaluate_parser.add_argument("--output", type=Path, required=True)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    {"orchestrate": _orchestrate, "fit": _fit_stage, "predict": _predict_stage, "commit": _commit_stage, "evaluate": _evaluate_stage}[args.stage](args)


if __name__ == "__main__":
    main()

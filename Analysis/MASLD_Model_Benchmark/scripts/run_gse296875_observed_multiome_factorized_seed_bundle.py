#!/usr/bin/env python3
"""Run one complete 25-surface seed bundle under the frozen requirements."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import time
from typing import Any, Mapping, Sequence

from masld_bench.artifacts import freeze_tree, reject_symlink_components, verify_frozen_tree, write_json_exclusive


CONTRACT_SCHEMA = "masld-bench-observed-multiome-factorized-full-rectangle-contract-v1"
CONTRACT_ID = "gse296875_observed_multiome_factorized_full_rectangle_20260825"
CONTRACT_ARTIFACT_SHA256 = "cb9db39bb00874113b4c1f8ef797e39e79e9b94e96fda855ddb6c09768ad8f2c"
SEEDS = (20260824, 20260825, 20260826, 20260827, 20260828)
MODEL_IDS = ("masked_modality", "observed_atac_glm", "observed_atac_only", "rna_only", "shuffled_modality")


class SeedBundleError(ValueError):
    """Raised when a seed bundle, surface, or process output file differs."""


def _digest(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SeedBundleError("JSON object required")
    return value


def _tree(root: Path, record: Mapping[str, Any], label: str, *, verify: bool = True) -> Path:
    path = reject_symlink_components(root / str(record.get("path", "")), label=label).resolve(strict=True)
    path.relative_to(root)
    if verify:
        verify_frozen_tree(path)
        if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
            raise SeedBundleError(f"{label} drifted")
    return path


def _record(path: Path, expected_sha256: str | None = None) -> dict[str, str]:
    observed = _digest(path / "ARTIFACTS.json")
    if expected_sha256 is not None and observed != expected_sha256:
        raise SeedBundleError("child artifact drifted")
    return {"path": str(path), "artifacts_sha256": observed}


def _index(rows: Sequence[Mapping[str, Any]], keys: Sequence[str]) -> dict[tuple[Any, ...], Mapping[str, Any]]:
    output = {tuple(row[key] for key in keys): row for row in rows}
    if len(output) != len(rows):
        raise SeedBundleError("child artifact keys are duplicated")
    return output


def validate(root: Path, config: Mapping[str, Any], contract_path: Path, seed: int) -> dict[str, Any]:
    if seed not in SEEDS or config.get("schema_version") != CONTRACT_SCHEMA or config.get("contract_id") != CONTRACT_ID:
        raise SeedBundleError("seed bundle identity differs")
    verify_frozen_tree(contract_path)
    if _digest(contract_path / "ARTIFACTS.json") != CONTRACT_ARTIFACT_SHA256 or _json(contract_path / "receipt.json").get("full_rectangle_execution_authorized") is not True:
        raise SeedBundleError("full rectangle execution is not authorized")
    sources = config.get("source_bindings")
    if not isinstance(sources, dict):
        raise SeedBundleError("seed bundle source bindings differ")
    for record in sources.values():
        path = (root / str(record["path"])).resolve(strict=True)
        path.relative_to(root)
        if _digest(path) != record["sha256"]:
            raise SeedBundleError("seed bundle source drifted")
    parents = config.get("parents")
    if not isinstance(parents, dict):
        raise SeedBundleError("seed bundle parents differ")
    resolved = {key: _tree(root, record, key, verify=key != "evaluator_outcomes") for key, record in parents.items()}
    if _json(resolved["verified_evaluator_outcomes"] / "receipt.json").get("promotion_gate_passed") is not True:
        raise SeedBundleError("evaluator verifier differs")
    model = _json(resolved["model_inputs"] / "receipt.json")
    labels = _json(resolved["training_labels"] / "receipt.json")
    targets = _json(resolved["target_features"] / "receipt.json")
    evaluator = _json(resolved["evaluator_outcomes"] / "receipt.json")
    return {
        "resolved": resolved,
        "model": _index(model["child_artifacts"], ("held_genomic_fold",)),
        "labels": _index(labels["child_artifacts"], ("held_genomic_fold", "held_donor_fold")),
        "targets": _index(targets["children"], ("role", "held_genomic_fold")),
        "evaluator": _index(evaluator["child_artifacts"], ("held_genomic_fold",)),
    }


def planned_surfaces(seed: int) -> list[dict[str, int]]:
    if seed not in SEEDS:
        raise SeedBundleError("seed is outside the frozen roster")
    return [{"held_genomic_fold": genomic, "held_donor_fold": donor, "seed": seed} for genomic in range(5) for donor in range(5)]


def _invoke(root: Path, surface_root: Path, label: str, stage: str, binding: Path, output: Path) -> None:
    runner = root / "scripts/run_gse296875_observed_multiome_factorized_surface_smoke.py"
    completed = subprocess.run([sys.executable, str(runner), stage, "--binding", str(binding), "--output", str(output)], check=False, capture_output=True, text=True)
    (surface_root / f"{label}.stdout.log").write_text(completed.stdout, encoding="utf-8")
    (surface_root / f"{label}.stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise SeedBundleError(f"surface {stage} failed with exit code {completed.returncode}")


def _child(parent: Path, record: Mapping[str, Any], *, verify: bool) -> Path:
    path = (parent / str(record["path"])).resolve(strict=True)
    path.relative_to(parent)
    if verify:
        verify_frozen_tree(path)
        if _digest(path / "ARTIFACTS.json") != record["artifacts_sha256"]:
            raise SeedBundleError("surface child drifted")
    return path


def _run_surface(root: Path, config: Mapping[str, Any], inputs: Mapping[str, Any], seed: int, genomic: int, donor: int, output: Path) -> dict[str, Any]:
    resolved = inputs["resolved"]
    model_record = inputs["model"][(genomic,)]
    label_record = inputs["labels"][(genomic, donor)]
    training_target_record = inputs["targets"][("training", genomic)]
    held_target_record = inputs["targets"][("held", genomic)]
    evaluator_record = inputs["evaluator"][(genomic,)]
    model = _child(resolved["model_inputs"], model_record, verify=True)
    labels = _child(resolved["training_labels"], label_record, verify=True)
    training_targets = _child(resolved["target_features"], training_target_record, verify=True)
    held_targets = _child(resolved["target_features"], held_target_record, verify=True)
    evaluator = _child(resolved["evaluator_outcomes"], evaluator_record, verify=False)
    label_receipt = _json(labels / "receipt.json")
    surface = {
        "held_genomic_fold": genomic,
        "held_donor_fold": donor,
        "seed": seed,
        "outer_training_donors": int(label_receipt["outer_training_donors"]),
        "outer_training_donor_lineage_units": int(label_receipt["outer_training_donor_lineage_units"]),
        "held_donors": 39 - int(label_receipt["outer_training_donors"]),
        "held_donor_lineage_units": 195 - int(label_receipt["outer_training_donor_lineage_units"]),
        "training_targets": 1000,
        "held_targets": 1000,
    }
    output.mkdir(parents=True, exist_ok=False)
    bindings = output / "bindings"
    bindings.mkdir()
    fit_binding = bindings / "fit.json"
    write_json_exclusive(fit_binding, {"model_input": _record(model, model_record["artifacts_sha256"]), "training_label": _record(labels, label_record["artifacts_sha256"]), "training_target_features": _record(training_targets, training_target_record["artifacts_sha256"]), "surface": surface, "preprocessing": config["preprocessing"], "fit": config["fit"]})
    state = output / "fit_state"
    _invoke(root, output, "fit", "fit", fit_binding, state)
    predict_binding = bindings / "predict.json"
    write_json_exclusive(predict_binding, {"model_input": _record(model), "frozen_fit_state": _record(state), "held_target_features": _record(held_targets, held_target_record["artifacts_sha256"]), "surface": surface, "fit": config["fit"]})
    predictions = output / "predictions"
    _invoke(root, output, "predict", "predict", predict_binding, predictions)
    commit_binding = bindings / "commit.json"
    write_json_exclusive(commit_binding, {"frozen_prediction_bundle": _record(predictions)})
    commit = output / "prediction_commit"
    _invoke(root, output, "commit", "commit", commit_binding, commit)
    committed_ns = int(_json(commit / "commit.json")["committed_unix_time_ns"])
    if time.time_ns() <= committed_ns:
        raise SeedBundleError("surface evaluator launch preceded commit")
    evaluate_binding = bindings / "evaluate.json"
    write_json_exclusive(evaluate_binding, {"prediction_commit": _record(commit), "evaluator_outcome": {"path": str(evaluator), "artifacts_sha256": evaluator_record["artifacts_sha256"]}, "surface": surface, "evaluation": config["evaluation"]})
    evaluation = output / "evaluator"
    _invoke(root, output, "evaluate", "evaluate", evaluate_binding, evaluation)
    evaluator_receipt = _json(evaluation / "receipt.json")
    if evaluator_receipt.get("unit_metrics_written") is not True or evaluator_receipt.get("partial_ranking_authorized") is not False:
        raise SeedBundleError("surface evaluator disposition differs")
    receipt = {"schema_version": "masld-bench-observed-multiome-factorized-production-surface-receipt-v1", "surface": surface, "models": list(MODEL_IDS), "prediction_committed_before_evaluator": True, "unit_metrics_written": True, "partial_ranking_authorized": False, "promotion_authorized": False, "sealed_outcomes_read": False}
    write_json_exclusive(output / "receipt.json", receipt)
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_factorized_production_surface", "held_genomic_fold": genomic, "held_donor_fold": donor, "seed": seed, "partial_ranking_authorized": False})
    return {"held_genomic_fold": genomic, "held_donor_fold": donor, "seed": seed, "path": output.name, "artifacts_sha256": digest}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--contract-artifact", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    config_path.relative_to(root)
    contract_path = args.contract_artifact.resolve(strict=True)
    contract_path.relative_to(root)
    config = _json(config_path)
    inputs = validate(root, config, contract_path, args.seed)
    output = reject_symlink_components(args.output, label="factorized seed bundle output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    surfaces = []
    for planned in planned_surfaces(args.seed):
        name = f"surface_g{planned['held_genomic_fold']}_d{planned['held_donor_fold']}_s{args.seed}"
        surfaces.append(_run_surface(root, config, inputs, args.seed, planned["held_genomic_fold"], planned["held_donor_fold"], output / name))
    receipt = {"schema_version": "masld-bench-observed-multiome-factorized-seed-bundle-receipt-v1", "contract_id": CONTRACT_ID, "dataset_id": "gse296875", "seed": args.seed, "surfaces": surfaces, "surface_count": len(surfaces), "models_per_surface": len(MODEL_IDS), "prediction_matrices": len(surfaces) * len(MODEL_IDS), "all_surfaces_complete": True, "partial_ranking_authorized": False, "promotion_authorized": False, "sealed_outcomes_read": False, "next_gate": "all_seed_bundles_then_independent_aggregation"}
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [config_path, Path(__file__).resolve(strict=True), root / "scripts/run_gse296875_observed_multiome_factorized_surface_smoke.py", root / "tests/unit/test_observed_multiome_factorized_seed_bundle.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_factorized_seed_bundle", "contract_id": CONTRACT_ID, "seed": args.seed, "surface_count": len(surfaces), "partial_ranking_authorized": False, "sealed_outcomes_accessed": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()

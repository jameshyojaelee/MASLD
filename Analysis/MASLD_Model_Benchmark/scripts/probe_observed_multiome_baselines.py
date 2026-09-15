#!/usr/bin/env python3
"""Run a biological-free masked-profile fixture for five modality controls."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys
from typing import Any

import numpy as np

from masld_bench.adapters.observed_multiome_baselines import MODEL_IDS, export, fit, predict
from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive


SCHEMA = "masld-bench-observed-multiome-baseline-fixture-v1"


class BaselineFixtureError(ValueError):
    """Raised when the synthetic fixture or its census binding differs."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise BaselineFixtureError("fixture config must be a JSON object")
    return value


def _synthetic(config: dict[str, Any]) -> dict[str, Any]:
    dimensions = config["synthetic_dimensions"]
    generator = np.random.default_rng(int(config["seed"]))
    n_train = int(dimensions["training_rows"])
    n_query = int(dimensions["query_rows"])
    n_rna = int(dimensions["rna_features"])
    n_sequence = int(dimensions["sequence_features"])
    n_atac = int(dimensions["unmasked_atac_features"])
    n_target = int(dimensions["masked_target_bins"])
    if min(n_train, n_query, n_rna, n_sequence, n_atac, n_target) < 2:
        raise BaselineFixtureError("synthetic dimensions are too small")
    train_rna = generator.normal(size=(n_train, n_rna))
    query_rna = generator.normal(size=(n_query, n_rna))
    train_sequence = generator.normal(size=(n_train, n_sequence))
    query_sequence = generator.normal(size=(n_query, n_sequence))
    train_atac = generator.poisson(3.0, size=(n_train, n_atac)).astype(float)
    query_atac = generator.poisson(3.0, size=(n_query, n_atac)).astype(float)
    weights = generator.normal(scale=0.12, size=(n_atac, n_target))
    rates = np.exp(np.clip(0.2 + train_atac @ weights, -2.0, 3.0))
    train_target = generator.poisson(rates).astype(float)
    return {
        "train_row_ids": [f"train-{index:03d}" for index in range(n_train)],
        "query_row_ids": [f"query-{index:03d}" for index in range(n_query)],
        "train_strata": [f"lineage-{index % 2}" for index in range(n_train)],
        "query_strata": [f"lineage-{index % 2}" for index in range(n_query)],
        "train_rna": train_rna,
        "query_rna": query_rna,
        "train_sequence": train_sequence,
        "query_sequence": query_sequence,
        "train_atac": train_atac,
        "query_atac": query_atac,
        "train_target": train_target,
    }


def run(root: Path, config: dict[str, Any], output: Path) -> dict[str, Any]:
    if config.get("schema_version") != SCHEMA:
        raise BaselineFixtureError("fixture schema differs")
    if tuple(config.get("model_ids", ())) != MODEL_IDS:
        raise BaselineFixtureError("fixture model roster differs")
    firewall = config.get("outcome_firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise BaselineFixtureError("fixture outcome firewall is open")
    census = config.get("census_artifact")
    if not isinstance(census, dict):
        raise BaselineFixtureError("census binding is absent")
    census_path = (root / str(census["path"])).resolve(strict=True)
    census_path.relative_to(root)
    manifest = verify_frozen_tree(census_path)
    if _digest(census_path / "ARTIFACTS.json") != census.get("artifacts_sha256"):
        raise BaselineFixtureError("census artifact binding drifted")
    receipt = _load(census_path / "receipt.json")
    if receipt.get("census_requirement_satisfied") is not True or receipt.get("execution_authorized") is not False:
        raise BaselineFixtureError("census disposition differs")

    values = _synthetic(config)
    fit_config = config["fit"]
    receipts = []
    for model_id in MODEL_IDS:
        state = fit(
            model_id,
            row_ids=values["train_row_ids"],
            strata=values["train_strata"],
            rna=values["train_rna"],
            sequence=values["train_sequence"],
            observed_atac=values["train_atac"],
            target_counts=values["train_target"],
            seed=int(config["seed"]),
            ridge_alpha=float(fit_config["ridge_alpha"]),
            poisson_alpha=float(fit_config["poisson_alpha"]),
            poisson_max_iter=int(fit_config["poisson_max_iter"]),
            poisson_tolerance=float(fit_config["poisson_tolerance"]),
        )
        prediction = predict(
            state,
            row_ids=values["query_row_ids"],
            strata=values["query_strata"],
            rna=values["query_rna"],
            sequence=values["query_sequence"],
            observed_atac=values["query_atac"],
        )
        model_path = output / "states" / model_id
        export(state, model_path)
        prediction_path = output / "predictions" / f"{model_id}.npy"
        prediction_path.parent.mkdir(parents=True, exist_ok=True)
        with prediction_path.open("xb") as handle:
            np.save(handle, prediction, allow_pickle=False)
        receipts.append(
            {
                "model_id": model_id,
                "prediction_shape": list(prediction.shape),
                "prediction_sha256": _digest(prediction_path),
                "state_npz_sha256": _digest(model_path / "state.npz"),
                "finite_nonnegative": bool(np.isfinite(prediction).all() and np.all(prediction >= 0)),
            }
        )
    if any(not record["finite_nonnegative"] for record in receipts):
        raise BaselineFixtureError("a fixture prediction is invalid")
    result = {
        "schema_version": "masld-bench-observed-multiome-baseline-fixture-receipt-v1",
        "fixture_id": config["fixture_id"],
        "census_artifacts_sha256": census["artifacts_sha256"],
        "model_count": len(receipts),
        "models": receipts,
        "biological_data_read": False,
        "development_outcomes_read": False,
        "sealed_outcomes_read": False,
        "benchmark_metrics_calculated": False,
        "candidate_execution_authorized": False,
        "baseline_implementation_fixture_passed": True,
        "runtime": {
            "python_executable": str(Path(sys.executable).resolve(strict=True)),
            "python_version": platform.python_version(),
            "numpy_version": np.__version__,
        },
    }
    write_json_exclusive(output / "receipt.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    config_path.relative_to(root)
    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    config = _load(config_path)
    result = run(root, config, output)
    source_paths = [
        config_path,
        Path(__file__).resolve(strict=True),
        root / "src/masld_bench/adapters/observed_multiome_baselines.py",
        root / "tests/unit/test_observed_multiome_baselines.py",
    ]
    (output / "source.sha256").write_text(
        "".join(f"{_digest(path)}  {path}\n" for path in source_paths),
        encoding="utf-8",
    )
    digest = freeze_tree(
        output,
        metadata={
            "artifact_class": "observed_multiome_baseline_fixture",
            "fixture_id": result["fixture_id"],
            "model_count": result["model_count"],
            "candidate_execution_authorized": False,
            "sealed_outcomes_accessed": False,
        },
    )
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()

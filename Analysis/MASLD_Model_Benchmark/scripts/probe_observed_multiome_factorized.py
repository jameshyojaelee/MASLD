#!/usr/bin/env python3
"""Run a biological-free held-target fixture for factorized baselines."""

from __future__ import annotations

import argparse
from hashlib import sha256
import inspect
import json
from pathlib import Path
import platform
import sys
from typing import Any

import numpy as np

from masld_bench.adapters.observed_multiome_factorized import MODEL_IDS, export, fit, predict
from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive


SCHEMA = "masld-bench-observed-multiome-factorized-fixture-v1"


class FactorizedFixtureError(ValueError):
    """Raised when the target-generalizing fixture or parent binding differs."""


def _load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise FactorizedFixtureError("factorized fixture config must be an object")
    return value


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _synthetic(config: dict[str, Any]) -> dict[str, Any]:
    d = config["synthetic_dimensions"]
    generator = np.random.default_rng(int(config["seed"]))
    n_train = int(d["training_rows"])
    n_query = int(d["query_rows"])
    n_train_target = int(d["training_targets"])
    n_held_target = int(d["held_targets"])
    rna = generator.normal(size=(n_train + n_query, int(d["rna_features"])))
    atac = generator.poisson(3.0, size=(n_train + n_query, int(d["unmasked_atac_features"]))).astype(float)
    sequence = generator.normal(size=(n_train_target + n_held_target, int(d["sequence_features"])))
    covariates = generator.normal(size=(n_train_target + n_held_target, int(d["target_covariates"])))
    row_latent = rna[:n_train, :2]
    target_latent = sequence[:n_train_target, :2]
    rate = np.exp(np.clip(0.2 + row_latent @ target_latent.T * 0.15, -2.0, 2.0))
    counts = generator.poisson(rate).astype(float)
    return {
        "train_rows": [f"train-row-{index:03d}" for index in range(n_train)],
        "query_rows": [f"query-row-{index:03d}" for index in range(n_query)],
        "train_strata": [f"lineage-{index % 2}" for index in range(n_train)],
        "query_strata": [f"lineage-{index % 2}" for index in range(n_query)],
        "train_targets": [f"train-target-{index:03d}" for index in range(n_train_target)],
        "held_targets": [f"held-target-{index:03d}" for index in range(n_held_target)],
        "train_rna": rna[:n_train],
        "query_rna": rna[n_train:],
        "train_atac": atac[:n_train],
        "query_atac": atac[n_train:],
        "train_sequence": sequence[:n_train_target],
        "held_sequence": sequence[n_train_target:],
        "train_covariates": covariates[:n_train_target],
        "held_covariates": covariates[n_train_target:],
        "counts": counts,
    }


def run(root: Path, config: dict[str, Any], output: Path) -> dict[str, Any]:
    if config.get("schema_version") != SCHEMA or tuple(config.get("model_ids", ())) != MODEL_IDS:
        raise FactorizedFixtureError("factorized fixture identity differs")
    firewall = config.get("outcome_firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise FactorizedFixtureError("factorized fixture firewall is open")
    parent = config.get("target_axis_artifact")
    if not isinstance(parent, dict):
        raise FactorizedFixtureError("target-axis parent is absent")
    parent_path = (root / str(parent.get("path", ""))).resolve(strict=True)
    parent_path.relative_to(root)
    verify_frozen_tree(parent_path)
    if _digest(parent_path / "ARTIFACTS.json") != parent.get("artifacts_sha256"):
        raise FactorizedFixtureError("target-axis parent drifted")
    parent_receipt = _load(parent_path / "receipt.json")
    if parent_receipt.get("production_target_axis_frozen") is not True:
        raise FactorizedFixtureError("target-axis parent is not frozen")
    if "target_counts" in inspect.signature(predict).parameters:
        raise FactorizedFixtureError("prediction accepts outcomes")

    values = _synthetic(config)
    fit_config = config["fit"]
    receipts = []
    for model_id in MODEL_IDS:
        state = fit(
            model_id,
            row_ids=values["train_rows"],
            strata=values["train_strata"],
            target_ids=values["train_targets"],
            rna=values["train_rna"],
            observed_atac=values["train_atac"],
            sequence_features=values["train_sequence"],
            target_covariates=values["train_covariates"],
            target_counts=values["counts"],
            seed=int(config["seed"]),
            ridge_alpha=float(fit_config["ridge_alpha"]),
            poisson_alpha=float(fit_config["poisson_alpha"]),
            poisson_max_iter=int(fit_config["poisson_max_iter"]),
            poisson_tolerance=float(fit_config["poisson_tolerance"]),
            maximum_expanded_pairs=int(fit_config["maximum_expanded_pairs"]),
        )
        prediction = predict(
            state,
            row_ids=values["query_rows"],
            strata=values["query_strata"],
            target_ids=values["held_targets"],
            rna=values["query_rna"],
            observed_atac=values["query_atac"],
            sequence_features=values["held_sequence"],
            target_covariates=values["held_covariates"],
            maximum_expanded_pairs=int(fit_config["maximum_expanded_pairs"]),
        )
        state_path = output / "states" / model_id
        export(state, state_path)
        prediction_path = output / "predictions" / f"{model_id}.npy"
        prediction_path.parent.mkdir(parents=True, exist_ok=True)
        with prediction_path.open("xb") as handle:
            np.save(handle, prediction, allow_pickle=False)
        receipts.append(
            {
                "model_id": model_id,
                "prediction_shape": list(prediction.shape),
                "prediction_sha256": _digest(prediction_path),
                "coefficient_shape": list(state.coefficients.shape),
                "target_specific_coefficients": False,
                "training_target_count": len(state.training_target_ids),
                "held_target_count": len(values["held_targets"]),
                "training_held_target_overlap": 0,
            }
        )
    receipt = {
        "schema_version": "masld-bench-observed-multiome-factorized-fixture-receipt-v1",
        "fixture_id": config["fixture_id"],
        "model_count": len(receipts),
        "models": receipts,
        "held_target_prediction_passed": True,
        "target_specific_coefficients": False,
        "biological_data_read": False,
        "development_outcome_read": False,
        "sealed_outcome_read": False,
        "benchmark_metric_calculated": False,
        "biological_execution_authorized": False,
        "runtime": {
            "python_executable": str(Path(sys.executable).resolve(strict=True)),
            "python_version": platform.python_version(),
            "numpy_version": np.__version__,
        },
    }
    write_json_exclusive(output / "receipt.json", receipt)
    return receipt


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
    receipt = run(root, _load(config_path), output)
    sources = [
        config_path,
        Path(__file__).resolve(strict=True),
        root / "src/masld_bench/adapters/observed_multiome_factorized.py",
        root / "tests/unit/test_observed_multiome_factorized.py",
    ]
    (output / "source.sha256").write_text(
        "".join(f"{_digest(path)}  {path}\n" for path in sources),
        encoding="utf-8",
    )
    digest = freeze_tree(
        output,
        metadata={
            "artifact_class": "observed_multiome_factorized_fixture",
            "fixture_id": receipt["fixture_id"],
            "model_count": receipt["model_count"],
            "biological_execution_authorized": False,
            "sealed_outcomes_accessed": False,
        },
    )
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()

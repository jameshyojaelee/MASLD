#!/usr/bin/env python3
"""Freeze an outcome-free roster of completed scBasset fold predictions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile
from typing import Any

from masld_bench.artifacts import (
    freeze_tree,
    publish_directory_noreplace,
    verify_frozen_tree,
    write_json_exclusive,
)
from masld_bench.hashing import sha256_file


EXPECTED_FOLDS = (1, 2, 3, 4)
EXPECTED_SEED = 20260824


class ScBassetPredictionLockError(ValueError):
    """Raised when a chain/input/prediction binding differs."""


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScBassetPredictionLockError(f"invalid {label}: {path}") from error
    if not isinstance(value, dict):
        raise ScBassetPredictionLockError(f"{label} must be an object")
    return value


def _verify_root(path: Path, label: str) -> tuple[Path, str]:
    root = path.resolve(strict=True)
    verify_frozen_tree(root)
    return root, sha256_file(root / "ARTIFACTS.json")


def build_lock(
    *, chain_root: Path, input_roots: list[Path], output: Path
) -> dict[str, Any]:
    if output.exists():
        raise ScBassetPredictionLockError("prediction-lock output already exists")
    chain_root, chain_sha256 = _verify_root(chain_root, "chain receipt")
    chain = _read_json(chain_root / "bundle_receipt.json", "chain receipt")
    if (
        chain.get("schema_version") != "masld-bench-scbasset-chain-receipt-v1"
        or chain.get("bundle_id") != "scbasset-fold1to4-chain-v1"
        or chain.get("status") != "succeeded"
        or chain.get("held_donor_atac_used") is not False
        or chain.get("genomic_test_atac_used") is not False
    ):
        raise ScBassetPredictionLockError("chain receipt contract differs")
    dispositions = chain.get("dispositions")
    if not isinstance(dispositions, list) or [
        int(item.get("fold", -1)) for item in dispositions
    ] != list(EXPECTED_FOLDS):
        raise ScBassetPredictionLockError("chain fold roster differs")

    input_by_fold: dict[int, tuple[Path, str, dict[str, Any]]] = {}
    for path in input_roots:
        root, artifacts_sha256 = _verify_root(path, "fold input")
        summary = _read_json(root / "inputs/summary.json", "input summary")
        fold = int(summary.get("donor_test_fold", -1))
        if (
            fold not in EXPECTED_FOLDS
            or fold in input_by_fold
            or summary.get("status") != "pass"
            or summary.get("split_id") != f"donor{fold}_genomic{fold}"
            or summary.get("donor_valid_fold") != (fold + 1) % 5
            or summary.get("held_donor_atac_exported") is not False
            or summary.get("genomic_test_atac_exported") is not False
        ):
            raise ScBassetPredictionLockError("fold input contract differs")
        input_by_fold[fold] = (root, artifacts_sha256, summary)
    if set(input_by_fold) != set(EXPECTED_FOLDS):
        raise ScBassetPredictionLockError("input fold roster differs")

    records = []
    for fold, disposition in zip(EXPECTED_FOLDS, dispositions):
        if (
            disposition.get("train_returncode") != 0
            or disposition.get("predict_returncode") != 0
        ):
            raise ScBassetPredictionLockError(f"fold {fold} is not successful")
        model_root, model_sha256 = _verify_root(
            Path(str(disposition.get("model_path", ""))), f"fold {fold} model"
        )
        prediction_root, prediction_sha256 = _verify_root(
            Path(str(disposition.get("prediction_path", ""))),
            f"fold {fold} prediction",
        )
        if (
            model_sha256 != disposition.get("model_artifacts_sha256")
            or prediction_sha256 != disposition.get("prediction_artifacts_sha256")
        ):
            raise ScBassetPredictionLockError(f"fold {fold} receipt hash differs")
        input_root, input_sha256, summary = input_by_fold[fold]
        model = _read_json(model_root / "model/model_receipt.json", "model receipt")
        prediction = _read_json(
            prediction_root / "predictions/prediction_receipt.json",
            "prediction receipt",
        )
        if (
            model.get("status") != "pass"
            or model.get("model_id") != "scbasset"
            or model.get("split_id") != f"donor{fold}_genomic{fold}"
            or model.get("seed") != EXPECTED_SEED
            or model.get("input_artifacts_sha256") != input_sha256
            or model.get("held_donor_atac_used") is not False
            or model.get("genomic_test_atac_used") is not False
            or prediction.get("status") != "pass"
            or prediction.get("role") != "valid"
            or prediction.get("split_id") != f"donor{fold}_genomic{fold}"
            or prediction.get("input_artifacts_sha256") != input_sha256
            or prediction.get("model_artifacts_sha256") != model_sha256
            or prediction.get("held_donor_atac_used") is not False
            or prediction.get("held_cell_embedding_available") is not False
            or prediction.get("training_cells") != summary.get("training_cells")
        ):
            raise ScBassetPredictionLockError(f"fold {fold} model/prediction differs")
        records.append(
            {
                "outer_fold": fold,
                "valid_fold": (fold + 1) % 5,
                "split_id": f"donor{fold}_genomic{fold}",
                "seed": EXPECTED_SEED,
                "input_root": input_root.as_posix(),
                "input_artifacts_sha256": input_sha256,
                "model_root": model_root.as_posix(),
                "model_artifacts_sha256": model_sha256,
                "prediction_root": prediction_root.as_posix(),
                "prediction_artifacts_sha256": prediction_sha256,
                "held_donors": int(prediction["held_donors"]),
                "training_cells": int(prediction["training_cells"]),
                "regions": int(prediction["regions"]),
            }
        )

    output_parent = output.parent.resolve(strict=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output_parent))
    try:
        payload = {
            "schema_version": "masld-bench-scbasset-crossfold-prediction-lock-v1",
            "status": "pass",
            "dataset_id": "gse296875",
            "task_id": "rna_conditioned_atac",
            "model_id": "scbasset",
            "seed": EXPECTED_SEED,
            "folds": list(EXPECTED_FOLDS),
            "chain_root": chain_root.as_posix(),
            "chain_artifacts_sha256": chain_sha256,
            "records": records,
            "outcomes_read": False,
            "outcome_paths_recorded": False,
            "champion_claim_allowed": False,
        }
        write_json_exclusive(stage / "prediction_lock.json", payload)
        digest = freeze_tree(
            stage,
            {
                "artifact_class": "scbasset_crossfold_prediction_lock",
                "dataset_id": "gse296875",
                "model_id": "scbasset",
                "seed": EXPECTED_SEED,
                "folds": list(EXPECTED_FOLDS),
                "chain_artifacts_sha256": chain_sha256,
                "outcomes_read": False,
                "status": "pass",
            },
        )
        publish_directory_noreplace(stage, output)
    except Exception:
        raise
    return {
        "output": output.as_posix(),
        "artifacts_sha256": digest,
        "records": len(records),
        "status": "pass",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chain-root", type=Path, required=True)
    parser.add_argument(
        "--input-root",
        dest="input_roots",
        type=Path,
        action="append",
        required=True,
    )
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    print(json.dumps(build_lock(**vars(arguments)), sort_keys=True))


if __name__ == "__main__":
    main()

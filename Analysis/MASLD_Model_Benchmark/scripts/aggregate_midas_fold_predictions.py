#!/usr/bin/env python3
"""Aggregate five frozen single-fold MIDAS predictions without reading outcomes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
from typing import Any, Mapping, Sequence


class MIDASAggregationError(ValueError):
    """Raised when frozen MIDAS fold artifacts cannot be combined exactly."""


def parse_fold_roots(values: Sequence[str]) -> dict[int, Path]:
    roots: dict[int, Path] = {}
    for value in values:
        fold_text, separator, root_text = value.partition("=")
        if not separator or not fold_text.isdigit():
            raise MIDASAggregationError("fold-root syntax differs")
        fold = int(fold_text)
        if fold not in range(5) or fold in roots:
            raise MIDASAggregationError("fold-root census differs")
        roots[fold] = Path(root_text)
    if set(roots) != set(range(5)):
        raise MIDASAggregationError("exactly five folds are required")
    return roots


def aggregate(
    fold_roots: Mapping[int, Path], output: Path, expected_seed: int
) -> dict[str, Any]:
    if output.exists() or expected_seed < 0:
        raise MIDASAggregationError("output or seed contract differs")
    output.mkdir(parents=True, mode=0o750)
    combined: dict[str, Any] = {}
    for fold in range(5):
        root = fold_roots[fold]
        receipt_path = root / "predictions/receipt.json"
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if (
            receipt.get("status") != "pass"
            or receipt.get("model_id") != "midas_inductive"
            or receipt.get("exact_distribution") != "scmidas 0.3.0"
            or receipt.get("pairing_topology") != "same_nucleus_training"
            or receipt.get("outer_unit") != "donor"
            or set(receipt.get("folds", {})) != {str(fold)}
            or receipt.get("base_seed") != expected_seed
            or receipt.get("query_rna_used_for_training") is not False
            or receipt.get("held_atac_input_exposed") is not False
            or receipt.get("outcomes_read") is not False
            or receipt.get("test_outcomes_read") is not False
        ):
            raise MIDASAggregationError(f"fold {fold} receipt differs")
        record = receipt["folds"][str(fold)]
        source = root / "predictions" / f"fold_{fold}"
        if source.is_symlink() or not source.is_dir():
            raise MIDASAggregationError(f"fold {fold} source differs")
        bundle = json.loads((source / "prediction_bundle.json").read_text())
        if (
            bundle.get("model_id") != "midas_inductive"
            or bundle.get("metadata", {}).get("held_out_fold") != fold
            or bundle.get("metadata", {}).get("seed") != expected_seed
            or bundle.get("metadata", {}).get("query_rna_used_for_training") is not False
            or bundle.get("metadata", {}).get("held_atac_input_exposed") is not False
        ):
            raise MIDASAggregationError(f"fold {fold} bundle differs")
        shutil.copytree(source, output / f"fold_{fold}", symlinks=False)
        combined[str(fold)] = record
    receipt = {
        "schema_version": "masld-bench-midas-inductive-aggregate-v1",
        "status": "pass",
        "model_id": "midas_inductive",
        "task_id": "rna_conditioned_atac",
        "dataset_id": "gse296875",
        "pairing_topology": "same_nucleus_training",
        "outer_unit": "donor",
        "folds": combined,
        "base_seed": expected_seed,
        "fold_count": 5,
        "prediction_frozen_before_outcomes": True,
        "query_rna_used_for_training": False,
        "held_atac_input_exposed": False,
        "outcomes_read": False,
        "test_outcomes_read": False,
        "sealed_inference_eligible": False,
        "champion_claim_allowed": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fold-root", action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-seed", type=int, required=True)
    arguments = parser.parse_args()
    aggregate(
        parse_fold_roots(arguments.fold_root), arguments.output, arguments.expected_seed
    )


if __name__ == "__main__":
    main()

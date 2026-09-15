#!/usr/bin/env python3
"""Validate and combine prediction-only sequence-control outer-fold bundles."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


FOLDS = tuple(f"fold-{index}" for index in range(5))
MODELS = ("sequence_cnn_control", "sequence_transformer_control")
SEEDS = (1103, 2909, 4721, 6673, 8111)
CONTEXTS = ("HepG2_control", "HepG2_PAOA")


class SequenceAggregationError(RuntimeError):
    """Raised when prediction-only output files differ from the frozen requirements."""


def file_sha256(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SequenceAggregationError("JSON artifact must be an object")
    return value


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return tuple(reader.fieldnames or ()), [dict(row) for row in reader]


def write_gzip_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
    buffer = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buffer, mtime=0) as zipped:
        with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
            writer = csv.DictWriter(
                text, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
    path.write_bytes(buffer.getvalue())


def aggregate(arguments: argparse.Namespace) -> dict[str, Any]:
    if arguments.output.exists():
        raise SequenceAggregationError("aggregate output exists")
    parts = arguments.parts.resolve(strict=True)
    if any(path.name == "training_targets.tsv" for path in parts.rglob("*")):
        raise SequenceAggregationError("raw training target leaked into frozen parts")
    rows: list[dict[str, str]] = []
    fields: tuple[str, ...] | None = None
    checkpoint_hashes: dict[tuple[str, str], set[str]] = {}
    checkpoint_count = 0
    training_receipt_count = 0
    part_receipts: list[dict[str, Any]] = []
    for fold in FOLDS:
        part = parts / fold
        receipt = read_json(part / "receipt.json")
        if (
            receipt.get("status") != "pass_prediction_only_outer_fold_bundle"
            or receipt.get("outer_test_fold") != fold
            or receipt.get("fit_count") != 10
            or receipt.get("checkpoint_count") != 10
            or receipt.get("training_receipt_count") != 10
            or receipt.get("held_test_outcomes_visible_to_adapter") is not False
            or receipt.get("raw_training_targets_retained") is not False
            or receipt.get("benchmark_metrics_calculated")
            or receipt.get("evaluator_invoked")
            or receipt.get("model_ranking_performed")
        ):
            raise SequenceAggregationError("outer-fold receipt differs")
        part_receipts.append(receipt)
        part_fields, part_rows = read_tsv(part / "predictions.tsv")
        if fields is None:
            fields = part_fields
        if part_fields != fields or len(part_rows) != receipt["prediction_rows"]:
            raise SequenceAggregationError("outer-fold prediction schema differs")
        rows.extend(part_rows)
        for model in MODELS:
            checkpoints = sorted(
                (part / "checkpoints" / model).glob("seed-*.safetensors")
            )
            receipts = sorted(
                (part / "training_receipts" / model).glob("seed-*.json")
            )
            if len(checkpoints) != 5 or len(receipts) != 5:
                raise SequenceAggregationError("fit artifact census differs")
            checkpoint_count += len(checkpoints)
            training_receipt_count += len(receipts)
            observed_hashes: set[str] = set()
            for checkpoint, fit_receipt_path in zip(
                checkpoints, receipts, strict=True
            ):
                fit_receipt = read_json(fit_receipt_path)
                observed_hash = file_sha256(checkpoint)
                if (
                    fit_receipt.get("status") != "pass_prediction_only_fit"
                    or fit_receipt.get("outer_test_fold") != fold
                    or fit_receipt.get("model_id") != model
                    or fit_receipt.get("checkpoint_sha256") != observed_hash
                    or fit_receipt.get("held_test_outcomes_visible_to_adapter")
                    is not False
                    or fit_receipt.get("benchmark_metrics_calculated")
                    or fit_receipt.get("evaluator_invoked")
                    or fit_receipt.get("model_ranked")
                ):
                    raise SequenceAggregationError("fit receipt differs")
                observed_hashes.add(observed_hash)
            if len(observed_hashes) != 5:
                raise SequenceAggregationError("seed checkpoints are not distinct")
            checkpoint_hashes[(model, fold)] = observed_hashes
    if fields is None:
        raise SequenceAggregationError("no prediction parts")
    required_fields = {
        "element_id",
        "source_locus_group_id",
        "long_range_block_id",
        "outer_fold",
        "model_id",
        "implementation_id",
        "seed",
        "context_id",
        "predicted_reference_activity",
        "predicted_alternative_activity",
        "predicted_alt_minus_ref_activity",
        "prediction_units",
        "checkpoint_sha256",
        "benchmark_metrics_calculated",
    }
    keys = [
        (row["model_id"], row["seed"], row["element_id"], row["context_id"])
        for row in rows
    ]
    if (
        len(rows) != 20660
        or len(set(keys)) != 20660
        or not required_fields.issubset(fields)
        or set(row["model_id"] for row in rows) != set(MODELS)
        or set(int(row["seed"]) for row in rows) != set(SEEDS)
        or set(row["context_id"] for row in rows) != set(CONTEXTS)
        or set(row["outer_fold"] for row in rows) != set(FOLDS)
        or any(row["benchmark_metrics_calculated"] != "false" for row in rows)
        or checkpoint_count != 50
        or training_receipt_count != 50
        or len(checkpoint_hashes) != 10
    ):
        raise SequenceAggregationError("campaign prediction contract differs")
    elements = {row["element_id"] for row in rows}
    if len(elements) != 1033:
        raise SequenceAggregationError("campaign element census differs")
    rows.sort(
        key=lambda row: (
            row["model_id"],
            int(row["seed"]),
            row["element_id"],
            row["context_id"],
        )
    )
    arguments.output.mkdir(mode=0o750)
    prediction_dir = arguments.output / "predictions"
    prediction_dir.mkdir(mode=0o750)
    prediction_path = prediction_dir / "predictions.tsv.gz"
    write_gzip_tsv(prediction_path, fields, rows)
    receipt = {
        "schema_version": "masld-bench-sequence-control-production-receipt-v1",
        "status": "pass_prediction_only_production_campaign",
        "dataset_id": "gse281364",
        "models": list(MODELS),
        "outer_folds": list(FOLDS),
        "genuine_training_seeds": list(SEEDS),
        "fit_count": 50,
        "inner_selection_fit_count": 50,
        "final_refit_count": 50,
        "checkpoint_count": checkpoint_count,
        "training_receipt_count": training_receipt_count,
        "prediction_rows": len(rows),
        "elements": len(elements),
        "contexts": list(CONTEXTS),
        "predictions_sha256": file_sha256(prediction_path),
        "checkpoint_hashes_distinct_within_model_outer_fold": True,
        "held_test_outcomes_visible_to_training_adapter": False,
        "raw_training_targets_retained": False,
        "outcomes_copied_to_output": False,
        "benchmark_metrics_calculated": False,
        "evaluator_invoked": False,
        "model_ranking_performed": False,
        "shortlist_performed": False,
        "promotion_performed": False,
        "champion_eligible": False,
        "external_claim_eligible": False,
        "part_elapsed_seconds": {
            receipt["outer_test_fold"]: receipt["elapsed_seconds"]
            for receipt in part_receipts
        },
    }
    (arguments.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parts", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    aggregate(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Combine three read-only scVI/scANVI seed predictions without scoring."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from masld_bench.adapters.scvi_scanvi_inductive import ROSTER, validate_probabilities
from scripts.fit_predict_scvi_scanvi_study_50000 import (
    DATASET_VIEW_ID,
    MODEL_IDS,
    SEEDS,
    SPLIT_ID,
    write_json_exclusive,
)


class AggregateSCVIError(ValueError):
    """Raised when read-only seed outputs cannot be combined exactly."""


FIELDS = (
    "row_id",
    "donor_id",
    "dataset",
    "outer_fold",
    "predicted_class",
    *(f"probability::{label}" for label in ROSTER),
)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise AggregateSCVIError(f"prediction schema differs: {path}")
        return list(reader)


def _write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def run(*, seed_roots: Sequence[Path], output: Path) -> None:
    if output.exists():
        raise AggregateSCVIError("refusing to overwrite seed ensemble")
    if len(seed_roots) != len(SEEDS):
        raise AggregateSCVIError("three seed roots are required")
    roots_by_seed: dict[int, Path] = {}
    for root in seed_roots:
        receipt = json.loads(
            (root / "prediction_receipt.json").read_text(encoding="utf-8")
        )
        seed = int(receipt.get("seed", -1))
        if (
            seed not in SEEDS
            or seed in roots_by_seed
            or receipt.get("status") != "pass_development_predictions"
            or receipt.get("dataset_view_id") != DATASET_VIEW_ID
            or receipt.get("split_id") != SPLIT_ID
            or receipt.get("models") != list(MODEL_IDS)
            or receipt.get("metrics_calculated") is not False
            or receipt.get("common_head_fit") is not False
            or receipt.get("held_query_training") is not False
            or receipt.get("held_query_adaptation") is not False
            or receipt.get("sealed_outcomes_read") is not False
        ):
            raise AggregateSCVIError(f"seed receipt differs: {root}")
        roots_by_seed[seed] = root
    if set(roots_by_seed) != set(SEEDS):
        raise AggregateSCVIError("seed roster differs")

    output.mkdir(mode=0o750)
    prediction_root = output / "predictions"
    prediction_root.mkdir()
    roster_fields = [f"probability::{label}" for label in ROSTER]
    row_inventory: list[tuple[str, str, str, str]] | None = None
    for model_id in MODEL_IDS:
        seed_values: list[np.ndarray] = []
        seed_rows: list[list[dict[str, str]]] = []
        for seed in SEEDS:
            rows = _read_tsv(
                roots_by_seed[seed]
                / "predictions"
                / f"{model_id}--seed-{seed}.tsv"
            )
            if len(rows) != 50_000:
                raise AggregateSCVIError("seed prediction row count differs")
            inventory = [
                (row["row_id"], row["donor_id"], row["dataset"], row["outer_fold"])
                for row in rows
            ]
            if row_inventory is None:
                row_inventory = inventory
            elif inventory != row_inventory:
                raise AggregateSCVIError("seed prediction row inventory differs")
            values = validate_probabilities(
                np.asarray(
                    [[float(row[field]) for field in roster_fields] for row in rows],
                    dtype=np.float64,
                ),
                rows=len(rows),
            )
            seed_rows.append(rows)
            seed_values.append(values)
            _write_tsv(
                prediction_root / f"{model_id}--seed-{seed}.tsv", rows
            )
        ensemble = validate_probabilities(
            np.mean(np.stack(seed_values, axis=0), axis=0), rows=50_000
        )
        template = seed_rows[0]
        _write_tsv(
            prediction_root / f"{model_id}.tsv",
            [
                {
                    "row_id": template[index]["row_id"],
                    "donor_id": template[index]["donor_id"],
                    "dataset": template[index]["dataset"],
                    "outer_fold": template[index]["outer_fold"],
                    "predicted_class": ROSTER[int(np.argmax(ensemble[index]))],
                    **{
                        field: format(float(ensemble[index, class_index]), ".17g")
                        for class_index, field in enumerate(roster_fields)
                    },
                }
                for index in range(len(template))
            ],
        )

    write_json_exclusive(
        output / "prediction_receipt.json",
        {
            "schema_version": "masld-bench-cell-baseline-study-screen-v1",
            "status": "pass_development_predictions",
            "dataset_view_id": DATASET_VIEW_ID,
            "split_id": SPLIT_ID,
            "rows": 50_000,
            "donors": 102,
            "studies": 7,
            "outer_folds": 5,
            "inner_folds": 5,
            "models": list(MODEL_IDS),
            "parameters": {
                "seeds": list(SEEDS),
                "final_prediction": "unweighted_mean_probability_ensemble",
            },
            "metrics_calculated": False,
            "common_head_fit": False,
            "held_query_training": False,
            "held_query_adaptation": False,
            "prediction_tables_contain_observed_labels": False,
            "histology_read": False,
            "sealed_outcomes_read": False,
        },
    )
    write_json_exclusive(
        output / "ensemble_receipt.json",
        {
            "schema_version": "masld-bench-scvi-scanvi-seed-ensemble-v1",
            "status": "pass_development_predictions",
            "dataset_view_id": DATASET_VIEW_ID,
            "split_id": SPLIT_ID,
            "rows": 50_000,
            "models": list(MODEL_IDS),
            "seeds": list(SEEDS),
            "ensemble": "unweighted_mean_probability",
            "metrics_calculated": False,
            "common_head_fit": False,
            "held_query_training": False,
            "held_query_adaptation": False,
            "sealed_outcomes_read": False,
        },
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-root", required=True, action="append", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    run(seed_roots=arguments.seed_root, output=arguments.output)
    print(json.dumps({"output": str(arguments.output.resolve()), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

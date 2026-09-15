#!/usr/bin/env python3
"""Materialise the evaluator-only GSE49541 fibrosis outcome after predictions freeze.

The deposited fibrosis group is an outcome and belongs to the evaluator alone.
This step exists so the ordering is a fact on disk rather than a claim: it will
not emit a label table until it has verified the committed SHA-256 of every
already-frozen prediction bundle, so no label can have informed a prediction
that this output file vouches for.

Only the participant identifier and the deposited group are projected.  Exact
stage, age, and sex stay structurally missing here as everywhere else.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

GROUPS = ("mild_f0_f1", "advanced_f3_f4")
EXPECTED_GROUP_COUNTS = {"mild_f0_f1": 40, "advanced_f3_f4": 32}
SERIES = "GSE49541"
COHORT_FAMILY_ID = "gse31803_gse49541_fibrosis_array"
PROJECTED_COLUMNS = ("participant_id", "fibrosis_stage_group")
NEVER_PROJECTED = ("exact_fibrosis_stage", "age", "sex", "nash_status")


class EvaluatorLabelError(RuntimeError):
    """Raised when the evaluator label table would break its fixed ordering."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_predictions_are_frozen(
    roots: Sequence[Path], expected: Sequence[str]
) -> dict[str, str]:
    if not roots or len(roots) != len(expected):
        raise EvaluatorLabelError("prediction roster and SHA roster differ")
    observed: dict[str, str] = {}
    for root, digest in zip(roots, expected, strict=True):
        manifest = root / "ARTIFACTS.json"
        complete = root / "COMPLETE"
        if not manifest.is_file() or not complete.is_file():
            raise EvaluatorLabelError(f"prediction bundle is not frozen: {root}")
        actual = sha256_file(manifest)
        if actual != digest:
            raise EvaluatorLabelError(f"prediction ARTIFACTS SHA differs: {root}")
        payload = json.loads(
            (root / "prediction_bundle.json").read_text(encoding="utf-8")
        )
        metadata = payload.get("metadata", {})
        if (
            metadata.get("external_labels_read") is not False
            or metadata.get("prediction_frozen_before_evaluator_label_join") is not True
        ):
            raise EvaluatorLabelError(f"prediction bundle is not label-blind: {root}")
        observed[str(payload["model_id"])] = actual
    return dict(sorted(observed.items()))


def build_labels(participants_path: Path) -> tuple[list[dict[str, str]], dict[str, Any]]:
    with participants_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise EvaluatorLabelError("participant table has no header")
        missing = sorted(set(PROJECTED_COLUMNS) - set(reader.fieldnames))
        if missing:
            raise EvaluatorLabelError(f"participant table lacks columns: {missing}")
        rows = [row for row in reader if row.get("series") == SERIES]
    if {row["cohort_family_id"] for row in rows} != {COHORT_FAMILY_ID}:
        raise EvaluatorLabelError("participant table carries a foreign cohort family")
    labels = [
        {
            "row_id": row["participant_id"],
            "fibrosis_stage_group": row["fibrosis_stage_group"],
        }
        for row in rows
    ]
    labels.sort(key=lambda row: row["row_id"])
    row_ids = [row["row_id"] for row in labels]
    if len(set(row_ids)) != len(row_ids):
        raise EvaluatorLabelError("a participant identifier repeats")
    if any(row["fibrosis_stage_group"] not in GROUPS for row in labels):
        raise EvaluatorLabelError("a deposited fibrosis group is off its roster")
    counts = Counter(row["fibrosis_stage_group"] for row in labels)
    if counts != Counter(EXPECTED_GROUP_COUNTS):
        raise EvaluatorLabelError(
            "deposited fibrosis-group census differs from 40 mild and 32 advanced"
        )
    # Structurally missing fields exist in the participant table but must not be
    # carried into the evaluator table under any name.
    for row in labels:
        if set(NEVER_PROJECTED) & set(row):
            raise EvaluatorLabelError("evaluator label table projected a forbidden column")
    receipt = {
        "schema_version": "masld-bench-gse49541-evaluator-labels-v1",
        "status": "pass_evaluator_only_outcome_materialized_after_prediction_freeze",
        "series": SERIES,
        "cohort_family_id": COHORT_FAMILY_ID,
        "participants": len(labels),
        "group_counts": {key: int(counts[key]) for key in GROUPS},
        "projected_columns": list(PROJECTED_COLUMNS),
        "structurally_missing_fields_projected": [],
        "structurally_missing_imputed": False,
        "outcome_is_evaluator_only": True,
        "outcome_visible_to_any_model": False,
    }
    return labels, receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--participants", required=True, type=Path)
    parser.add_argument("--prediction", action="append", required=True)
    parser.add_argument("--prediction-artifacts-sha256", action="append", required=True)
    parser.add_argument("--output-labels", required=True, type=Path)
    parser.add_argument("--output-receipt", required=True, type=Path)
    arguments = parser.parse_args()
    frozen = verify_predictions_are_frozen(
        [Path(value) for value in arguments.prediction],
        arguments.prediction_artifacts_sha256,
    )
    labels, receipt = build_labels(arguments.participants)
    receipt["prediction_bundle_artifacts_sha256_verified_before_label_write"] = frozen
    receipt["participants_table_sha256"] = sha256_file(arguments.participants)
    with arguments.output_labels.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("row_id", "fibrosis_stage_group"),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(labels)
    with arguments.output_receipt.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

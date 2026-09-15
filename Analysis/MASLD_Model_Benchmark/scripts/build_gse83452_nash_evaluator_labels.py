#!/usr/bin/env python3
"""Materialise the evaluator-only GSE83452 baseline NASH outcome after predictions freeze.

The deposited NASH status is an outcome and belongs to the evaluator alone.  This
step exists so the ordering is a fact on disk rather than a claim: it will not
emit a label table until it has verified the committed SHA-256 of every
already-frozen prediction bundle, so no label can have informed a prediction
that this output file vouches for.

Only the participant identifier and the deposited baseline status are projected.
Age, sex and intervention stay out of the evaluator table: age and sex belong to
a separate lane the training source cannot support, and intervention is a
transport stratum, not an outcome.

Four of the 152 baseline participants carry ``undefined``.  The TaskSpec's
missingness policy keeps undefined as undefined and excludes it from the NASH
endpoint; it is never mapped to no-NASH.  All 152 stay in the audit and 148 are
emitted as scorable.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

CLASSES = ("no_nash", "nash")
UNDEFINED = "undefined"
SERIES = "GSE83452"
COHORT_FAMILY_ID = "antwerp_inserm_shared"
BASELINE_TIMEPOINT = "baseline"
EXPECTED_BASELINE_RECORDS = 152
EXPECTED_ENDPOINT_EVALUABLE = 148
EXPECTED_CLASS_COUNTS = {"nash": 104, "no_nash": 44}
EXPECTED_UNDEFINED = 4
PROJECTED_COLUMNS = ("participant_id", "nash_status")
NEVER_PROJECTED = ("age", "sex", "intervention", "scan_date")


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
        # Keyed by path, not by model_id: this lane freezes one bundle per model
        # per arm, and keying by model_id would silently record five verifications
        # where ten happened.
        key = f"{root.parent.parent.name}/{payload['model_id']}"
        if key in observed:
            raise EvaluatorLabelError(f"prediction bundle is listed twice: {key}")
        observed[key] = actual
    return dict(sorted(observed.items()))


def build_labels(records_path: Path) -> tuple[list[dict[str, str]], dict[str, Any]]:
    with records_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise EvaluatorLabelError("record table has no header")
        missing = sorted(
            {"participant_id", "nash_status", "timepoint", "series", "cohort_family_id"}
            - set(reader.fieldnames)
        )
        if missing:
            raise EvaluatorLabelError(f"record table lacks columns: {missing}")
        rows = [row for row in reader if row.get("series") == SERIES]
    if {row["cohort_family_id"] for row in rows} != {COHORT_FAMILY_ID}:
        raise EvaluatorLabelError("record table carries a foreign cohort family")
    baseline = [row for row in rows if row["timepoint"] == BASELINE_TIMEPOINT]
    if len(baseline) != EXPECTED_BASELINE_RECORDS:
        raise EvaluatorLabelError("baseline record census differs")
    participants = [row["participant_id"] for row in baseline]
    if len(set(participants)) != len(participants):
        raise EvaluatorLabelError("a baseline participant identifier repeats")
    census = Counter(row["nash_status"] for row in baseline)
    if any(value not in set(CLASSES) | {UNDEFINED} for value in census):
        raise EvaluatorLabelError("a deposited NASH status is off its roster")
    if census.get(UNDEFINED, 0) != EXPECTED_UNDEFINED:
        raise EvaluatorLabelError("undefined NASH census differs")
    labels = [
        {"row_id": row["participant_id"], "nash_status": row["nash_status"]}
        for row in baseline
        if row["nash_status"] != UNDEFINED
    ]
    labels.sort(key=lambda row: row["row_id"])
    if len(labels) != EXPECTED_ENDPOINT_EVALUABLE:
        raise EvaluatorLabelError("endpoint-evaluable census differs")
    counts = Counter(row["nash_status"] for row in labels)
    if dict(counts) != EXPECTED_CLASS_COUNTS:
        raise EvaluatorLabelError(
            "deposited NASH census differs from 104 NASH and 44 no NASH"
        )
    for row in labels:
        if set(NEVER_PROJECTED) & set(row):
            raise EvaluatorLabelError("evaluator label table projected a forbidden column")
    receipt = {
        "schema_version": "masld-bench-gse83452-evaluator-labels-v1",
        "status": "pass_evaluator_only_outcome_materialized_after_prediction_freeze",
        "series": SERIES,
        "cohort_family_id": COHORT_FAMILY_ID,
        "timepoint": BASELINE_TIMEPOINT,
        "baseline_records": len(baseline),
        "endpoint_evaluable": len(labels),
        "class_counts": {key: int(counts[key]) for key in CLASSES},
        "undefined_retained_as_missing": int(census.get(UNDEFINED, 0)),
        "undefined_mapped_to_no_nash": False,
        "projected_columns": list(PROJECTED_COLUMNS),
        "never_projected_columns": list(NEVER_PROJECTED),
        "structurally_missing_imputed": False,
        "outcome_is_evaluator_only": True,
        "outcome_visible_to_any_model": False,
    }
    return labels, receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", required=True, type=Path)
    parser.add_argument("--prediction", action="append", required=True)
    parser.add_argument("--prediction-artifacts-sha256", action="append", required=True)
    parser.add_argument("--output-labels", required=True, type=Path)
    parser.add_argument("--output-receipt", required=True, type=Path)
    arguments = parser.parse_args()
    frozen = verify_predictions_are_frozen(
        [Path(value) for value in arguments.prediction],
        arguments.prediction_artifacts_sha256,
    )
    labels, receipt = build_labels(arguments.records)
    receipt["prediction_bundle_artifacts_sha256_verified_before_label_write"] = frozen
    receipt["records_table_sha256"] = sha256_file(arguments.records)
    with arguments.output_labels.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("row_id", "nash_status"),
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

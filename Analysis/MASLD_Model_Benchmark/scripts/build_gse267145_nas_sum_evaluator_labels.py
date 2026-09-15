#!/usr/bin/env python3
"""Materialise the evaluator-only GSE267145 activity sum after predictions freeze.

Ordering is a fact on disk rather than a claim: no label table is written until
every frozen prediction bundle's committed SHA-256 has been verified.  Only the
participant identifier and the deposited component sum are projected; the
individual components, fibrosis, stage and recorded sex stay with the source.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

PROJECTED = ("participant_id", "nash_crn_component_sum")
NEVER_PROJECTED = ("steatosis", "ballooning", "lobular_inflammation",
                   "lobular_necrosis", "fibrosis", "stage3", "stage5", "recorded_sex")


class EvaluatorLabelError(RuntimeError):
    """Raised when the label table would break its fixed ordering."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_predictions_are_frozen(roots: Sequence[Path], expected: Sequence[str]) -> dict[str, str]:
    if not roots or len(roots) != len(expected):
        raise EvaluatorLabelError("prediction roster and SHA roster differ")
    observed: dict[str, str] = {}
    for root, digest in zip(roots, expected, strict=True):
        if not (root / "ARTIFACTS.json").is_file() or not (root / "COMPLETE").is_file():
            raise EvaluatorLabelError(f"prediction bundle is not frozen: {root}")
        actual = sha256_file(root / "ARTIFACTS.json")
        if actual != digest:
            raise EvaluatorLabelError(f"prediction ARTIFACTS SHA differs: {root}")
        payload = json.loads((root / "prediction_bundle.json").read_text(encoding="utf-8"))
        meta = payload.get("metadata", {})
        if (
            meta.get("target_outcomes_read") is not False
            or meta.get("prediction_frozen_before_evaluator_label_join") is not True
        ):
            raise EvaluatorLabelError(f"prediction bundle is not outcome-blind: {root}")
        observed[str(payload["model_id"])] = actual
    return dict(sorted(observed.items()))


def build(path: Path, expected_participants: int) -> tuple[list[dict[str, str]], dict[str, Any]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise EvaluatorLabelError("endpoint table has no header")
        missing = sorted(set(PROJECTED) - set(reader.fieldnames))
        if missing:
            raise EvaluatorLabelError(f"endpoint table lacks {missing}")
        rows = list(reader)
    labels = [{k: row[k] for k in PROJECTED} for row in rows]
    labels = [{"row_id": r["participant_id"],
               "nash_crn_component_sum": r["nash_crn_component_sum"]} for r in labels]
    labels.sort(key=lambda r: r["row_id"])
    if len(labels) != expected_participants:
        raise EvaluatorLabelError("participant census differs")
    if len({r["row_id"] for r in labels}) != len(labels):
        raise EvaluatorLabelError("a participant identifier repeats")
    values = [int(r["nash_crn_component_sum"]) for r in labels]
    if any(v < 0 or v > 8 for v in values):
        raise EvaluatorLabelError("an activity sum is off its 0-to-8 scale")
    for row in labels:
        if set(NEVER_PROJECTED) & set(row):
            raise EvaluatorLabelError("label table projected a forbidden column")
    census: dict[str, int] = {}
    for v in values:
        census[str(v)] = census.get(str(v), 0) + 1
    receipt = {
        "schema_version": "masld-bench-gse267145-nas-evaluator-labels-v1",
        "status": "pass_evaluator_only_outcome_materialized_after_prediction_freeze",
        "series": "GSE267145",
        "participants": len(labels),
        "activity_sum_census": {k: census[k] for k in sorted(census, key=int)},
        "activity_sum_mean": sum(values) / len(values),
        "projected_columns": ["row_id", "nash_crn_component_sum"],
        "components_projected": False,
        "outcome_is_evaluator_only": True,
        "outcome_visible_to_any_model": False,
    }
    return labels, receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endpoints", required=True, type=Path)
    parser.add_argument("--expected-participants", required=True, type=int)
    parser.add_argument("--prediction", action="append", required=True)
    parser.add_argument("--prediction-artifacts-sha256", action="append", required=True)
    parser.add_argument("--output-labels", required=True, type=Path)
    parser.add_argument("--output-receipt", required=True, type=Path)
    a = parser.parse_args()
    frozen = verify_predictions_are_frozen(
        [Path(v) for v in a.prediction], a.prediction_artifacts_sha256)
    labels, receipt = build(a.endpoints, a.expected_participants)
    receipt["prediction_bundle_artifacts_sha256_verified_before_label_write"] = frozen
    receipt["endpoints_table_sha256"] = sha256_file(a.endpoints)
    with a.output_labels.open("x", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=("row_id", "nash_crn_component_sum"),
                           delimiter="\t", lineterminator="\n", extrasaction="raise")
        w.writeheader()
        w.writerows(labels)
    with a.output_receipt.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

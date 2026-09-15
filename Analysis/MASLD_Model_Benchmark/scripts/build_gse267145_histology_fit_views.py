#!/usr/bin/env python3
"""Build fold-scoped GSE267145 training labels without held outcomes or sex."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


class HistologyFitViewError(RuntimeError):
    """Raised when a fit view would cross the outer-test outcome separation."""


TRAINING_FIELDS = (
    "participant_id",
    "inner_validation_fold",
    "stage3",
    "nash_crn_component_sum",
    "fibrosis",
    "fibrosis_group3",
)
QUERY_FIELDS = ("participant_id", "outer_fold")
INNER_FIELDS = ("outer_fold", "participant_id", "inner_validation_fold")


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise HistologyFitViewError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def fibrosis_group(value: int) -> str:
    if value == 0:
        return "F0"
    if value == 1:
        return "F1"
    if value in {2, 3}:
        return "F2_3"
    raise HistologyFitViewError("fibrosis is outside the frozen source scale")


def build_fit_views(
    *,
    participant_axis: Path,
    outer_folds: Path,
    endpoints: Path,
    inner_roster: Path,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise HistologyFitViewError(f"refusing to overwrite fit views: {output}")
    participant_fields, participants = read_tsv(participant_axis)
    fold_fields, folds = read_tsv(outer_folds)
    endpoint_fields, endpoint_rows = read_tsv(endpoints)
    inner_fields, inner_rows = read_tsv(inner_roster)
    if (
        "participant_id" not in participant_fields
        or fold_fields != ("participant_id", "outer_fold")
        or inner_fields != INNER_FIELDS
        or not {
            "participant_id",
            "outer_fold",
            "stage3",
            "nash_crn_component_sum",
            "fibrosis",
        }
        <= set(endpoint_fields)
    ):
        raise HistologyFitViewError("input schemas differ")
    participant_ids = [row["participant_id"] for row in participants]
    if (
        len(participant_ids) != 99
        or len(set(participant_ids)) != 99
        or participant_ids != [row["participant_id"] for row in folds]
        or participant_ids != [row["participant_id"] for row in endpoint_rows]
        or [row["outer_fold"] for row in folds]
        != [row["outer_fold"] for row in endpoint_rows]
    ):
        raise HistologyFitViewError("participant, fold, and endpoint axes differ")
    inner_lookup: dict[tuple[int, str], int] = {}
    for row in inner_rows:
        key = (int(row["outer_fold"]), row["participant_id"])
        if key in inner_lookup:
            raise HistologyFitViewError("inner participant assignment is duplicated")
        inner_lookup[key] = int(row["inner_validation_fold"])
    output.mkdir(parents=True)
    endpoint_by_id = {row["participant_id"]: row for row in endpoint_rows}
    receipts: list[dict[str, Any]] = []
    for outer_fold in range(5):
        fold_root = output / f"outer_{outer_fold}"
        fold_root.mkdir()
        training_ids = [
            row["participant_id"] for row in folds if int(row["outer_fold"]) != outer_fold
        ]
        query_ids = [
            row["participant_id"] for row in folds if int(row["outer_fold"]) == outer_fold
        ]
        training: list[dict[str, Any]] = []
        for participant in training_ids:
            row = endpoint_by_id[participant]
            key = (outer_fold, participant)
            if key not in inner_lookup or inner_lookup[key] not in range(4):
                raise HistologyFitViewError("inner fit assignment is missing or invalid")
            fibrosis = int(row["fibrosis"])
            nas = int(row["nash_crn_component_sum"])
            if row["stage3"] not in {"NOR", "NAFL", "NASH"} or nas not in range(9):
                raise HistologyFitViewError("training endpoint differs")
            training.append(
                {
                    "participant_id": participant,
                    "inner_validation_fold": inner_lookup[key],
                    "stage3": row["stage3"],
                    "nash_crn_component_sum": nas,
                    "fibrosis": fibrosis,
                    "fibrosis_group3": fibrosis_group(fibrosis),
                }
            )
        query = [
            {"participant_id": participant, "outer_fold": outer_fold}
            for participant in query_ids
        ]
        if set(training_ids) & set(query_ids) or len(training_ids) + len(query_ids) != 99:
            raise HistologyFitViewError("outer training/query partition differs")
        write_tsv(fold_root / "training_endpoints.tsv", TRAINING_FIELDS, training)
        write_tsv(fold_root / "query_participants.tsv", QUERY_FIELDS, query)
        receipt = {
            "outer_fold": outer_fold,
            "training_participants": len(training),
            "query_participants": len(query),
            "training_stage3": dict(sorted(Counter(row["stage3"] for row in training).items())),
            "training_fibrosis": dict(
                sorted(Counter(str(row["fibrosis"]) for row in training).items())
            ),
            "inner_fold_counts": {
                str(key): value
                for key, value in sorted(
                    Counter(row["inner_validation_fold"] for row in training).items()
                )
            },
            "held_endpoint_values_included": False,
            "recorded_sex_included": False,
            "source_stage5_included": False,
        }
        (fold_root / "receipt.json").write_text(
            json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        receipts.append(receipt)
    summary = {
        "schema_version": "masld-bench-gse267145-histology-fit-views-v1",
        "status": "passed",
        "participants": 99,
        "outer_folds": 5,
        "training_label_fields": list(TRAINING_FIELDS),
        "query_fields": list(QUERY_FIELDS),
        "held_endpoint_values_included": False,
        "recorded_sex_included": False,
        "source_stage5_included": False,
        "participant_age_state": "structurally_missing",
        "folds": receipts,
    }
    (output / "receipt.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecular", required=True, type=Path)
    parser.add_argument("--molecular-artifacts-sha256", required=True)
    parser.add_argument("--folds", required=True, type=Path)
    parser.add_argument("--folds-artifacts-sha256", required=True)
    parser.add_argument("--outcomes", required=True, type=Path)
    parser.add_argument("--outcomes-artifacts-sha256", required=True)
    parser.add_argument("--inner-validation", required=True, type=Path)
    parser.add_argument("--inner-validation-artifacts-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    for root, expected in (
        (arguments.molecular, arguments.molecular_artifacts_sha256),
        (arguments.folds, arguments.folds_artifacts_sha256),
        (arguments.outcomes, arguments.outcomes_artifacts_sha256),
        (arguments.inner_validation, arguments.inner_validation_artifacts_sha256),
    ):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise HistologyFitViewError(f"input ARTIFACTS SHA differs: {root}")
    result = build_fit_views(
        participant_axis=arguments.molecular / "participant_axis.tsv",
        outer_folds=arguments.folds / "participant_outer_folds.tsv",
        endpoints=arguments.outcomes / "participant_endpoints.tsv",
        inner_roster=arguments.inner_validation / "inner_fold_roster.tsv",
        output=arguments.output,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

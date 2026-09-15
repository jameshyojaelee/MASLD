#!/usr/bin/env python3
"""Freeze study-held-out outer folds for the 50,000-cell Atlas screen."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import canonical_sha256, sha256_file


DATASET_VIEW_ID = "resource_atlas_frozen_screen_50000_v1"
SPLIT_ID = "resource_atlas_study_outer_5fold_v1"
INNER_SEED = 20260824
EXPECTED_ROWS = 50_000
EXPECTED_DONORS = 102
EXPECTED_DATASET_COUNTS = {
    "GSE136103": 4_433,
    "GSE174748": 2_538,
    "GSE185477": 2_026,
    "GSE189600": 149,
    "GSE202379": 21_572,
    "GSE244832": 12_241,
    "Liver_Atlas": 7_041,
}
EXPECTED_DATASET_DONORS = {
    "GSE136103": 10,
    "GSE174748": 4,
    "GSE185477": 3,
    "GSE189600": 2,
    "GSE202379": 46,
    "GSE244832": 18,
    "Liver_Atlas": 19,
}
EXPECTED_LABELS = (
    "cholangiocyte",
    "endothelial",
    "hepatocyte",
    "immune",
    "mesenchymal_stromal",
)
REQUIRED_SELECTION_FIELDS = {
    "row_id",
    "donor_id",
    "dataset",
    "broad_label",
}


class StudyFoldError(ValueError):
    """Raised when a source or split does not meet the frozen requirements."""


def assign_outer_folds(dataset_counts: Mapping[str, int]) -> dict[str, int]:
    """Give the four largest studies independent folds and group the rest."""

    if len(dataset_counts) < 5 or any(not key or value < 1 for key, value in dataset_counts.items()):
        raise StudyFoldError("outer study assignment requires at least five nonempty studies")
    ranked = sorted(dataset_counts, key=lambda key: (-dataset_counts[key], key))
    return {
        dataset: rank if rank < 4 else 4
        for rank, dataset in enumerate(ranked)
    }


def assign_inner_donors(
    donors: Iterable[str], *, outer_fold: int, seed: int = INNER_SEED
) -> dict[str, int]:
    """Assign training donors to five near-balanced, label-free inner folds."""

    unique = sorted(set(donors))
    if len(unique) < 5 or not 0 <= outer_fold < 5:
        raise StudyFoldError("inner assignment requires five donors and a valid outer fold")
    ranked = sorted(
        unique,
        key=lambda donor: hashlib.sha256(
            f"{seed}\0outer={outer_fold}\0inner\0{donor}".encode("utf-8")
        ).digest(),
    )
    return {donor: rank % 5 for rank, donor in enumerate(ranked)}


def _read_selection(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not REQUIRED_SELECTION_FIELDS.issubset(
            reader.fieldnames
        ):
            raise StudyFoldError("selection table lacks required split fields")
        rows = list(reader)
    if not rows:
        raise StudyFoldError("selection table is empty")
    return rows


def _write_tsv(path: Path, fields: tuple[str, ...], rows: Iterable[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def build(
    *,
    source: Path,
    output: Path,
    expected_source_artifacts_sha256: str,
    expected_rows: int = EXPECTED_ROWS,
    expected_donors: int = EXPECTED_DONORS,
    expected_dataset_counts: Mapping[str, int] = EXPECTED_DATASET_COUNTS,
    expected_dataset_donors: Mapping[str, int] = EXPECTED_DATASET_DONORS,
) -> None:
    if output.exists():
        raise StudyFoldError("refusing to overwrite split output")
    if sha256_file(source / "ARTIFACTS.json") != expected_source_artifacts_sha256:
        raise StudyFoldError("source ARTIFACTS SHA-256 differs")
    source_manifest = verify_frozen_tree(source)
    source_metadata = source_manifest.get("metadata", {})
    if (
        source_metadata.get("artifact_class")
        != "resource_atlas_frozen_screen_fixture"
        or source_metadata.get("subset_id") != DATASET_VIEW_ID
        or source_metadata.get("row_count") != expected_rows
        or source_metadata.get("sealed_outcomes_read") is not False
        or source_metadata.get("histology_read") is not False
    ):
        raise StudyFoldError("source frozen-screen metadata differs")

    rows = _read_selection(source / "selection.tsv")
    if len(rows) != expected_rows:
        raise StudyFoldError("selection row count differs")
    row_ids = [row["row_id"] for row in rows]
    if any(not row_id for row_id in row_ids) or len(set(row_ids)) != len(row_ids):
        raise StudyFoldError("row identities are empty or duplicated")

    dataset_counts = Counter(row["dataset"] for row in rows)
    donor_datasets: dict[str, set[str]] = defaultdict(set)
    dataset_donors: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        donor = row["donor_id"]
        dataset = row["dataset"]
        if not donor or not dataset:
            raise StudyFoldError("donor and dataset identities must be nonempty")
        donor_datasets[donor].add(dataset)
        dataset_donors[dataset].add(donor)
    if any(len(values) != 1 for values in donor_datasets.values()):
        raise StudyFoldError("one donor identity occurs in multiple studies")
    observed_dataset_donors = {
        dataset: len(donors) for dataset, donors in dataset_donors.items()
    }
    if dict(dataset_counts) != dict(expected_dataset_counts):
        raise StudyFoldError("study cell counts differ from the frozen contract")
    if observed_dataset_donors != dict(expected_dataset_donors):
        raise StudyFoldError("study donor counts differ from the frozen contract")
    if len(donor_datasets) != expected_donors:
        raise StudyFoldError("total donor count differs")

    outer = assign_outer_folds(dataset_counts)
    donors_outer = {
        donor: outer[next(iter(datasets))]
        for donor, datasets in donor_datasets.items()
    }
    if set(outer.values()) != set(range(5)):
        raise StudyFoldError("outer folds are incomplete")

    output.mkdir(mode=0o750)
    ranked_datasets = sorted(dataset_counts, key=lambda key: (-dataset_counts[key], key))
    study_rows = [
        {
            "dataset_rank_by_cell_count": rank,
            "dataset": dataset,
            "cells": dataset_counts[dataset],
            "donors": len(dataset_donors[dataset]),
            "outer_fold": outer[dataset],
        }
        for rank, dataset in enumerate(ranked_datasets)
    ]
    _write_tsv(
        output / "study_outer_folds.tsv",
        ("dataset_rank_by_cell_count", "dataset", "cells", "donors", "outer_fold"),
        study_rows,
    )

    split_rows = [
        {
            "row_id": row["row_id"],
            "donor_id": row["donor_id"],
            "dataset": row["dataset"],
            "outer_fold": outer[row["dataset"]],
        }
        for row in rows
    ]
    _write_tsv(
        output / "row_outer_folds.tsv",
        ("row_id", "donor_id", "dataset", "outer_fold"),
        split_rows,
    )

    inner_rows: list[dict[str, Any]] = []
    for held_outer in range(5):
        training_donors = sorted(
            donor for donor, fold in donors_outer.items() if fold != held_outer
        )
        inner = assign_inner_donors(training_donors, outer_fold=held_outer)
        for donor in training_donors:
            inner_rows.append(
                {
                    "held_outer_fold": held_outer,
                    "donor_id": donor,
                    "dataset": next(iter(donor_datasets[donor])),
                    "inner_fold": inner[donor],
                }
            )
    _write_tsv(
        output / "inner_donor_folds.tsv",
        ("held_outer_fold", "donor_id", "dataset", "inner_fold"),
        inner_rows,
    )

    coverage_rows: list[dict[str, Any]] = []
    coverage_complete = True
    for held_outer in range(5):
        held = [row for row in rows if outer[row["dataset"]] == held_outer]
        for label in EXPECTED_LABELS:
            selected = [row for row in held if row["broad_label"] == label]
            coverage_rows.append(
                {
                    "outer_fold": held_outer,
                    "broad_label": label,
                    "cells": len(selected),
                    "donors": len({row["donor_id"] for row in selected}),
                }
            )
            coverage_complete &= bool(selected)
    _write_tsv(
        output / "post_assignment_class_coverage.tsv",
        ("outer_fold", "broad_label", "cells", "donors"),
        coverage_rows,
    )
    if not coverage_complete:
        raise StudyFoldError("a frozen outer fold lacks an evaluation class")

    manifest = {
        "schema_version": "masld-bench-study-folds-v1",
        "split_id": SPLIT_ID,
        "dataset_view_id": DATASET_VIEW_ID,
        "source": {
            "path": source.resolve(strict=True).as_posix(),
            "artifacts_sha256": expected_source_artifacts_sha256,
        },
        "biological_unit": "donor",
        "outer_assignment": {
            "folds": 5,
            "policy": "four_largest_studies_independent_then_all_remaining_studies_grouped",
            "ranking_keys": ["descending_cell_count", "ascending_dataset_id_tiebreak"],
            "target_labels_used": False,
            "histology_used": False,
            "sealed_outcomes_used": False,
            "mapping": dict(sorted(outer.items())),
        },
        "inner_assignment": {
            "folds": 5,
            "policy": "per_outer_training_donor_sha256_rank_then_round_robin",
            "seed": INNER_SEED,
            "target_labels_used": False,
        },
        "post_assignment_validation": {
            "target_labels_read": True,
            "target_labels_used_to_change_assignment": False,
            "all_five_classes_in_every_outer_fold": coverage_complete,
        },
        "counts": {
            "rows": len(rows),
            "donors": len(donor_datasets),
            "studies": len(dataset_counts),
        },
    }
    manifest["split_sha256"] = canonical_sha256(manifest)
    write_json_exclusive(output / "split_manifest.json", manifest)
    freeze_tree(
        output,
        {
            "artifact_class": "cell_state_study_outer_split",
            "dataset_view_id": DATASET_VIEW_ID,
            "split_id": SPLIT_ID,
            "rows": len(rows),
            "donors": len(donor_datasets),
            "studies": len(dataset_counts),
            "target_labels_used_for_assignment": False,
            "target_labels_read_for_post_assignment_validation": True,
            "histology_read": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--source", required=True, type=Path)
    value.add_argument("--output", required=True, type=Path)
    value.add_argument("--expected-source-artifacts-sha256", required=True)
    return value


def main() -> int:
    arguments = parser().parse_args()
    build(
        source=arguments.source,
        output=arguments.output,
        expected_source_artifacts_sha256=arguments.expected_source_artifacts_sha256,
    )
    print(json.dumps({"output": arguments.output.resolve().as_posix(), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Assemble outcome-safe ChromBPNet train, validation, and test regions."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Mapping, Sequence


PRIMARY_CONTIGS = tuple([f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY"])
CCRE_FIELDS = (
    "contig",
    "output_start",
    "output_end",
    "window_id",
    "genomic_fold",
    "window_class",
    "ccre_class",
    "ccre_id",
    "ccre_start",
    "ccre_end",
    "input_start",
    "input_end",
    "selection_hash",
)
SPLIT_FIELDS = (
    "split_id",
    "donor_train_folds",
    "donor_valid_fold",
    "donor_test_fold",
    "genomic_train_folds",
    "genomic_valid_fold",
    "genomic_test_fold",
)
GENOMIC_FOLD_FIELDS = ("genomic_fold", "contigs", "total_bp")


class FoldRegionError(ValueError):
    """Raised when a fold-region assembly does not meet its leakage requirements."""


def read_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise FoldRegionError(f"{path.name} fields differ")
        return [dict(row) for row in reader]


def read_narrowpeak(path: Path) -> list[list[str]]:
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 10:
                raise FoldRegionError(
                    f"training peak row {line_number} does not have ten fields"
                )
            try:
                start, end, summit = int(fields[1]), int(fields[2]), int(fields[9])
            except ValueError as error:
                raise FoldRegionError("training peak geometry is not integer-valued") from error
            if start < 0 or end <= start or summit < 0:
                raise FoldRegionError("training peak geometry differs")
            rows.append(fields)
    if not rows:
        raise FoldRegionError("training peak file is empty")
    return rows


def fold_contract(
    split_contract: Path, split_id: str
) -> tuple[dict[str, str], dict[int, tuple[str, ...]], dict[str, int]]:
    matches = [
        row
        for row in read_tsv(split_contract / "crossed_outer_splits.tsv", SPLIT_FIELDS)
        if row["split_id"] == split_id
    ]
    if len(matches) != 1:
        raise FoldRegionError("crossed split identifier is not unique")
    row = matches[0]
    contigs_by_fold = {
        int(item["genomic_fold"]): tuple(item["contigs"].split(","))
        for item in read_tsv(
            split_contract / "genomic_folds.tsv", GENOMIC_FOLD_FIELDS
        )
    }
    if set(contigs_by_fold) != set(range(5)):
        raise FoldRegionError("genomic fold roster differs")
    fold_by_contig = {
        contig: fold for fold, contigs in contigs_by_fold.items() for contig in contigs
    }
    if (
        set(fold_by_contig) != set(PRIMARY_CONTIGS)
        or sum(len(contigs) for contigs in contigs_by_fold.values())
        != len(PRIMARY_CONTIGS)
    ):
        raise FoldRegionError("primary contig roster differs")
    return row, contigs_by_fold, fold_by_contig


def fixed_ccre_narrowpeak(row: Mapping[str, str], role: str) -> list[str]:
    start, end = int(row["output_start"]), int(row["output_end"])
    if end - start != 1000:
        raise FoldRegionError("fixed cCRE output window length differs")
    return [
        row["contig"],
        str(start),
        str(end),
        f"fixed_{role}|{row['window_id']}|{row['ccre_class']}",
        "0",
        ".",
        "0",
        "-1",
        "-1",
        "500",
    ]


def assemble(
    *,
    training_peaks: Path,
    split_contract: Path,
    split_id: str,
    output_peaks: Path,
    output_fold: Path,
    output_summary: Path,
    expected_ccre_per_fold: int = 16_000,
) -> dict[str, object]:
    if any(path.exists() for path in (output_peaks, output_fold, output_summary)):
        raise FoldRegionError("refusing to overwrite fold-region output")
    split, contigs_by_fold, fold_by_contig = fold_contract(split_contract, split_id)
    train_folds = tuple(int(value) for value in split["genomic_train_folds"].split(","))
    valid_fold = int(split["genomic_valid_fold"])
    test_fold = int(split["genomic_test_fold"])
    if (
        len(train_folds) != 3
        or set(train_folds) | {valid_fold, test_fold} != set(range(5))
        or len(set(train_folds) | {valid_fold, test_fold}) != 5
    ):
        raise FoldRegionError("crossed genomic split does not partition five folds")

    assembled: list[tuple[str, list[str]]] = []
    seen_centers: set[tuple[str, int]] = set()
    train_rows = read_narrowpeak(training_peaks)
    for fields in train_rows:
        contig = fields[0]
        if fold_by_contig.get(contig) not in train_folds:
            raise FoldRegionError("MACS2 training peak is outside genomic training folds")
        center = int(fields[1]) + int(fields[9])
        key = (contig, center)
        if key in seen_centers:
            raise FoldRegionError("training peak center is duplicated")
        seen_centers.add(key)
        fields = list(fields)
        fields[3] = f"train_macs2|{fields[3]}"
        assembled.append(("train", fields))

    ccre_rows = read_tsv(
        split_contract / "ccre_evaluation_windows.tsv", CCRE_FIELDS
    )
    ccre_counts = {fold: 0 for fold in range(5)}
    role_counts = {"train": len(train_rows), "valid": 0, "test": 0}
    role_class_counts: dict[str, dict[str, int]] = {"valid": {}, "test": {}}
    for row in ccre_rows:
        fold = int(row["genomic_fold"])
        if fold_by_contig.get(row["contig"]) != fold:
            raise FoldRegionError("cCRE row fold differs from its contig fold")
        ccre_counts[fold] += 1
        if fold not in {valid_fold, test_fold}:
            continue
        role = "valid" if fold == valid_fold else "test"
        fields = fixed_ccre_narrowpeak(row, role)
        center = int(fields[1]) + int(fields[9])
        key = (fields[0], center)
        if key in seen_centers:
            raise FoldRegionError("assembled region center is duplicated")
        seen_centers.add(key)
        assembled.append((role, fields))
        role_counts[role] += 1
        class_counts = role_class_counts[role]
        class_counts[row["ccre_class"]] = class_counts.get(row["ccre_class"], 0) + 1
    if any(value != expected_ccre_per_fold for value in ccre_counts.values()):
        raise FoldRegionError(f"cCRE fold census differs: {ccre_counts!r}")
    if role_counts["valid"] != expected_ccre_per_fold or role_counts[
        "test"
    ] != expected_ccre_per_fold:
        raise FoldRegionError("fixed held-region census differs")

    contig_rank = {contig: index for index, contig in enumerate(PRIMARY_CONTIGS)}
    assembled.sort(
        key=lambda item: (
            contig_rank[item[1][0]],
            int(item[1][1]) + int(item[1][9]),
            item[0],
            item[1][3],
        )
    )
    with output_peaks.open("x", encoding="utf-8", newline="") as handle:
        for _role, fields in assembled:
            handle.write("\t".join(fields) + "\n")
    fold_payload = {
        "train": [contig for fold in train_folds for contig in contigs_by_fold[fold]],
        "valid": list(contigs_by_fold[valid_fold]),
        "test": list(contigs_by_fold[test_fold]),
    }
    output_fold.write_text(
        json.dumps(fold_payload, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    summary = {
        "schema_version": "masld-bench-chrombpnet-fold-regions-v1",
        "split_id": split_id,
        "donor_train_folds": split["donor_train_folds"],
        "donor_valid_fold": int(split["donor_valid_fold"]),
        "donor_test_fold": int(split["donor_test_fold"]),
        "genomic_train_folds": list(train_folds),
        "genomic_valid_fold": valid_fold,
        "genomic_test_fold": test_fold,
        "role_counts": role_counts,
        "held_role_ccre_class_counts": role_class_counts,
        "training_region_source": "MACS2_training_donor_and_genomic_folds_only",
        "validation_region_source": "fixed_ENCODE_v4_cCRE_contract",
        "test_region_source": "fixed_ENCODE_v4_cCRE_contract",
        "held_atac_used_for_region_selection": False,
        "model_training_may_use_test_role": False,
        "assembled_regions": len(assembled),
    }
    output_summary.write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--training-peaks", type=Path, required=True)
    parser.add_argument("--split-contract", type=Path, required=True)
    parser.add_argument("--split-id", required=True)
    parser.add_argument("--output-peaks", type=Path, required=True)
    parser.add_argument("--output-fold", type=Path, required=True)
    parser.add_argument("--output-summary", type=Path, required=True)
    arguments = parser.parse_args()
    result = assemble(
        training_peaks=arguments.training_peaks,
        split_contract=arguments.split_contract,
        split_id=arguments.split_id,
        output_peaks=arguments.output_peaks,
        output_fold=arguments.output_fold,
        output_summary=arguments.output_summary,
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

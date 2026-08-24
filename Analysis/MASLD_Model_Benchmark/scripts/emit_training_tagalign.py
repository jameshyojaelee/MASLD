#!/usr/bin/env python3
"""Emit donor-safe, genomic-test-excluded Tn5 tagAlign records."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from scripts.build_training_fold_pseudobulk import (
    close_deterministic_gzip,
    open_deterministic_gzip,
    parse_fragment_line,
    sha256_file,
)


PSEUDOBULK_FIELDS = (
    "donor_test_fold",
    "donor_valid_fold",
    "donor_train_folds",
    "lineage_id",
    "donors",
    "nuclei",
    "unique_fragments",
    "read_support",
    "tn5_insertions",
    "nonzero_positions",
    "max_pending_positions",
    "fragment_path",
    "fragment_size_bytes",
    "fragment_sha256",
    "bigwig_path",
    "bigwig_size_bytes",
    "bigwig_sha256",
)


class TagAlignError(ValueError):
    """Raised when a train-only tagAlign violates its contract."""


def read_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise TagAlignError(f"{path.name} fields differ")
        return [dict(row) for row in reader]


def genomic_fold_map(split_contract: Path) -> dict[str, int]:
    rows = read_tsv(
        split_contract / "genomic_folds.tsv",
        ("genomic_fold", "contigs", "total_bp"),
    )
    mapping = {
        contig: int(row["genomic_fold"])
        for row in rows
        for contig in row["contigs"].split(",")
    }
    if len(rows) != 5 or set(mapping.values()) != set(range(5)):
        raise TagAlignError("genomic fold contract differs")
    return mapping


def emit(
    *,
    source: Path,
    source_size_bytes: int,
    source_sha256: str,
    fold_by_contig: Mapping[str, int],
    genomic_test_fold: int,
    genomic_valid_fold: int,
    output: Path,
) -> dict[str, object]:
    if (
        output.exists()
        or not 0 <= genomic_test_fold < 5
        or not 0 <= genomic_valid_fold < 5
        or genomic_test_fold == genomic_valid_fold
    ):
        raise TagAlignError("invalid output or genomic fold contract")
    if source.stat().st_size != source_size_bytes or sha256_file(source) != source_sha256:
        raise TagAlignError("pseudobulk fragment source differs")
    text, hashed, raw = open_deterministic_gzip(output)
    included_fragments = 0
    included_read_support = 0
    excluded_fragments = {"test": 0, "valid": 0}
    excluded_read_support = {"test": 0, "valid": 0}
    previous = None
    with gzip.open(source, "rt", encoding="utf-8", newline="") as handle:
        try:
            for line in handle:
                key, fields = parse_fragment_line(line)
                if previous is not None and key < previous:
                    raise TagAlignError("pseudobulk fragment source is not sorted")
                previous = key
                contig, start_text, end_text, cell_id, support_text = fields
                support = int(support_text)
                fold = fold_by_contig.get(contig)
                if fold in {genomic_test_fold, genomic_valid_fold}:
                    role = "test" if fold == genomic_test_fold else "valid"
                    excluded_fragments[role] += 1
                    excluded_read_support[role] += support
                    continue
                if contig not in fold_by_contig:
                    raise TagAlignError("fragment contig is outside the fold contract")
                start, end = int(start_text), int(end_text)
                text.write(
                    "\t".join(
                        (contig, start_text, str(start + 1), f"{cell_id}/L", "1000", "+")
                    )
                    + "\n"
                )
                text.write(
                    "\t".join(
                        (contig, str(end - 1), end_text, f"{cell_id}/R", "1000", "-")
                    )
                    + "\n"
                )
                included_fragments += 1
                included_read_support += support
            output_sha256 = close_deterministic_gzip(text, hashed, raw)
        finally:
            if not text.closed:
                text.close()
                raw.close()
    return {
        "included_unique_fragments": included_fragments,
        "included_read_support": included_read_support,
        "included_tn5_events": 2 * included_fragments,
        "excluded_genomic_test_fragments": excluded_fragments["test"],
        "excluded_genomic_test_read_support": excluded_read_support["test"],
        "excluded_genomic_valid_fragments": excluded_fragments["valid"],
        "excluded_genomic_valid_read_support": excluded_read_support["valid"],
        "genomic_test_fold": genomic_test_fold,
        "genomic_valid_fold": genomic_valid_fold,
        "allowed_genomic_folds": [
            value
            for value in range(5)
            if value not in {genomic_test_fold, genomic_valid_fold}
        ],
        "read_support_used_as_signal_weight": False,
        "left_insertion_position": "fragment_start",
        "right_insertion_position": "fragment_end_minus_1_BED_exclusive",
        "output_sha256": output_sha256,
        "output_size_bytes": output.stat().st_size,
    }


def emit_by_genomic_fold(
    *,
    source: Path,
    source_size_bytes: int,
    source_sha256: str,
    fold_by_contig: Mapping[str, int],
    genomic_test_fold: int,
    genomic_valid_fold: int,
    output_directory: Path,
) -> dict[str, object]:
    if (
        output_directory.exists()
        or not 0 <= genomic_test_fold < 5
        or not 0 <= genomic_valid_fold < 5
        or genomic_test_fold == genomic_valid_fold
    ):
        raise TagAlignError("invalid output or genomic fold contract")
    if source.stat().st_size != source_size_bytes or sha256_file(source) != source_sha256:
        raise TagAlignError("pseudobulk fragment source differs")
    allowed_folds = [
        value
        for value in range(5)
        if value not in {genomic_test_fold, genomic_valid_fold}
    ]
    output_directory.mkdir(mode=0o750)
    states = {}
    for fold in allowed_folds:
        path = output_directory / f"genomic_fold_{fold}.tagAlign.gz"
        text, hashed, raw = open_deterministic_gzip(path)
        states[fold] = {
            "path": path,
            "text": text,
            "hashed": hashed,
            "raw": raw,
            "unique_fragments": 0,
            "read_support": 0,
        }
    excluded_fragments = {"test": 0, "valid": 0}
    excluded_read_support = {"test": 0, "valid": 0}
    previous = None
    completed = False
    try:
        with gzip.open(source, "rt", encoding="utf-8", newline="") as handle:
            for line in handle:
                key, fields = parse_fragment_line(line)
                if previous is not None and key < previous:
                    raise TagAlignError("pseudobulk fragment source is not sorted")
                previous = key
                contig, start_text, end_text, cell_id, support_text = fields
                support = int(support_text)
                fold = fold_by_contig.get(contig)
                if fold in {genomic_test_fold, genomic_valid_fold}:
                    role = "test" if fold == genomic_test_fold else "valid"
                    excluded_fragments[role] += 1
                    excluded_read_support[role] += support
                    continue
                if fold not in states:
                    raise TagAlignError("fragment contig is outside the fold contract")
                start, end = int(start_text), int(end_text)
                text = states[fold]["text"]
                text.write(
                    "\t".join(
                        (contig, start_text, str(start + 1), f"{cell_id}/L", "1000", "+")
                    )
                    + "\n"
                )
                text.write(
                    "\t".join(
                        (contig, str(end - 1), end_text, f"{cell_id}/R", "1000", "-")
                    )
                    + "\n"
                )
                states[fold]["unique_fragments"] += 1
                states[fold]["read_support"] += support
        completed = True
    finally:
        for state in states.values():
            text = state["text"]
            raw = state["raw"]
            if completed:
                state["sha256"] = close_deterministic_gzip(
                    text, state["hashed"], raw
                )
            elif not text.closed:
                text.close()
                raw.close()
    outputs = []
    for fold in allowed_folds:
        state = states[fold]
        path = state["path"]
        outputs.append(
            {
                "genomic_fold": fold,
                "unique_fragments": state["unique_fragments"],
                "read_support": state["read_support"],
                "tn5_events": 2 * state["unique_fragments"],
                "path": path.name,
                "sha256": state["sha256"],
                "size_bytes": path.stat().st_size,
            }
        )
    return {
        "schema_version": "masld-bench-training-tagalign-by-genomic-fold-v1",
        "included_unique_fragments": sum(
            row["unique_fragments"] for row in outputs
        ),
        "included_read_support": sum(row["read_support"] for row in outputs),
        "included_tn5_events": sum(row["tn5_events"] for row in outputs),
        "excluded_genomic_test_fragments": excluded_fragments["test"],
        "excluded_genomic_test_read_support": excluded_read_support["test"],
        "excluded_genomic_valid_fragments": excluded_fragments["valid"],
        "excluded_genomic_valid_read_support": excluded_read_support["valid"],
        "genomic_test_fold": genomic_test_fold,
        "genomic_valid_fold": genomic_valid_fold,
        "allowed_genomic_folds": allowed_folds,
        "read_support_used_as_signal_weight": False,
        "left_insertion_position": "fragment_start",
        "right_insertion_position": "fragment_end_minus_1_BED_exclusive",
        "outputs": outputs,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pseudobulk", type=Path, required=True)
    parser.add_argument("--pseudobulk-artifacts-sha256", required=True)
    parser.add_argument("--split-contract", type=Path, required=True)
    parser.add_argument("--split-artifacts-sha256", required=True)
    parser.add_argument("--donor-test-fold", type=int, required=True)
    parser.add_argument("--genomic-test-fold", type=int, required=True)
    parser.add_argument("--genomic-valid-fold", type=int, required=True)
    parser.add_argument("--lineage", required=True)
    output_group = parser.add_mutually_exclusive_group(required=True)
    output_group.add_argument("--output", type=Path)
    output_group.add_argument("--output-directory", type=Path)
    parser.add_argument("--summary", type=Path, required=True)
    arguments = parser.parse_args()
    if (
        sha256_file(arguments.pseudobulk / "ARTIFACTS.json")
        != arguments.pseudobulk_artifacts_sha256
        or sha256_file(arguments.split_contract / "ARTIFACTS.json")
        != arguments.split_artifacts_sha256
    ):
        raise TagAlignError("input ARTIFACTS SHA-256 differs")
    rows = read_tsv(
        arguments.pseudobulk / "training_pseudobulk_manifest.tsv",
        PSEUDOBULK_FIELDS,
    )
    matches = [
        row
        for row in rows
        if int(row["donor_test_fold"]) == arguments.donor_test_fold
        and row["lineage_id"] == arguments.lineage
    ]
    if len(matches) != 1:
        raise TagAlignError("pseudobulk donor-fold/lineage row is not unique")
    row = matches[0]
    source = arguments.pseudobulk / row["fragment_path"]
    emit_arguments = {
        "source": source,
        "source_size_bytes": int(row["fragment_size_bytes"]),
        "source_sha256": row["fragment_sha256"],
        "fold_by_contig": genomic_fold_map(arguments.split_contract),
        "genomic_test_fold": arguments.genomic_test_fold,
        "genomic_valid_fold": arguments.genomic_valid_fold,
    }
    if arguments.output_directory is not None:
        result = emit_by_genomic_fold(
            **emit_arguments, output_directory=arguments.output_directory
        )
    else:
        result = emit(**emit_arguments, output=arguments.output)
    result.update(
        {
            "donor_test_fold": arguments.donor_test_fold,
            "donor_valid_fold": int(row["donor_valid_fold"]),
            "donor_train_folds": row["donor_train_folds"],
            "lineage_id": arguments.lineage,
            "source_fragment_sha256": row["fragment_sha256"],
        }
    )
    arguments.summary.write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

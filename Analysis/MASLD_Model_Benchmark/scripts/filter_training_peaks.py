#!/usr/bin/env python3
"""Filter MACS2 peaks without using held donor or genomic-test outcomes."""

from __future__ import annotations

import argparse
from bisect import bisect_left
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence


INPUT_LENGTH = 2114
OUTPUT_LENGTH = 1000
MAX_JITTER = 500
PEAKS_PER_ALLOWED_FOLD = 25_000
PRIMARY_CONTIGS = tuple([f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY"])


class TrainingPeakError(ValueError):
    """Raised when peak filtering does not meet the training-only requirements."""


def genomic_fold_map(split_contract: Path) -> dict[str, int]:
    path = split_contract / "genomic_folds.tsv"
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != (
            "genomic_fold",
            "contigs",
            "total_bp",
        ):
            raise TrainingPeakError("genomic-fold fields differ")
        rows = [dict(row) for row in reader]
    mapping = {
        contig: int(row["genomic_fold"])
        for row in rows
        for contig in row["contigs"].split(",")
    }
    if len(rows) != 5 or set(mapping.values()) != set(range(5)):
        raise TrainingPeakError("genomic-fold contract differs")
    return mapping


def merge_intervals(
    values: Iterable[tuple[int, int]],
) -> tuple[tuple[int, int], ...]:
    merged: list[list[int]] = []
    for start, end in sorted(values):
        if start < 0 or end <= start:
            raise TrainingPeakError("blacklist interval is invalid")
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        elif end > merged[-1][1]:
            merged[-1][1] = end
    return tuple((start, end) for start, end in merged)


def overlaps(
    intervals: Sequence[tuple[int, int]], start: int, end: int
) -> bool:
    starts = [value[0] for value in intervals]
    index = bisect_left(starts, end)
    return index > 0 and intervals[index - 1][1] > start


def peak_windows(
    start: int, summit_offset: int
) -> tuple[tuple[int, int], tuple[int, int]]:
    center = start + summit_offset
    sequence_half = INPUT_LENGTH // 2 + MAX_JITTER
    signal_half = OUTPUT_LENGTH // 2 + MAX_JITTER
    return (
        (center - sequence_half, center + sequence_half),
        (center - signal_half, center + signal_half),
    )


def selection_digest(fields: Sequence[str]) -> str:
    return sha256("\0".join(fields).encode("utf-8")).hexdigest()


def build(
    *,
    narrowpeak: Path,
    blacklist: Path,
    fasta: Path,
    fold_by_contig: Mapping[str, int],
    genomic_test_fold: int,
    genomic_valid_fold: int,
    output: Path,
    summary_path: Path,
    peaks_per_allowed_fold: int = PEAKS_PER_ALLOWED_FOLD,
    minimum_total_peaks: int = 1000,
) -> dict[str, object]:
    import pyfaidx

    if (
        output.exists()
        or summary_path.exists()
        or not 0 <= genomic_test_fold < 5
        or not 0 <= genomic_valid_fold < 5
        or genomic_test_fold == genomic_valid_fold
    ):
        raise TrainingPeakError("invalid output or genomic fold contract")
    blacklist_raw: dict[str, list[tuple[int, int]]] = {}
    with gzip.open(blacklist, "rt", encoding="utf-8", newline="") as handle:
        for line in handle:
            contig, start_text, end_text = line.rstrip("\n").split("\t")
            blacklist_raw.setdefault(contig, []).append(
                (int(start_text), int(end_text))
            )
    blacklist_by_contig = {
        contig: merge_intervals(values) for contig, values in blacklist_raw.items()
    }
    reference = pyfaidx.Fasta(
        str(fasta),
        as_raw=True,
        sequence_always_upper=True,
        rebuild=False,
        read_ahead=1_000_000,
    )
    counts = {
        "input": 0,
        "held_genomic_test": 0,
        "held_genomic_valid": 0,
        "blacklist_signal_window": 0,
        "edge_or_ambiguous_sequence": 0,
    }
    eligible_by_fold: dict[int, list[tuple[tuple[object, ...], list[str]]]] = {
        fold: []
        for fold in range(5)
        if fold not in {genomic_test_fold, genomic_valid_fold}
    }
    try:
        with narrowpeak.open("r", encoding="utf-8", newline="") as handle:
            for line_number, line in enumerate(handle, start=1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 10:
                    raise TrainingPeakError(
                        f"narrowPeak row {line_number} does not have ten fields"
                    )
                counts["input"] += 1
                contig = fields[0]
                fold = fold_by_contig.get(contig)
                if fold is None:
                    raise TrainingPeakError("MACS2 peak contig is outside the contract")
                if fold == genomic_test_fold:
                    counts["held_genomic_test"] += 1
                    continue
                if fold == genomic_valid_fold:
                    counts["held_genomic_valid"] += 1
                    continue
                start, end, summit = int(fields[1]), int(fields[2]), int(fields[9])
                if start < 0 or end <= start or summit < 0:
                    raise TrainingPeakError("MACS2 peak geometry differs")
                sequence_window, signal_window = peak_windows(start, summit)
                blacklist_intervals = blacklist_by_contig.get(contig, ())
                if overlaps(blacklist_intervals, *signal_window):
                    counts["blacklist_signal_window"] += 1
                    continue
                sequence_start, sequence_end = sequence_window
                if sequence_start < 0 or sequence_end > len(reference[contig]):
                    counts["edge_or_ambiguous_sequence"] += 1
                    continue
                sequence = str(reference[contig][sequence_start:sequence_end]).upper()
                if len(sequence) != INPUT_LENGTH + 2 * MAX_JITTER or set(
                    sequence
                ).difference("ACGT"):
                    counts["edge_or_ambiguous_sequence"] += 1
                    continue
                rank = (
                    -float(fields[8]),
                    -float(fields[7]),
                    -float(fields[6]),
                    selection_digest(fields),
                )
                eligible_by_fold[fold].append((rank, fields))
    finally:
        reference.close()

    selected = []
    fold_census = {}
    for fold, values in sorted(eligible_by_fold.items()):
        values.sort(key=lambda row: row[0])
        retained = values[:peaks_per_allowed_fold]
        selected.extend(row[1] for row in retained)
        fold_census[str(fold)] = {
            "eligible": len(values),
            "retained": len(retained),
        }
    if len(selected) < minimum_total_peaks:
        raise TrainingPeakError("too few training peaks remain")
    selected.sort(
        key=lambda row: (
            PRIMARY_CONTIGS.index(row[0]),
            int(row[1]),
            int(row[2]),
        )
    )
    with output.open("x", encoding="utf-8", newline="") as handle:
        for fields in selected:
            handle.write("\t".join(fields) + "\n")
    summary = {
        **counts,
        "retained": len(selected),
        "genomic_test_fold": genomic_test_fold,
        "genomic_valid_fold": genomic_valid_fold,
        "allowed_genomic_folds": sorted(eligible_by_fold),
        "peaks_per_allowed_fold_cap": peaks_per_allowed_fold,
        "fold_census": fold_census,
        "held_atac_test_used_for_peak_selection": False,
        "held_atac_valid_used_for_peak_selection": False,
        "blacklist_policy": "exclude_max_jittered_output_window_overlap",
        "sequence_policy": "unambiguous_max_jittered_2114bp_input",
    }
    summary_path.write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--narrowpeak", type=Path, required=True)
    parser.add_argument("--blacklist", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--split-contract", type=Path, required=True)
    parser.add_argument("--genomic-test-fold", type=int, required=True)
    parser.add_argument("--genomic-valid-fold", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument(
        "--peaks-per-allowed-fold", type=int, default=PEAKS_PER_ALLOWED_FOLD
    )
    arguments = parser.parse_args()
    result = build(
        narrowpeak=arguments.narrowpeak,
        blacklist=arguments.blacklist,
        fasta=arguments.fasta,
        fold_by_contig=genomic_fold_map(arguments.split_contract),
        genomic_test_fold=arguments.genomic_test_fold,
        genomic_valid_fold=arguments.genomic_valid_fold,
        output=arguments.output,
        summary_path=arguments.summary,
        peaks_per_allowed_fold=arguments.peaks_per_allowed_fold,
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

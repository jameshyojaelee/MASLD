#!/usr/bin/env python3
"""Prepare train-only sequence-control coordinates and loss weighting."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import random
from typing import Mapping, Sequence


PRIMARY_CONTIGS = tuple([f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY"])


class SequenceControlPreparationError(ValueError):
    """Raised when sequence-control preparation violates the split contract."""


def read_fold(path: Path) -> dict[str, list[str]]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if set(raw) != {"train", "valid", "test"} or any(
        not isinstance(raw[role], list) for role in raw
    ):
        raise SequenceControlPreparationError("fold JSON fields differ")
    fold = {role: [str(value) for value in raw[role]] for role in raw}
    flattened = [
        contig for role in ("train", "valid", "test") for contig in fold[role]
    ]
    if (
        not all(fold.values())
        or len(flattened) != len(set(flattened))
        or set(flattened) != set(PRIMARY_CONTIGS)
    ):
        raise SequenceControlPreparationError("fold JSON is not a primary partition")
    return fold


def _read_narrowpeak(path: Path) -> list[tuple[str, ...]]:
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = tuple(line.rstrip("\n").split("\t"))
            if len(fields) != 10:
                raise SequenceControlPreparationError(
                    f"{path.name} row {line_number} differs"
                )
            try:
                start = int(fields[1])
                end = int(fields[2])
                summit = int(fields[9])
            except ValueError as exc:
                raise SequenceControlPreparationError(
                    f"{path.name} row {line_number} coordinates differ"
                ) from exc
            if start < 0 or end <= start or summit < 0:
                raise SequenceControlPreparationError(
                    f"{path.name} row {line_number} coordinates differ"
                )
            rows.append(fields)
    if not rows:
        raise SequenceControlPreparationError(f"{path.name} is empty")
    return rows


def _center(row: Sequence[str]) -> int:
    return int(row[1]) + int(row[9])


def _edge_filter(
    rows: Sequence[tuple[str, ...]], chrom_sizes: Mapping[str, int], width: int
) -> list[tuple[str, ...]]:
    kept = []
    half = width // 2
    for row in rows:
        size = chrom_sizes.get(row[0])
        center = _center(row)
        if size is not None and center - half >= 0 and center + half <= size:
            kept.append(row)
    return kept


def _counts(
    rows: Sequence[tuple[str, ...]], bigwig, output_length: int
):
    import numpy as np

    values = []
    half = output_length // 2
    for row in rows:
        center = _center(row)
        signal = bigwig.values(row[0], center - half, center + half, numpy=True)
        if signal is None or len(signal) != output_length:
            raise SequenceControlPreparationError("bigWig output window differs")
        values.append(float(np.nan_to_num(signal).sum()))
    result = np.asarray(values, dtype=np.float64)
    if result.size == 0 or not np.isfinite(result).all() or np.any(result < 0):
        raise SequenceControlPreparationError("training counts are invalid")
    return result


def _role_census(
    rows: Sequence[tuple[str, ...]], fold: Mapping[str, Sequence[str]]
) -> dict[str, int]:
    role_by_contig = {
        contig: role for role, contigs in fold.items() for contig in contigs
    }
    counts = Counter(role_by_contig.get(row[0], "outside") for row in rows)
    if counts.get("outside", 0):
        raise SequenceControlPreparationError("region is outside the fold contract")
    return {role: counts.get(role, 0) for role in ("train", "valid", "test")}


def _write_rows(path: Path, rows: Sequence[Sequence[str]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        for row in rows:
            handle.write("\t".join(row) + "\n")


def prepare(
    *,
    bigwig_path: Path,
    peaks_path: Path,
    nonpeaks_path: Path,
    fold_path: Path,
    output: Path,
    seed: int,
    negative_sampling_ratio: float,
    outlier_threshold: float,
    input_length: int,
    output_length: int,
    max_jitter: int,
) -> dict[str, object]:
    import numpy as np
    import pyBigWig

    if (
        output.exists()
        or seed < 0
        or not 0.0 <= negative_sampling_ratio <= 1.0
        or not 0.5 < outlier_threshold < 1.0
        or input_length != 2114
        or output_length != 1000
        or max_jitter < 0
    ):
        raise SequenceControlPreparationError("preparation arguments differ")
    for path in (bigwig_path, peaks_path, nonpeaks_path, fold_path):
        if path.is_symlink() or not path.is_file():
            raise SequenceControlPreparationError("preparation input is not regular")
    random.seed(seed)
    np.random.seed(seed)
    fold = read_fold(fold_path)
    peaks = _read_narrowpeak(peaks_path)
    nonpeaks = _read_narrowpeak(nonpeaks_path)
    role_by_contig = {
        contig: role for role, contigs in fold.items() for contig in contigs
    }
    if any(row[0] not in role_by_contig for row in peaks + nonpeaks):
        raise SequenceControlPreparationError("region is outside the fold contract")
    train_peaks = [row for row in peaks if role_by_contig[row[0]] == "train"]
    held_peaks = [row for row in peaks if role_by_contig[row[0]] != "train"]
    train_nonpeaks = [row for row in nonpeaks if role_by_contig[row[0]] == "train"]
    held_nonpeaks = [row for row in nonpeaks if role_by_contig[row[0]] != "train"]

    bigwig = pyBigWig.open(str(bigwig_path))
    try:
        chrom_sizes = {str(key): int(value) for key, value in bigwig.chroms().items()}
        train_peaks = _edge_filter(
            train_peaks, chrom_sizes, input_length + 2 * max_jitter
        )
        held_peaks = _edge_filter(held_peaks, chrom_sizes, input_length)
        train_nonpeaks = _edge_filter(train_nonpeaks, chrom_sizes, input_length)
        held_nonpeaks = _edge_filter(held_nonpeaks, chrom_sizes, input_length)
        peak_counts = _counts(train_peaks, bigwig, output_length)
        nonpeak_counts = _counts(train_nonpeaks, bigwig, output_length)
    finally:
        bigwig.close()

    negative_count = int(negative_sampling_ratio * len(peak_counts))
    if negative_count > len(nonpeak_counts):
        raise SequenceControlPreparationError("not enough training nonpeaks")
    if negative_count:
        sampled = np.random.choice(
            nonpeak_counts, replace=False, size=negative_count
        )
        combined = np.concatenate((peak_counts, sampled))
    else:
        combined = peak_counts
    upper = float(np.quantile(combined, outlier_threshold))
    lower = float(np.quantile(combined, 1.0 - outlier_threshold))
    peak_keep = (peak_counts < upper) & (peak_counts > lower)
    nonpeak_keep = (nonpeak_counts < upper) & (nonpeak_counts > lower)
    filtered_train_peaks = [
        row for row, keep in zip(train_peaks, peak_keep) if bool(keep)
    ]
    filtered_train_nonpeaks = [
        row for row, keep in zip(train_nonpeaks, nonpeak_keep) if bool(keep)
    ]
    if len(filtered_train_nonpeaks) > len(filtered_train_peaks):
        indices = np.random.RandomState(1).choice(
            len(filtered_train_nonpeaks),
            size=len(filtered_train_peaks),
            replace=False,
        )
        filtered_train_nonpeaks = [
            filtered_train_nonpeaks[int(index)] for index in indices
        ]
    retained_counts = combined[(combined <= upper) & (combined >= lower)]
    if retained_counts.size == 0:
        raise SequenceControlPreparationError("no counts survive train-only filtering")
    counts_loss_weight = round(
        max(float(np.median(retained_counts) / 10.0), 1.0), 2
    )
    if not filtered_train_peaks or not filtered_train_nonpeaks:
        raise SequenceControlPreparationError("train-only filtering emptied a class")

    output.mkdir(mode=0o750)
    artifact_fold = output / "actual_fold.json"
    artifact_fold.write_text(
        json.dumps(fold, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    filtered_peaks = filtered_train_peaks + held_peaks
    filtered_nonpeaks = filtered_train_nonpeaks + held_nonpeaks
    _write_rows(output / "sequence_control.filtered.peaks.bed", filtered_peaks)
    _write_rows(output / "sequence_control.filtered.nonpeaks.bed", filtered_nonpeaks)
    parameters = {
        "counts_loss_weight": str(counts_loss_weight),
        "filters": "1",
        "n_dil_layers": "1",
        "inputlen": str(input_length),
        "outputlen": str(output_length),
        "max_jitter": str(max_jitter),
        "chr_fold_path": artifact_fold.name,
        "negative_sampling_ratio": str(negative_sampling_ratio),
    }
    parameter_path = output / "sequence_control_model_params.actual_fold.tsv"
    with parameter_path.open("x", encoding="utf-8") as handle:
        for key, value in parameters.items():
            handle.write(f"{key}\t{value}\n")
    summary = {
        "schema_version": "masld-bench-sequence-control-training-preparation-v1",
        "status": "pass",
        "seed": seed,
        "input_length": input_length,
        "output_length": output_length,
        "max_jitter": max_jitter,
        "negative_sampling_ratio": negative_sampling_ratio,
        "outlier_threshold": outlier_threshold,
        "lower_count_threshold": lower,
        "upper_count_threshold": upper,
        "counts_loss_weight": counts_loss_weight,
        "input_peak_role_census": _role_census(peaks, fold),
        "input_nonpeak_role_census": _role_census(nonpeaks, fold),
        "filtered_peak_role_census": _role_census(filtered_peaks, fold),
        "filtered_nonpeak_role_census": _role_census(filtered_nonpeaks, fold),
        "hyperparameter_fit_outcome_roles": ["train"],
        "validation_outcomes_used_for_hyperparameter_fit": False,
        "test_outcomes_used_for_hyperparameter_fit": False,
        "held_peak_rows_preserved_after_outcome_filtering": True,
        "held_nonpeak_rows_preserved_after_outcome_filtering": True,
        "bias_model_loaded": False,
        "upstream_equation": (
            "ChromBPNet_1.0.1_median_train_profile_count_divided_by_10_"
            "with_minimum_1"
        ),
        "final_model_parameters": parameters,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bigwig", type=Path, required=True)
    parser.add_argument("--peaks", type=Path, required=True)
    parser.add_argument("--nonpeaks", type=Path, required=True)
    parser.add_argument("--fold", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260824)
    parser.add_argument("--negative-sampling-ratio", type=float, default=0.1)
    parser.add_argument("--outlier-threshold", type=float, default=0.9999)
    parser.add_argument("--input-length", type=int, default=2114)
    parser.add_argument("--output-length", type=int, default=1000)
    parser.add_argument("--max-jitter", type=int, default=500)
    arguments = parser.parse_args()
    result = prepare(
        bigwig_path=arguments.bigwig,
        peaks_path=arguments.peaks,
        nonpeaks_path=arguments.nonpeaks,
        fold_path=arguments.fold,
        output=arguments.output,
        seed=arguments.seed,
        negative_sampling_ratio=arguments.negative_sampling_ratio,
        outlier_threshold=arguments.outlier_threshold,
        input_length=arguments.input_length,
        output_length=arguments.output_length,
        max_jitter=arguments.max_jitter,
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

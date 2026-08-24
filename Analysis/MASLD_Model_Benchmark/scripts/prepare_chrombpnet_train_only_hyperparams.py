#!/usr/bin/env python3
"""Fit ChromBPNet hyperparameters without genomic-validation outcomes."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
import random
from types import SimpleNamespace
from typing import Mapping, Sequence


PRIMARY_CONTIGS = tuple([f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY"])
PARAMETER_KEYS = (
    "counts_loss_weight",
    "filters",
    "n_dil_layers",
    "bias_model_path",
    "inputlen",
    "outputlen",
    "max_jitter",
    "chr_fold_path",
    "negative_sampling_ratio",
)
DATA_PARAMETER_KEYS = (
    "counts_sum_min_thresh",
    "counts_sum_max_thresh",
    "trainings_pts_post_thresh",
)


class TrainOnlyHyperparameterError(ValueError):
    """Raised when train-only ChromBPNet hyperparameter fitting leaks outcomes."""


def read_fold(path: Path) -> dict[str, list[str]]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if set(raw) != {"train", "valid", "test"} or any(
        not isinstance(raw[role], list) for role in raw
    ):
        raise TrainOnlyHyperparameterError("actual fold JSON fields differ")
    fold = {role: [str(value) for value in raw[role]] for role in raw}
    flattened = [contig for role in ("train", "valid", "test") for contig in fold[role]]
    if (
        len(fold["train"]) == 0
        or len(fold["valid"]) == 0
        or len(fold["test"]) == 0
        or len(flattened) != len(set(flattened))
        or set(flattened) != set(PRIMARY_CONTIGS)
    ):
        raise TrainOnlyHyperparameterError("actual fold JSON is not a primary partition")
    return fold


def build_hyperparameter_fold(actual: Mapping[str, Sequence[str]]) -> dict[str, list[str]]:
    return {
        "train": list(actual["train"]),
        "valid": [],
        "test": list(actual["valid"]) + list(actual["test"]),
    }


def read_two_column_rows(path: Path) -> list[tuple[str, str]]:
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise TrainOnlyHyperparameterError(
                    f"parameter row {line_number} differs"
                )
            rows.append((fields[0], fields[1]))
    return rows


def read_parameter_rows(path: Path) -> list[tuple[str, str]]:
    rows = read_two_column_rows(path)
    if tuple(key for key, _value in rows) != PARAMETER_KEYS:
        raise TrainOnlyHyperparameterError("ChromBPNet model parameter keys differ")
    return rows


def write_actual_fold_parameters(
    source: Path,
    actual_fold_binding: str,
    bias_model_binding: str,
    output: Path,
) -> dict[str, str]:
    if (
        output.exists()
        or Path(actual_fold_binding).name != actual_fold_binding
        or Path(bias_model_binding).name != bias_model_binding
    ):
        raise TrainOnlyHyperparameterError("invalid artifact-relative parameter binding")
    rows = read_parameter_rows(source)
    rewritten = []
    for key, value in rows:
        if key == "chr_fold_path":
            value = actual_fold_binding
        elif key == "bias_model_path":
            value = bias_model_binding
        rewritten.append((key, value))
    with output.open("x", encoding="utf-8", newline="") as handle:
        for key, value in rewritten:
            handle.write(f"{key}\t{value}\n")
    return dict(rewritten)


def read_narrowpeak(path: Path) -> list[tuple[str, ...]]:
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = tuple(line.rstrip("\n").split("\t"))
            if len(fields) != 10:
                raise TrainOnlyHyperparameterError(
                    f"{path.name} row {line_number} differs"
                )
            rows.append(fields)
    return rows


def role_census(
    rows: Sequence[tuple[str, ...]], actual: Mapping[str, Sequence[str]]
) -> dict[str, int]:
    role_by_contig = {
        contig: role for role, contigs in actual.items() for contig in contigs
    }
    counts = Counter(role_by_contig.get(row[0], "outside") for row in rows)
    if counts.get("outside", 0):
        raise TrainOnlyHyperparameterError("narrowPeak row is outside the fold contract")
    return {role: counts.get(role, 0) for role in ("train", "valid", "test")}


def held_rows(
    rows: Sequence[tuple[str, ...]], actual: Mapping[str, Sequence[str]]
) -> Counter[tuple[object, ...]]:
    held = set(actual["valid"]) | set(actual["test"])
    return Counter(narrowpeak_semantic_row(row) for row in rows if row[0] in held)


def narrowpeak_semantic_row(row: tuple[str, ...]) -> tuple[object, ...]:
    """Normalize narrowPeak numeric payloads while preserving row identity."""

    if len(row) != 10:
        raise TrainOnlyHyperparameterError("narrowPeak semantic row differs")
    try:
        summit = int(row[9])
    except ValueError as error:
        raise TrainOnlyHyperparameterError(
            "narrowPeak semantic numeric field differs"
        ) from error
    numeric_payload: list[float | None] = []
    for value in row[6:9]:
        if value == ".":
            numeric_payload.append(None)
            continue
        try:
            numeric = float(value)
        except ValueError as error:
            raise TrainOnlyHyperparameterError(
                "narrowPeak semantic numeric field differs"
            ) from error
        if not math.isfinite(numeric):
            raise TrainOnlyHyperparameterError(
                "narrowPeak semantic numeric field is non-finite"
            )
        numeric_payload.append(numeric)
    return (*row[:6], *numeric_payload, summit)


def prepare(
    *,
    genome: Path,
    bigwig: Path,
    peaks: Path,
    nonpeaks: Path,
    actual_fold: Path,
    bias_model: Path,
    output: Path,
    seed: int,
    negative_sampling_ratio: float,
    outlier_threshold: float,
    input_length: int,
    output_length: int,
    max_jitter: int,
    filters: int,
    dilation_layers: int,
) -> dict[str, object]:
    if output.exists() or not 0 <= negative_sampling_ratio <= 1:
        raise TrainOnlyHyperparameterError("invalid output or sampling ratio")
    actual = read_fold(actual_fold.resolve(strict=True))
    output.mkdir(mode=0o750)
    artifact_actual_fold = output / "actual_fold.json"
    artifact_actual_fold.write_text(
        json.dumps(actual, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    hyperparameter_fold = output / "hyperparameter_fold.json"
    hyperparameter_fold.write_text(
        json.dumps(build_hyperparameter_fold(actual), sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )
    prefix = output / "train_only_"

    import numpy as np
    import tensorflow as tf
    from chrombpnet.helpers.hyperparameters import find_chrombpnet_hyperparams

    random.seed(seed)
    np.random.seed(seed)
    tf.random.set_seed(seed)
    arguments = SimpleNamespace(
        genome=str(genome.resolve(strict=True)),
        bigwig=str(bigwig.resolve(strict=True)),
        peaks=str(peaks.resolve(strict=True)),
        nonpeaks=str(nonpeaks.resolve(strict=True)),
        negative_sampling_ratio=negative_sampling_ratio,
        outlier_threshold=outlier_threshold,
        max_jitter=max_jitter,
        chr_fold_path=str(hyperparameter_fold.resolve(strict=True)),
        inputlen=input_length,
        outputlen=output_length,
        filters=filters,
        n_dilation_layers=dilation_layers,
        bias_model_path=str(bias_model.resolve(strict=True)),
        output_prefix=str(prefix),
    )
    find_chrombpnet_hyperparams.main(arguments)

    filtered_peaks = output / "train_only_filtered.peaks.bed"
    filtered_nonpeaks = output / "train_only_filtered.nonpeaks.bed"
    upstream_parameters = output / "train_only_chrombpnet_model_params.tsv"
    actual_parameters = output / "chrombpnet_model_params.actual_fold.tsv"
    final_parameters = write_actual_fold_parameters(
        upstream_parameters,
        artifact_actual_fold.name,
        "train_only_bias_model_scaled.h5",
        actual_parameters,
    )
    input_peak_rows = read_narrowpeak(peaks)
    input_nonpeak_rows = read_narrowpeak(nonpeaks)
    filtered_peak_rows = read_narrowpeak(filtered_peaks)
    filtered_nonpeak_rows = read_narrowpeak(filtered_nonpeaks)
    if held_rows(input_peak_rows, actual) != held_rows(filtered_peak_rows, actual):
        raise TrainOnlyHyperparameterError("held peak rows changed during fitting")
    if held_rows(input_nonpeak_rows, actual) != held_rows(
        filtered_nonpeak_rows, actual
    ):
        raise TrainOnlyHyperparameterError("held nonpeak rows changed during fitting")
    data_parameter_rows = read_two_column_rows(
        output / "train_only_chrombpnet_data_params.tsv"
    )
    if tuple(key for key, _value in data_parameter_rows) != DATA_PARAMETER_KEYS:
        raise TrainOnlyHyperparameterError("ChromBPNet data parameter keys differ")
    data_parameters = dict(data_parameter_rows)
    summary = {
        "schema_version": "masld-bench-chrombpnet-train-only-hyperparams-v1",
        "status": "pass",
        "seed": seed,
        "negative_sampling_ratio": negative_sampling_ratio,
        "outlier_threshold": outlier_threshold,
        "input_length": input_length,
        "output_length": output_length,
        "max_jitter": max_jitter,
        "filters": filters,
        "dilation_layers": dilation_layers,
        "actual_fold": actual,
        "hyperparameter_fold": build_hyperparameter_fold(actual),
        "input_peak_role_census": role_census(input_peak_rows, actual),
        "filtered_peak_role_census": role_census(filtered_peak_rows, actual),
        "input_nonpeak_role_census": role_census(input_nonpeak_rows, actual),
        "filtered_nonpeak_role_census": role_census(filtered_nonpeak_rows, actual),
        "data_parameters": data_parameters,
        "final_model_parameters": final_parameters,
        "hyperparameter_fit_outcome_roles": ["train"],
        "validation_outcomes_used_for_hyperparameter_fit": False,
        "test_outcomes_used_for_hyperparameter_fit": False,
        "held_coordinate_rows_preserved": True,
        "final_parameter_paths_are_artifact_relative": True,
        "upstream_helper": "chrombpnet_1.0.1_find_chrombpnet_hyperparams",
    }
    (output / "summary.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--bigwig", type=Path, required=True)
    parser.add_argument("--peaks", type=Path, required=True)
    parser.add_argument("--nonpeaks", type=Path, required=True)
    parser.add_argument("--actual-fold", type=Path, required=True)
    parser.add_argument("--bias-model", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260824)
    parser.add_argument("--negative-sampling-ratio", type=float, default=0.1)
    parser.add_argument("--outlier-threshold", type=float, default=0.9999)
    parser.add_argument("--input-length", type=int, default=2114)
    parser.add_argument("--output-length", type=int, default=1000)
    parser.add_argument("--max-jitter", type=int, default=500)
    parser.add_argument("--filters", type=int, default=512)
    parser.add_argument("--dilation-layers", type=int, default=8)
    arguments = parser.parse_args()
    result = prepare(
        genome=arguments.genome,
        bigwig=arguments.bigwig,
        peaks=arguments.peaks,
        nonpeaks=arguments.nonpeaks,
        actual_fold=arguments.actual_fold,
        bias_model=arguments.bias_model,
        output=arguments.output,
        seed=arguments.seed,
        negative_sampling_ratio=arguments.negative_sampling_ratio,
        outlier_threshold=arguments.outlier_threshold,
        input_length=arguments.input_length,
        output_length=arguments.output_length,
        max_jitter=arguments.max_jitter,
        filters=arguments.filters,
        dilation_layers=arguments.dilation_layers,
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

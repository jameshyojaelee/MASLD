#!/usr/bin/env python3
"""Predict fixed held cCRE windows with sequence-only BPNet."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import random
from typing import Any, Iterable, Sequence

import numpy as np

try:
    from .bpnet_contract import (
        BPNetContractError,
        DEFAULT_SEED,
        INPUT_LENGTH,
        OUTPUT_LENGTH,
        one_hot_dna,
        read_fold,
        reverse_complement_one_hot,
        softmax,
        validate_model_inputs,
    )
except ImportError:
    from bpnet_contract import (
        BPNetContractError,
        DEFAULT_SEED,
        INPUT_LENGTH,
        OUTPUT_LENGTH,
        one_hot_dna,
        read_fold,
        reverse_complement_one_hot,
        softmax,
        validate_model_inputs,
    )


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
OUTPUT_FIELDS = (
    "window_id",
    "selection_hash",
    "contig",
    "output_start",
    "output_end",
    "genomic_fold",
    "role",
    "ccre_class",
    "forward_log1p_count",
    "reverse_log1p_count",
    "strand_averaged_mass",
)


def read_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise BPNetContractError(f"{path.name} fields differ")
        return [dict(row) for row in reader]


def batches(values: Sequence[Any], batch_size: int) -> Iterable[Sequence[Any]]:
    for start in range(0, len(values), batch_size):
        yield values[start : start + batch_size]


def count_mass(log1p_counts: np.ndarray) -> np.ndarray:
    values = np.asarray(log1p_counts, dtype=np.float64).reshape(-1)
    if not np.isfinite(values).all() or np.any(values > 50):
        raise BPNetContractError("BPNet log1p-count output is invalid")
    return np.maximum(np.expm1(values), 0.0)


def predict(
    *,
    model_path: Path,
    fasta_path: Path,
    windows_path: Path,
    fold_path: Path,
    output: Path,
    batch_size: int,
    seed: int,
    expected_windows: int = 32_000,
) -> dict[str, object]:
    if output.exists() or batch_size < 1 or seed < 0 or expected_windows < 1:
        raise BPNetContractError("invalid output, batch size, seed, or window count")
    validate_model_inputs((model_path, fasta_path, windows_path, fold_path))

    import h5py
    import pyfaidx
    import tensorflow as tf

    random.seed(seed)
    np.random.seed(seed)
    tf.random.set_seed(seed)
    fold = read_fold(fold_path)
    role_by_contig = {
        contig: role for role, contigs in fold.items() for contig in contigs
    }
    windows = [
        row
        for row in read_tsv(windows_path, CCRE_FIELDS)
        if role_by_contig.get(row["contig"]) in {"valid", "test"}
    ]
    if len(windows) != expected_windows or len(
        {row["window_id"] for row in windows}
    ) != len(windows):
        raise BPNetContractError("fixed held-window census differs")

    model = tf.keras.models.load_model(model_path, compile=False)
    if model.input_shape != (None, INPUT_LENGTH, 4):
        raise BPNetContractError("BPNet input geometry differs")
    if [tuple(value) for value in model.output_shape] != [
        (None, OUTPUT_LENGTH, 1),
        (None, 1),
    ]:
        raise BPNetContractError("BPNet output geometry differs")
    if len(model.inputs) != 1:
        raise BPNetContractError("BPNet inference model accepts bias/control input")

    output.mkdir(mode=0o750)
    reference = pyfaidx.Fasta(
        str(fasta_path),
        as_raw=True,
        sequence_always_upper=True,
        rebuild=False,
        read_ahead=1_000_000,
    )
    profile_path = output / "profile_probabilities.h5"
    counts_path = output / "regional_counts.tsv"
    string_type = h5py.string_dtype(encoding="utf-8")
    role_counts = {"valid": 0, "test": 0}
    count_min = float("inf")
    count_max = float("-inf")
    try:
        with h5py.File(profile_path, "x") as profiles, counts_path.open(
            "x", encoding="utf-8", newline=""
        ) as counts_handle:
            profile_values = profiles.create_dataset(
                "profile_probability",
                shape=(len(windows), OUTPUT_LENGTH),
                dtype="f4",
                chunks=(min(batch_size, len(windows)), OUTPUT_LENGTH),
                compression="gzip",
                compression_opts=1,
            )
            profiles.create_dataset(
                "window_id",
                data=np.asarray([row["window_id"] for row in windows], dtype=object),
                dtype=string_type,
            )
            profiles.create_dataset(
                "selection_hash",
                data=np.asarray(
                    [row["selection_hash"] for row in windows], dtype=object
                ),
                dtype=string_type,
            )
            profiles.attrs["schema_version"] = (
                "masld-bench-bpnet-ccre-profile-probability-v1"
            )
            profiles.attrs["strand_averaging"] = (
                "mean_probability_after_reverse_axis_restoration"
            )
            writer = csv.DictWriter(
                counts_handle, fieldnames=OUTPUT_FIELDS, delimiter="\t"
            )
            writer.writeheader()
            offset = 0
            for batch in batches(windows, batch_size):
                sequences = []
                for row in batch:
                    start, end = int(row["input_start"]), int(row["input_end"])
                    if start < 0 or end - start != INPUT_LENGTH:
                        raise BPNetContractError("input window geometry differs")
                    sequence = str(reference[row["contig"]][start:end]).upper()
                    if len(sequence) != INPUT_LENGTH:
                        raise BPNetContractError("reference window length differs")
                    sequences.append(sequence)
                forward = one_hot_dna(sequences)
                reverse = reverse_complement_one_hot(forward)
                forward_profile_raw, forward_count_raw = model.predict(
                    forward, batch_size=batch_size, verbose=0
                )
                reverse_profile_raw, reverse_count_raw = model.predict(
                    reverse, batch_size=batch_size, verbose=0
                )
                forward_probability = softmax(forward_profile_raw)
                reverse_probability = softmax(reverse_profile_raw)[:, ::-1]
                averaged_probability = (
                    0.5 * (forward_probability + reverse_probability)
                ).astype(np.float32)
                if not np.allclose(
                    averaged_probability.sum(axis=1), 1.0, rtol=0, atol=1e-6
                ):
                    raise BPNetContractError(
                        "strand-averaged profile does not sum to one"
                    )
                forward_count = np.asarray(forward_count_raw).reshape(-1)
                reverse_count = np.asarray(reverse_count_raw).reshape(-1)
                averaged_mass = 0.5 * (
                    count_mass(forward_count) + count_mass(reverse_count)
                )
                profile_values[offset : offset + len(batch)] = averaged_probability
                for index, row in enumerate(batch):
                    role = role_by_contig[row["contig"]]
                    role_counts[role] += 1
                    value = float(averaged_mass[index])
                    count_min = min(count_min, value)
                    count_max = max(count_max, value)
                    writer.writerow(
                        {
                            "window_id": row["window_id"],
                            "selection_hash": row["selection_hash"],
                            "contig": row["contig"],
                            "output_start": row["output_start"],
                            "output_end": row["output_end"],
                            "genomic_fold": row["genomic_fold"],
                            "role": role,
                            "ccre_class": row["ccre_class"],
                            "forward_log1p_count": format(
                                float(forward_count[index]), ".17g"
                            ),
                            "reverse_log1p_count": format(
                                float(reverse_count[index]), ".17g"
                            ),
                            "strand_averaged_mass": format(value, ".17g"),
                        }
                    )
                offset += len(batch)
    finally:
        reference.close()

    expected_each = expected_windows // 2
    if expected_windows % 2 or role_counts != {
        "valid": expected_each,
        "test": expected_each,
    }:
        raise BPNetContractError("held prediction role census differs")
    summary = {
        "schema_version": "masld-bench-bpnet-ccre-predictions-v1",
        "status": "pass",
        "model_id": "bpnet",
        "source_revision": "f4593eedca51741f25b8c5cc0c0d647faf8c8a3c",
        "windows": len(windows),
        "role_counts": role_counts,
        "input_length": INPUT_LENGTH,
        "output_length": OUTPUT_LENGTH,
        "batch_size": batch_size,
        "seed": seed,
        "reverse_complement_tta": True,
        "profile_strand_averaging": (
            "mean_probability_after_reverse_axis_restoration"
        ),
        "regional_mass_semantics": "arithmetic_mean_of_expm1_log1p_count",
        "minimum_predicted_count": count_min,
        "maximum_predicted_count": count_max,
        "observed_atac_input_exposed": False,
        "benchmark_metrics_calculated": False,
        "regional_counts_path": counts_path.name,
        "profile_probabilities_path": profile_path.name,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--windows", type=Path, required=True)
    parser.add_argument("--fold", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--expected-windows", type=int, default=32_000)
    arguments = parser.parse_args()
    try:
        result = predict(
            model_path=arguments.model,
            fasta_path=arguments.fasta,
            windows_path=arguments.windows,
            fold_path=arguments.fold,
            output=arguments.output,
            batch_size=arguments.batch_size,
            seed=arguments.seed,
            expected_windows=arguments.expected_windows,
        )
    except BPNetContractError as error:
        raise SystemExit(str(error)) from error
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

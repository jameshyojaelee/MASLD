#!/usr/bin/env python3
"""Predict fixed cCRE profiles with ChromBPNet without exposing ATAC outcomes."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import random
from typing import Any, Iterable, Sequence

import numpy as np


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
    "forward_logcount",
    "reverse_logcount",
    "strand_averaged_mass",
)
PRIMARY_CONTIGS = tuple([f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY"])


class ChromBPNetPredictionError(ValueError):
    """Raised when fixed-window ChromBPNet inference does not meet its requirements."""


def one_hot_dna(sequences: Sequence[str]) -> np.ndarray:
    if not sequences or len({len(sequence) for sequence in sequences}) != 1:
        raise ChromBPNetPredictionError("DNA sequence batch has invalid geometry")
    encoded = np.frombuffer("".join(sequences).encode("ascii"), dtype=np.uint8)
    lookup = np.full(256, -1, dtype=np.int8)
    for index, base in enumerate(b"ACGT"):
        lookup[base] = index
    indices = lookup[encoded]
    if np.any(indices < 0):
        raise ChromBPNetPredictionError("DNA sequence batch contains an ambiguous base")
    shape = (len(sequences), len(sequences[0]))
    return np.eye(4, dtype=np.float32)[indices.reshape(shape)]


def reverse_complement_one_hot(values: np.ndarray) -> np.ndarray:
    matrix = np.asarray(values)
    if matrix.ndim != 3 or matrix.shape[2] != 4:
        raise ChromBPNetPredictionError("one-hot sequence tensor has invalid geometry")
    return matrix[:, ::-1, ::-1]


def softmax(values: np.ndarray) -> np.ndarray:
    logits = np.asarray(values, dtype=np.float64)
    if logits.ndim != 2 or not np.isfinite(logits).all():
        raise ChromBPNetPredictionError("profile logits are invalid")
    shifted = logits - logits.max(axis=1, keepdims=True)
    exponentiated = np.exp(shifted)
    probabilities = exponentiated / exponentiated.sum(axis=1, keepdims=True)
    if not np.isfinite(probabilities).all():
        raise ChromBPNetPredictionError("profile probabilities are invalid")
    return probabilities.astype(np.float32)


def regional_mass_from_logscore(values: np.ndarray, transform: str) -> np.ndarray:
    logcounts = np.asarray(values, dtype=np.float64).reshape(-1)
    if not np.isfinite(logcounts).all() or np.any(logcounts > 50):
        raise ChromBPNetPredictionError("log-count prediction is invalid")
    if transform == "log1p_absolute":
        return np.maximum(np.expm1(logcounts), 0.0)
    if transform == "log_component":
        return np.exp(logcounts)
    raise ChromBPNetPredictionError("regional-mass transform is invalid")


def read_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise ChromBPNetPredictionError(f"{path.name} fields differ")
        return [dict(row) for row in reader]


def read_fold(path: Path) -> tuple[dict[str, list[str]], dict[str, str]]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if set(raw) != {"train", "valid", "test"} or any(
        not isinstance(raw[role], list) for role in raw
    ):
        raise ChromBPNetPredictionError("fold JSON fields differ")
    fold = {role: [str(value) for value in raw[role]] for role in raw}
    role_by_contig = {
        contig: role for role, contigs in fold.items() for contig in contigs
    }
    if set(role_by_contig) != set(PRIMARY_CONTIGS) or sum(
        len(values) for values in fold.values()
    ) != len(PRIMARY_CONTIGS):
        raise ChromBPNetPredictionError("fold JSON is not a primary partition")
    return fold, role_by_contig


def batches(values: Sequence[Any], batch_size: int) -> Iterable[Sequence[Any]]:
    for start in range(0, len(values), batch_size):
        yield values[start : start + batch_size]


def predict(
    *,
    model_path: Path,
    fasta_path: Path,
    windows_path: Path,
    fold_path: Path,
    output: Path,
    batch_size: int,
    seed: int,
    count_transform: str,
) -> dict[str, object]:
    if (
        output.exists()
        or batch_size < 1
        or seed < 0
        or count_transform not in {"log1p_absolute", "log_component"}
    ):
        raise ChromBPNetPredictionError("invalid output or batch size")
    for path in (model_path, fasta_path, windows_path, fold_path):
        if path.is_symlink() or not path.is_file():
            raise ChromBPNetPredictionError("inference input is not a regular file")

    import h5py
    import pyfaidx
    import tensorflow as tf

    from chrombpnet.training.utils.losses import multinomial_nll

    random.seed(seed)
    np.random.seed(seed)
    tf.random.set_seed(seed)
    _fold, role_by_contig = read_fold(fold_path)
    windows = [
        row
        for row in read_tsv(windows_path, CCRE_FIELDS)
        if role_by_contig.get(row["contig"]) in {"valid", "test"}
    ]
    if len(windows) != 32_000 or len({row["window_id"] for row in windows}) != len(
        windows
    ):
        raise ChromBPNetPredictionError("fixed held-window census differs")

    model = tf.keras.models.load_model(
        model_path,
        custom_objects={"multinomial_nll": multinomial_nll, "tf": tf},
        compile=False,
    )
    if model.input_shape != (None, 2114, 4):
        raise ChromBPNetPredictionError("ChromBPNet input geometry differs")
    output_shapes = [tuple(value) for value in model.output_shape]
    if output_shapes != [(None, 1000), (None, 1)]:
        raise ChromBPNetPredictionError("ChromBPNet output geometry differs")

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
    count_min = float("inf")
    count_max = float("-inf")
    role_counts = {"valid": 0, "test": 0}
    try:
        with h5py.File(profile_path, "x") as profiles, counts_path.open(
            "x", encoding="utf-8", newline=""
        ) as counts_handle:
            profile_values = profiles.create_dataset(
                "profile_probability",
                shape=(len(windows), 1000),
                dtype="f4",
                chunks=(min(batch_size, len(windows)), 1000),
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
                "masld-bench-chrombpnet-ccre-profile-probability-v1"
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
                    if start < 0 or end - start != 2114:
                        raise ChromBPNetPredictionError("input window geometry differs")
                    sequence = str(reference[row["contig"]][start:end]).upper()
                    if len(sequence) != 2114:
                        raise ChromBPNetPredictionError("reference window length differs")
                    sequences.append(sequence)
                forward = one_hot_dna(sequences)
                reverse = reverse_complement_one_hot(forward)
                forward_profile, forward_logcount_raw = model.predict(
                    forward, batch_size=batch_size, verbose=0
                )
                reverse_profile_raw, reverse_logcount_raw = model.predict(
                    reverse, batch_size=batch_size, verbose=0
                )
                forward_probability = softmax(forward_profile)
                reverse_probability = softmax(reverse_profile_raw)[:, ::-1]
                averaged_probability = (
                    0.5 * (forward_probability + reverse_probability)
                ).astype(np.float32)
                if not np.allclose(
                    averaged_probability.sum(axis=1), 1.0, rtol=0, atol=1e-6
                ):
                    raise ChromBPNetPredictionError(
                        "strand-averaged profile does not sum to one"
                    )
                forward_logcount = np.asarray(forward_logcount_raw).reshape(-1)
                reverse_logcount = np.asarray(reverse_logcount_raw).reshape(-1)
                averaged_count = 0.5 * (
                    regional_mass_from_logscore(forward_logcount, count_transform)
                    + regional_mass_from_logscore(reverse_logcount, count_transform)
                )
                profile_values[offset : offset + len(batch)] = averaged_probability
                for index, row in enumerate(batch):
                    role = role_by_contig[row["contig"]]
                    role_counts[role] += 1
                    value = float(averaged_count[index])
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
                            "forward_logcount": format(
                                float(forward_logcount[index]), ".17g"
                            ),
                            "reverse_logcount": format(
                                float(reverse_logcount[index]), ".17g"
                            ),
                            "strand_averaged_mass": format(value, ".17g"),
                        }
                    )
                offset += len(batch)
    finally:
        reference.close()
    if role_counts != {"valid": 16_000, "test": 16_000}:
        raise ChromBPNetPredictionError("held prediction role census differs")
    summary = {
        "schema_version": "masld-bench-chrombpnet-ccre-predictions-v1",
        "status": "pass",
        "windows": len(windows),
        "role_counts": role_counts,
        "input_length": 2114,
        "output_length": 1000,
        "batch_size": batch_size,
        "seed": seed,
        "reverse_complement_tta": True,
        "profile_strand_averaging": (
            "mean_probability_after_reverse_axis_restoration"
        ),
        "regional_mass_transform": count_transform,
        "regional_mass_semantics": (
            "predicted_count_after_log1p_inverse"
            if count_transform == "log1p_absolute"
            else "positive_sequence_component_mass_from_exp_log_component"
        ),
        "count_strand_averaging": "arithmetic_mean_in_positive_mass_space",
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
    parser.add_argument("--seed", type=int, default=20260824)
    parser.add_argument(
        "--count-transform",
        choices=("log1p_absolute", "log_component"),
        default="log1p_absolute",
    )
    arguments = parser.parse_args()
    result = predict(
        model_path=arguments.model,
        fasta_path=arguments.fasta,
        windows_path=arguments.windows,
        fold_path=arguments.fold,
        output=arguments.output,
        batch_size=arguments.batch_size,
        seed=arguments.seed,
        count_transform=arguments.count_transform,
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

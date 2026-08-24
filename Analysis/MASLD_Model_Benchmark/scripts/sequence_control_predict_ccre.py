#!/usr/bin/env python3
"""Emit outcome-blind fixed-cCRE predictions for a sequence control."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import random
from typing import Any, Iterable, Sequence

import numpy as np


MODEL_IDS = ("sequence_cnn_control", "sequence_transformer_control")
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
COUNT_FIELDS = (
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
STANDARDIZED_FIELDS = (
    "window_id",
    "selection_hash",
    "contig",
    "output_start",
    "output_end",
    "genomic_fold",
    "role",
    "stratum",
    "predicted",
)


class SequenceControlPredictionError(ValueError):
    """Raised when fixed-window control inference violates its contract."""


def one_hot_dna(sequences: Sequence[str]) -> np.ndarray:
    if not sequences or len({len(sequence) for sequence in sequences}) != 1:
        raise SequenceControlPredictionError("DNA sequence batch has invalid geometry")
    encoded = np.frombuffer("".join(sequences).encode("ascii"), dtype=np.uint8)
    lookup = np.full(256, -1, dtype=np.int8)
    for index, base in enumerate(b"ACGT"):
        lookup[base] = index
    indices = lookup[encoded]
    if np.any(indices < 0):
        raise SequenceControlPredictionError("DNA sequence batch contains ambiguity")
    shape = (len(sequences), len(sequences[0]))
    return np.eye(4, dtype=np.float32)[indices.reshape(shape)]


def reverse_complement_one_hot(values: np.ndarray) -> np.ndarray:
    matrix = np.asarray(values)
    if matrix.ndim != 3 or matrix.shape[2] != 4:
        raise SequenceControlPredictionError("one-hot tensor has invalid geometry")
    return matrix[:, ::-1, ::-1]


def softmax(values: np.ndarray) -> np.ndarray:
    logits = np.asarray(values, dtype=np.float64)
    if logits.ndim != 2 or not np.isfinite(logits).all():
        raise SequenceControlPredictionError("profile logits are invalid")
    shifted = logits - logits.max(axis=1, keepdims=True)
    exponentiated = np.exp(shifted)
    probabilities = exponentiated / exponentiated.sum(axis=1, keepdims=True)
    if not np.isfinite(probabilities).all():
        raise SequenceControlPredictionError("profile probabilities are invalid")
    return probabilities.astype(np.float32)


def mass_from_log1p(values: np.ndarray) -> np.ndarray:
    logcounts = np.asarray(values, dtype=np.float64).reshape(-1)
    if not np.isfinite(logcounts).all() or np.any(logcounts > 50):
        raise SequenceControlPredictionError("log-count prediction is invalid")
    return np.maximum(np.expm1(logcounts), 0.0)


def _read_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise SequenceControlPredictionError(f"{path.name} fields differ")
        return [dict(row) for row in reader]


def _read_fold(path: Path) -> dict[str, str]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if set(raw) != {"train", "valid", "test"} or any(
        not isinstance(raw[role], list) for role in raw
    ):
        raise SequenceControlPredictionError("fold JSON fields differ")
    role_by_contig = {
        str(contig): role for role in ("train", "valid", "test") for contig in raw[role]
    }
    if set(role_by_contig) != set(PRIMARY_CONTIGS) or sum(
        len(raw[role]) for role in raw
    ) != len(PRIMARY_CONTIGS):
        raise SequenceControlPredictionError("fold JSON is not a primary partition")
    return role_by_contig


def _batches(values: Sequence[Any], batch_size: int) -> Iterable[Sequence[Any]]:
    for start in range(0, len(values), batch_size):
        yield values[start : start + batch_size]


def _table_schema_sha256(fields: Sequence[str]) -> str:
    return sha256(("\t".join(fields) + "\n").encode("utf-8")).hexdigest()


def predict(
    *,
    model_id: str,
    model_path: Path,
    fasta_path: Path,
    windows_path: Path,
    fold_path: Path,
    output: Path,
    lineage: str,
    batch_size: int,
    seed: int,
    expected_windows_per_role: int,
) -> dict[str, object]:
    if (
        model_id not in MODEL_IDS
        or output.exists()
        or not lineage
        or batch_size < 1
        or seed < 0
        or expected_windows_per_role < 1
    ):
        raise SequenceControlPredictionError("invalid inference request")
    for path in (model_path, fasta_path, windows_path, fold_path):
        if path.is_symlink() or not path.is_file():
            raise SequenceControlPredictionError("inference input is not a regular file")

    import h5py
    import pyfaidx
    import tensorflow as tf

    from chrombpnet.training.utils.losses import multinomial_nll

    random.seed(seed)
    np.random.seed(seed)
    tf.random.set_seed(seed)
    role_by_contig = _read_fold(fold_path)
    windows = [
        row
        for row in _read_tsv(windows_path, CCRE_FIELDS)
        if role_by_contig.get(row["contig"]) in {"valid", "test"}
    ]
    expected_total = 2 * expected_windows_per_role
    if len(windows) != expected_total or len(
        {row["window_id"] for row in windows}
    ) != len(windows):
        raise SequenceControlPredictionError("fixed held-window census differs")

    model = tf.keras.models.load_model(
        model_path,
        custom_objects={"multinomial_nll": multinomial_nll, "tf": tf},
        compile=False,
    )
    if model.name != model_id or model.input_shape != (None, 2114, 4):
        raise SequenceControlPredictionError("sequence-control model identity differs")
    if [tuple(shape) for shape in model.output_shape] != [
        (None, 1000),
        (None, 1),
    ]:
        raise SequenceControlPredictionError("sequence-control output geometry differs")

    output.mkdir(mode=0o750)
    reference = pyfaidx.Fasta(
        str(fasta_path),
        as_raw=True,
        sequence_always_upper=True,
        rebuild=False,
        read_ahead=1_000_000,
    )
    profiles_path = output / "profile_probabilities.h5"
    counts_path = output / "regional_counts.tsv"
    standardized_path = output / "standardized_predictions.tsv"
    string_type = h5py.string_dtype(encoding="utf-8")
    role_counts = {"valid": 0, "test": 0}
    minimum_mass = float("inf")
    maximum_mass = float("-inf")
    try:
        with h5py.File(profiles_path, "x") as profiles, counts_path.open(
            "x", encoding="utf-8", newline=""
        ) as counts_handle, standardized_path.open(
            "x", encoding="utf-8", newline=""
        ) as standardized_handle:
            profile_values = profiles.create_dataset(
                "profile_probability",
                shape=(len(windows), 1000),
                dtype="f4",
                chunks=(min(batch_size, len(windows)), 1000),
                compression="gzip",
                compression_opts=1,
            )
            for field in ("window_id", "selection_hash"):
                profiles.create_dataset(
                    field,
                    data=np.asarray([row[field] for row in windows], dtype=object),
                    dtype=string_type,
                )
            profiles.attrs["schema_version"] = (
                "masld-bench-sequence-control-ccre-profile-v1"
            )
            profiles.attrs["model_id"] = model_id
            profiles.attrs["strand_averaging"] = (
                "mean_probability_after_reverse_axis_restoration"
            )
            counts_writer = csv.DictWriter(
                counts_handle, fieldnames=COUNT_FIELDS, delimiter="\t"
            )
            standardized_writer = csv.DictWriter(
                standardized_handle, fieldnames=STANDARDIZED_FIELDS, delimiter="\t"
            )
            counts_writer.writeheader()
            standardized_writer.writeheader()
            offset = 0
            for batch in _batches(windows, batch_size):
                sequences = []
                for row in batch:
                    start, end = int(row["input_start"]), int(row["input_end"])
                    if start < 0 or end - start != 2114:
                        raise SequenceControlPredictionError(
                            "input window geometry differs"
                        )
                    sequence = str(reference[row["contig"]][start:end]).upper()
                    if len(sequence) != 2114:
                        raise SequenceControlPredictionError(
                            "reference window length differs"
                        )
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
                    averaged_probability.sum(axis=1), 1.0, rtol=0, atol=1.0e-6
                ):
                    raise SequenceControlPredictionError(
                        "strand-averaged profile does not sum to one"
                    )
                forward_logcount = np.asarray(forward_count_raw).reshape(-1)
                reverse_logcount = np.asarray(reverse_count_raw).reshape(-1)
                averaged_mass = 0.5 * (
                    mass_from_log1p(forward_logcount)
                    + mass_from_log1p(reverse_logcount)
                )
                profile_values[offset : offset + len(batch)] = averaged_probability
                for index, row in enumerate(batch):
                    role = role_by_contig[row["contig"]]
                    role_counts[role] += 1
                    mass = float(averaged_mass[index])
                    minimum_mass = min(minimum_mass, mass)
                    maximum_mass = max(maximum_mass, mass)
                    common = {
                        "window_id": row["window_id"],
                        "selection_hash": row["selection_hash"],
                        "contig": row["contig"],
                        "output_start": row["output_start"],
                        "output_end": row["output_end"],
                        "genomic_fold": row["genomic_fold"],
                        "role": role,
                    }
                    counts_writer.writerow(
                        {
                            **common,
                            "ccre_class": row["ccre_class"],
                            "forward_logcount": format(
                                float(forward_logcount[index]), ".17g"
                            ),
                            "reverse_logcount": format(
                                float(reverse_logcount[index]), ".17g"
                            ),
                            "strand_averaged_mass": format(mass, ".17g"),
                        }
                    )
                    standardized_writer.writerow(
                        {
                            **common,
                            "stratum": lineage,
                            "predicted": format(mass, ".17g"),
                        }
                    )
                offset += len(batch)
    finally:
        reference.close()

    expected_role_counts = {
        "valid": expected_windows_per_role,
        "test": expected_windows_per_role,
    }
    if role_counts != expected_role_counts:
        raise SequenceControlPredictionError("held prediction role census differs")
    prediction_contract = {
        "schema_version": "masld-bench-sequence-control-prediction-contract-v1",
        "task_id": "rna_conditioned_atac",
        "model_id": model_id,
        "prediction_granularity": "fixed_ccre_window",
        "join_keys": ["window_id", "selection_hash"],
        "prediction_field": "predicted",
        "stratum": lineage,
        "donor_dependence": "sequence_only_donor_invariant",
        "evaluator_action": (
            "expand_each_window_prediction_across_the_eligible_held_donor_roster_"
            "then_assign_evaluator_salted_row_donor_and_block_hashes"
        ),
        "observed_atac_input_exposed": False,
        "model_inference_input_eligible": True,
        "ready_for_independent_evaluator_expansion": True,
        "prediction_bundle_created_by_model_adapter": False,
    }
    (output / "prediction_contract.json").write_text(
        json.dumps(prediction_contract, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )
    summary = {
        "schema_version": "masld-bench-sequence-control-ccre-predictions-v1",
        "status": "pass",
        "model_id": model_id,
        "lineage_id": lineage,
        "windows": len(windows),
        "role_counts": role_counts,
        "input_length": 2114,
        "output_length": 1000,
        "batch_size": batch_size,
        "seed": seed,
        "reverse_complement_tta": True,
        "regional_mass_semantics": "predicted_count_after_log1p_inverse",
        "count_strand_averaging": "arithmetic_mean_in_positive_mass_space",
        "minimum_predicted_count": minimum_mass,
        "maximum_predicted_count": maximum_mass,
        "observed_atac_input_exposed": False,
        "benchmark_metrics_calculated": False,
        "standardized_table_schema_sha256": _table_schema_sha256(
            STANDARDIZED_FIELDS
        ),
        "standardized_predictions_path": standardized_path.name,
        "regional_counts_path": counts_path.name,
        "profile_probabilities_path": profiles_path.name,
        "prediction_contract_path": "prediction_contract.json",
    }
    (output / "summary.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-id", choices=MODEL_IDS, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--windows", type=Path, required=True)
    parser.add_argument("--fold", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--lineage", required=True)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--seed", type=int, default=20260824)
    parser.add_argument("--expected-windows-per-role", type=int, default=16000)
    arguments = parser.parse_args()
    result = predict(
        model_id=arguments.model_id,
        model_path=arguments.model,
        fasta_path=arguments.fasta,
        windows_path=arguments.windows,
        fold_path=arguments.fold,
        output=arguments.output,
        lineage=arguments.lineage,
        batch_size=arguments.batch_size,
        seed=arguments.seed,
        expected_windows_per_role=arguments.expected_windows_per_role,
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

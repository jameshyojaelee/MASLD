#!/usr/bin/env python3
"""Frozen helpers for the exact project-trained BPNet control."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


INPUT_LENGTH = 2114
OUTPUT_LENGTH = 1000
DEFAULT_SEED = 20260824
PRIMARY_CONTIGS = tuple([f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY"])
FORBIDDEN_OUTCOME_TOKEN = "chrombpnet-hepatocyte-donor0-genomic0-outcomes"


class BPNetContractError(ValueError):
    """Raised when the BPNet matched-control requirements are not met."""


def architecture_parameters(counts_loss_weight: float) -> dict[str, object]:
    if not np.isfinite(counts_loss_weight) or counts_loss_weight <= 0:
        raise BPNetContractError("counts loss weight must be finite and positive")
    return {
        "input_len": INPUT_LENGTH,
        "output_profile_len": OUTPUT_LENGTH,
        "motif_module_params": {
            "filters": [64],
            "kernel_sizes": [21],
            "padding": "valid",
        },
        "syntax_module_params": {
            "num_dilation_layers": 8,
            "filters": 64,
            "kernel_size": 3,
            "padding": "valid",
            "pre_activation_residual_unit": True,
        },
        "profile_head_params": {
            "filters": 1,
            "kernel_size": 75,
            "padding": "valid",
        },
        "counts_head_params": {
            "units": [1],
            "dropouts": [0.0],
            "activations": ["linear"],
        },
        "profile_bias_module_params": {"kernel_sizes": [1]},
        "counts_bias_module_params": {},
        "loss_weights": [1.0, float(counts_loss_weight)],
        "counts_loss": "MSE",
    }


def single_task(
    *, bigwig: Path, peaks: Path, nonpeaks: Path, negative_ratio: float
) -> dict[int, dict[str, object]]:
    if not np.isfinite(negative_ratio) or negative_ratio < 0:
        raise BPNetContractError("negative ratio must be finite and nonnegative")
    return {
        0: {
            "signal": {"source": [str(bigwig.resolve(strict=True))]},
            "loci": {"source": [str(peaks.resolve(strict=True))]},
            "background_loci": {
                "source": [str(nonpeaks.resolve(strict=True))],
                "ratio": [float(negative_ratio)],
            },
            "bias": {"source": [], "smoothing": []},
        }
    }


def validate_model_inputs(paths: Iterable[Path]) -> None:
    for path in paths:
        if path.is_symlink():
            raise BPNetContractError(f"BPNet input is a symlink: {path}")
        resolved = path.resolve(strict=True)
        if not resolved.is_file():
            raise BPNetContractError(f"BPNet input is not a regular file: {path}")
        if FORBIDDEN_OUTCOME_TOKEN in resolved.as_posix() or "outcomes" in resolved.name:
            raise BPNetContractError("evaluator-only outcome artifact is forbidden")


def read_fold(path: Path, *, require_primary_partition: bool = True) -> dict[str, list[str]]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if set(raw) != {"train", "valid", "test"}:
        raise BPNetContractError("fold JSON fields differ")
    fold: dict[str, list[str]] = {}
    observed: list[str] = []
    for role in ("train", "valid", "test"):
        values = raw[role]
        if not isinstance(values, list) or not values:
            raise BPNetContractError(f"fold role {role!r} is invalid")
        fold[role] = [str(value) for value in values]
        observed.extend(fold[role])
    if len(observed) != len(set(observed)):
        raise BPNetContractError("fold roles overlap")
    if require_primary_partition and set(observed) != set(PRIMARY_CONTIGS):
        raise BPNetContractError("fold JSON is not a primary-contig partition")
    return fold


def normalize_narrowpeak_fields(
    fields: Sequence[str], *, row_number: int
) -> tuple[str, ...]:
    if len(fields) < 3:
        raise BPNetContractError(f"BED row {row_number} has fewer than three fields")
    chrom = fields[0]
    try:
        start, end = int(fields[1]), int(fields[2])
    except ValueError as error:
        raise BPNetContractError(f"BED row {row_number} has invalid coordinates") from error
    if start < 0 or end <= start:
        raise BPNetContractError(f"BED row {row_number} has invalid interval")
    if len(fields) >= 10:
        try:
            summit = int(fields[9])
        except ValueError as error:
            raise BPNetContractError(f"BED row {row_number} has invalid summit") from error
        if summit < 0 or summit >= end - start:
            raise BPNetContractError(f"BED row {row_number} summit is outside interval")
        return (
            chrom,
            str(start),
            str(end),
            fields[3] or f"bpnet_{row_number}",
            fields[4] if fields[4] else "0",
            fields[5] if fields[5] else ".",
            fields[6] if fields[6] else "0",
            fields[7] if fields[7] else "-1",
            fields[8] if fields[8] else "-1",
            str(summit),
        )
    return (
        chrom,
        str(start),
        str(end),
        fields[3] if len(fields) > 3 and fields[3] else f"bpnet_{row_number}",
        "0",
        ".",
        "0",
        "-1",
        "-1",
        str((end - start) // 2),
    )


def prepare_narrowpeak(
    source: Path,
    destination: Path,
    *,
    allowed_contigs: set[str],
    chrom_sizes: dict[str, int] | None = None,
    required_flank: int = 0,
    required_flank_by_contig: dict[str, int] | None = None,
) -> int:
    if destination.exists():
        raise BPNetContractError(f"refusing to overwrite {destination}")
    if (
        required_flank < 0
        or (required_flank and chrom_sizes is None)
        or (required_flank_by_contig is not None and chrom_sizes is None)
        or (
            required_flank_by_contig is not None
            and (
                set(required_flank_by_contig) != allowed_contigs
                or any(value < 0 for value in required_flank_by_contig.values())
            )
        )
    ):
        raise BPNetContractError("narrowPeak boundary contract is invalid")
    count = 0
    with source.open("r", encoding="utf-8") as source_handle, destination.open(
        "x", encoding="utf-8"
    ) as destination_handle:
        for row_number, line in enumerate(source_handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            normalized = normalize_narrowpeak_fields(fields, row_number=row_number)
            if normalized[0] not in allowed_contigs:
                continue
            if chrom_sizes is not None:
                contig = normalized[0]
                if contig not in chrom_sizes or chrom_sizes[contig] <= 0:
                    raise BPNetContractError(
                        f"missing or invalid chromosome size for {contig}"
                    )
                center = int(normalized[1]) + int(normalized[9])
                flank = (
                    required_flank_by_contig[contig]
                    if required_flank_by_contig is not None
                    else required_flank
                )
                if (
                    center - flank < 0
                    or center + flank > chrom_sizes[contig]
                ):
                    continue
            destination_handle.write("\t".join(normalized) + "\n")
            count += 1
    if count == 0:
        raise BPNetContractError(f"no admitted rows in {source}")
    return count


def read_chrom_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open("r", encoding="utf-8") as handle:
        for row_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2 or not fields[0] or fields[0] in sizes:
                raise BPNetContractError(
                    f"chromosome-size row {row_number} is invalid"
                )
            try:
                size = int(fields[1])
            except ValueError as error:
                raise BPNetContractError(
                    f"chromosome-size row {row_number} is invalid"
                ) from error
            if size <= 0:
                raise BPNetContractError(
                    f"chromosome-size row {row_number} is invalid"
                )
            sizes[fields[0]] = size
    if not sizes:
        raise BPNetContractError("chromosome-size file is empty")
    return sizes


def one_hot_dna(sequences: Sequence[str]) -> np.ndarray:
    if not sequences or len({len(sequence) for sequence in sequences}) != 1:
        raise BPNetContractError("DNA sequence batch has invalid geometry")
    encoded = np.frombuffer("".join(sequences).encode("ascii"), dtype=np.uint8)
    lookup = np.full(256, -1, dtype=np.int8)
    for index, base in enumerate(b"ACGT"):
        lookup[base] = index
    indices = lookup[encoded]
    if np.any(indices < 0):
        raise BPNetContractError("DNA sequence batch contains an ambiguous base")
    shape = (len(sequences), len(sequences[0]))
    return np.eye(4, dtype=np.float32)[indices.reshape(shape)]


def reverse_complement_one_hot(values: np.ndarray) -> np.ndarray:
    matrix = np.asarray(values)
    if matrix.ndim != 3 or matrix.shape[2] != 4:
        raise BPNetContractError("one-hot sequence tensor has invalid geometry")
    return matrix[:, ::-1, ::-1].copy()


def softmax(values: np.ndarray) -> np.ndarray:
    logits = np.asarray(values, dtype=np.float64)
    if logits.ndim == 3 and logits.shape[2] == 1:
        logits = logits[:, :, 0]
    if logits.ndim != 2 or not np.isfinite(logits).all():
        raise BPNetContractError("profile logits are invalid")
    shifted = logits - logits.max(axis=1, keepdims=True)
    exponentiated = np.exp(shifted)
    probabilities = exponentiated / exponentiated.sum(axis=1, keepdims=True)
    if not np.isfinite(probabilities).all():
        raise BPNetContractError("profile probabilities are invalid")
    return probabilities.astype(np.float32)

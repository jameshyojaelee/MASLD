#!/usr/bin/env python3
"""Shared token-coordinate and fold-safe projection helpers for GSE281364 DNA LMs."""

from __future__ import annotations

import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np


ALLELES = ("REF", "ALT", "REF_RC", "ALT_RC")
MODELS = ("dnabert2", "nucleotide_transformer", "hyenadna")
INPUT_LENGTH = 4_096
NT_PHASES = 6
NT_CONTEXT_BP = 4_086
NT_KMER = 6
PROJECTION_WIDTH = 256
STD_FLOOR = 1e-6
EIGENVALUE_ABSOLUTE_FLOOR = 1e-8
EIGENVALUE_RELATIVE_FLOOR = 1e-6


class NativeContractError(ValueError):
    """Raised when a token-coordinate or projection requirement differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def sequence_digest(sequence: str) -> str:
    return sha256(sequence.encode("ascii")).hexdigest()


def read_fixture(fixture: Path) -> tuple[list[dict[str, str]], dict[str, str]]:
    with (fixture / "fixture/sequence_manifest.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    records: dict[str, list[str]] = {}
    name: str | None = None
    with gzip.open(
        fixture / "fixture/common_4096.alleles.fa.gz", "rt", encoding="ascii"
    ) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                name = line[1:]
                if not name or name in records:
                    raise NativeContractError("fixture FASTA header differs")
                records[name] = []
            elif name is None:
                raise NativeContractError("fixture FASTA sequence precedes header")
            else:
                records[name].append(line)
    fasta = {record_id: "".join(parts) for record_id, parts in records.items()}
    if (
        len(rows) != 1_033
        or len(fasta) != 4_132
        or len({row["outer_locus_sequence_group_id"] for row in rows}) != 1_033
        or {row["outer_fold"] for row in rows} != {"0", "1", "2", "3", "4"}
    ):
        raise NativeContractError("fixture census or fold coverage differs")
    hash_fields = {
        "REF": "reference_sequence_sha256",
        "ALT": "alternative_sequence_sha256",
        "REF_RC": "reverse_complement_reference_sha256",
        "ALT_RC": "reverse_complement_alternative_sha256",
    }
    for row in rows:
        prefix = row["fasta_record_prefix"]
        for allele, field in hash_fields.items():
            sequence = fasta.get(f"{prefix}|{allele}")
            if (
                sequence is None
                or len(sequence) != INPUT_LENGTH
                or set(sequence) - set("ACGT")
                or sequence_digest(sequence) != row[field]
            ):
                raise NativeContractError("fixture sequence identity differs")
    return rows, fasta


def dnabert_encoding(tokenizer: object, sequence: str) -> tuple[list[int], list[tuple[int, int]]]:
    encoded = tokenizer(
        sequence,
        add_special_tokens=True,
        return_attention_mask=True,
        return_offsets_mapping=True,
        truncation=False,
    )
    ids = [int(value) for value in encoded["input_ids"]]
    mask = [int(value) for value in encoded["attention_mask"]]
    offsets = [(int(start), int(end)) for start, end in encoded["offset_mapping"]]
    if (
        len(ids) != len(mask)
        or len(ids) != len(offsets)
        or any(value != 1 for value in mask)
        or offsets[0] != (0, 0)
        or offsets[-1] != (0, 0)
        or any(value >= 4_096 for value in ids)
    ):
        raise NativeContractError("DNABERT-2 token structure differs")
    cursor = 0
    for start, end in offsets[1:-1]:
        if start != cursor or end <= start:
            raise NativeContractError("DNABERT-2 BPE offsets do not tile bases")
        cursor = end
    if cursor != len(sequence):
        raise NativeContractError("DNABERT-2 BPE offsets do not cover input")
    return ids, offsets


def affected_interval(
    reference_offsets: Sequence[tuple[int, int]],
    alternative_offsets: Sequence[tuple[int, int]],
    variant_index0: int,
) -> tuple[int, int]:
    start, end = variant_index0, variant_index0 + 1
    changed = True
    while changed:
        changed = False
        for offsets in (reference_offsets, alternative_offsets):
            for token_start, token_end in offsets:
                if token_end > start and token_start < end:
                    updated = min(start, token_start), max(end, token_end)
                    if updated != (start, end):
                        start, end = updated
                        changed = True
    return start, end


def overlapping_indices(
    offsets: Sequence[tuple[int, int]], start: int, end: int
) -> list[int]:
    return [
        index
        for index, (token_start, token_end) in enumerate(offsets)
        if token_end > start and token_start < end
    ]


def nt_phase_spans(length: int = INPUT_LENGTH) -> list[tuple[int, int]]:
    spans = [(phase, phase + NT_CONTEXT_BP) for phase in range(NT_PHASES)]
    if (
        length != INPUT_LENGTH
        or NT_CONTEXT_BP % NT_KMER
        or any(start < 0 or end > length for start, end in spans)
        or any(start + (length - end) != 10 for start, end in spans)
    ):
        raise NativeContractError("Nucleotide Transformer phase geometry differs")
    return spans


def nt_expected_tokens(sequence: str, phase: int) -> list[str]:
    start, end = nt_phase_spans(len(sequence))[phase]
    return ["<cls>"] + [
        sequence[index : index + NT_KMER]
        for index in range(start, end, NT_KMER)
    ]


def nt_tokenize(tokenizer: object, sequence: str) -> list[list[str]]:
    phased = [sequence[start:end] for start, end in nt_phase_spans(len(sequence))]
    encoded = tokenizer(
        phased,
        add_special_tokens=True,
        padding=False,
        truncation=False,
        return_attention_mask=True,
    )
    token_rows: list[list[str]] = []
    for phase, (ids, mask) in enumerate(
        zip(encoded["input_ids"], encoded["attention_mask"])
    ):
        active = [int(token_id) for token_id, keep in zip(ids, mask) if int(keep)]
        if len(active) != 682 or max(active) >= 4_105:
            raise NativeContractError("Nucleotide Transformer token IDs differ")
        tokens = tokenizer.convert_ids_to_tokens(active)
        if tokens != nt_expected_tokens(sequence, phase):
            raise NativeContractError("Nucleotide Transformer six-phase tokens differ")
        token_rows.append(tokens)
    return token_rows


def nt_pool_indices(phase: int, pool_start0: int, pool_end0: int) -> list[int]:
    start, _ = nt_phase_spans()[phase]
    indices = []
    for token_index, base_start in enumerate(
        range(start, start + NT_CONTEXT_BP, NT_KMER), start=1
    ):
        if base_start < pool_end0 and base_start + NT_KMER > pool_start0:
            indices.append(token_index)
    if not indices:
        raise NativeContractError("Nucleotide Transformer pool has no tokens")
    return indices


def validate_nt_pair(
    tokenizer: object,
    reference: str,
    alternative: str,
    variant_index0: int,
) -> tuple[list[list[str]], list[list[str]]]:
    reference_tokens = nt_tokenize(tokenizer, reference)
    alternative_tokens = nt_tokenize(tokenizer, alternative)
    for phase, (left, right) in enumerate(zip(reference_tokens, alternative_tokens)):
        differences = [index for index, pair in enumerate(zip(left, right)) if pair[0] != pair[1]]
        if len(differences) != 1:
            raise NativeContractError("Nucleotide Transformer allele token difference differs")
        token_index = differences[0]
        token_start = phase + (token_index - 1) * NT_KMER
        if not token_start <= variant_index0 < token_start + NT_KMER:
            raise NativeContractError("Nucleotide Transformer variant-token mapping differs")
    return reference_tokens, alternative_tokens


def fit_projection(
    embeddings: np.ndarray,
    folds: np.ndarray,
    held_out_fold: int,
    width: int = PROJECTION_WIDTH,
) -> dict[str, np.ndarray]:
    values = np.asarray(embeddings, dtype=np.float64)
    fold_values = np.asarray(folds, dtype=np.int64)
    if (
        values.ndim != 3
        or values.shape[1] != len(ALLELES)
        or fold_values.shape != (values.shape[0],)
        or held_out_fold not in range(5)
        or set(fold_values.tolist()) != set(range(5))
        or width < 1
        or width > values.shape[2]
        or not np.isfinite(values).all()
    ):
        raise NativeContractError("projection fit request differs")
    training = values[fold_values != held_out_fold].reshape(-1, values.shape[2])
    mean = training.mean(axis=0)
    standard_deviation = training.std(axis=0)
    standard_deviation = np.maximum(standard_deviation, STD_FLOOR)
    normalized = (training - mean) / standard_deviation
    _, singular_values, right_vectors = np.linalg.svd(normalized, full_matrices=False)
    components = right_vectors[:width].copy()
    for index in range(width):
        pivot = int(np.argmax(np.abs(components[index])))
        if components[index, pivot] < 0:
            components[index] *= -1
    eigenvalues = singular_values[:width] ** 2 / max(1, training.shape[0] - 1)
    eigenvalue_floor = max(
        EIGENVALUE_ABSOLUTE_FLOOR,
        EIGENVALUE_RELATIVE_FLOOR * float(eigenvalues[0]),
    )
    whitening_scale = np.sqrt(np.maximum(eigenvalues, eigenvalue_floor))
    return {
        "mean": mean,
        "standard_deviation": standard_deviation,
        "components": components,
        "whitening_scale": whitening_scale,
        "held_out_fold": np.asarray(held_out_fold, dtype=np.int64),
        "training_elements": np.asarray(np.sum(fold_values != held_out_fold), dtype=np.int64),
    }


def apply_projection(
    embeddings: np.ndarray, parameters: Mapping[str, np.ndarray]
) -> np.ndarray:
    values = np.asarray(embeddings, dtype=np.float64)
    mean = np.asarray(parameters["mean"], dtype=np.float64)
    standard_deviation = np.asarray(parameters["standard_deviation"], dtype=np.float64)
    components = np.asarray(parameters["components"], dtype=np.float64)
    whitening_scale = np.asarray(parameters["whitening_scale"], dtype=np.float64)
    if (
        values.ndim != 3
        or values.shape[1] != len(ALLELES)
        or values.shape[2] != mean.shape[0]
        or standard_deviation.shape != mean.shape
        or components.shape[1] != mean.shape[0]
        or whitening_scale.shape != (components.shape[0],)
    ):
        raise NativeContractError("projection transform request differs")
    normalized = (values - mean) / standard_deviation
    projected = normalized @ components.T
    projected /= whitening_scale
    if not np.isfinite(projected).all():
        raise NativeContractError("projected embeddings are not finite")
    return projected.astype(np.float32)


def head_features(projected_alleles: np.ndarray) -> np.ndarray:
    values = np.asarray(projected_alleles, dtype=np.float32)
    if values.ndim != 3 or values.shape[1] != len(ALLELES):
        raise NativeContractError("head feature request differs")
    reference = (values[:, 0] + values[:, 2]) / 2.0
    alternative = (values[:, 1] + values[:, 3]) / 2.0
    difference = alternative - reference
    features = np.concatenate(
        (reference, alternative, difference, np.abs(difference)), axis=1
    )
    if features.shape[1] != 4 * values.shape[2] or not np.isfinite(features).all():
        raise NativeContractError("head features differ")
    return features


def validate_contract_config(config: Mapping[str, object]) -> None:
    if (
        config.get("schema_version")
        != "masld-bench-gse281364-dna-lm-native-contract-v1"
        or config.get("dataset_id") != "gse281364"
        or tuple(config.get("models", {})) != MODELS
        or any(config.get("firewall", {}).values())
    ):
        raise NativeContractError("native contract configuration differs")
    projection = config.get("normalization_projection", {})
    if (
        projection.get("projection_width") != PROJECTION_WIDTH
        or projection.get("head_input_width") != 4 * PROJECTION_WIDTH
        or projection.get("head_fit") is not False
        or projection.get("reporter_outcomes_read") is not False
    ):
        raise NativeContractError("normalization/projection configuration differs")


def load_config(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    validate_contract_config(value)
    return value

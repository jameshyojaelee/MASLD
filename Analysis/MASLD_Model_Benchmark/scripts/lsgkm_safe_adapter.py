#!/usr/bin/env python3
"""Fail-closed adapters for project-generated LS-GKM models and scores."""

from __future__ import annotations

from hashlib import sha256
import gzip
import math
from pathlib import Path
import re
from typing import Iterable, Mapping, Sequence


ALPHABET = frozenset("ACGT")
COMPLEMENT = str.maketrans("ACGT", "TGCA")
EXPECTED_HEADER_KEYS = (
    "svm_type",
    "kernel_type",
    "L",
    "k",
    "d",
    "norc",
    "nr_class",
    "total_sv",
    "rho",
    "label",
    "nr_sv",
    "SV",
)
MAX_FASTA_BYTES = 256 * 1024 * 1024
MAX_MODEL_COMPRESSED_BYTES = 2 * 1024 * 1024 * 1024
MAX_MODEL_UNCOMPRESSED_BYTES = 8 * 1024 * 1024 * 1024
MAX_PREDICTION_BYTES = 2 * 1024 * 1024 * 1024
_IDENTIFIER = re.compile(r"^[!-~]+$")


class LSGKMSafetyError(RuntimeError):
    """Raised before unsafe or requirements-incompatible bytes reach native code."""


def file_sha256(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def reverse_complement(sequence: str) -> str:
    _validate_sequence(sequence, label="sequence", maximum_length=2047)
    return sequence.translate(COMPLEMENT)[::-1]


def canonical_kmer(sequence: str) -> str:
    _validate_sequence(sequence, label="k-mer", maximum_length=2047)
    reverse = reverse_complement(sequence)
    return min(sequence, reverse)


def _validate_sequence(
    sequence: str,
    *,
    label: str,
    maximum_length: int,
    exact_length: int | None = None,
) -> None:
    if not sequence:
        raise LSGKMSafetyError(f"{label} is empty")
    if len(sequence) > maximum_length:
        raise LSGKMSafetyError(f"{label} exceeds {maximum_length} bp")
    if exact_length is not None and len(sequence) != exact_length:
        raise LSGKMSafetyError(f"{label} is not exactly {exact_length} bp")
    if set(sequence) - ALPHABET:
        raise LSGKMSafetyError(f"{label} is not uppercase A/C/G/T")


def validate_fasta(
    path: str | Path,
    *,
    expected_length: int | None = None,
    maximum_length: int = 2047,
    maximum_records: int = 10_000_000,
    maximum_bytes: int = MAX_FASTA_BYTES,
) -> list[tuple[str, str]]:
    """Parse a plain FASTA with the stricter rules required before LS-GKM."""

    source = Path(path)
    if source.is_symlink() or not source.is_file():
        raise LSGKMSafetyError("FASTA must be a regular non-symlink file")
    if source.stat().st_size > maximum_bytes:
        raise LSGKMSafetyError("FASTA exceeds the configured byte bound")
    records: list[tuple[str, str]] = []
    seen: set[str] = set()
    identifier: str | None = None
    pieces: list[str] = []

    def finish() -> None:
        nonlocal identifier, pieces
        if identifier is None:
            return
        sequence = "".join(pieces)
        _validate_sequence(
            sequence,
            label=f"FASTA record {identifier}",
            maximum_length=maximum_length,
            exact_length=expected_length,
        )
        records.append((identifier, sequence))
        if len(records) > maximum_records:
            raise LSGKMSafetyError("FASTA record count exceeds configured bound")

    try:
        with source.open("r", encoding="ascii", errors="strict", newline="") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\r\n")
                if line.startswith(">"):
                    finish()
                    identifier = line[1:]
                    pieces = []
                    if (
                        not identifier
                        or len(identifier) > 2047
                        or _IDENTIFIER.fullmatch(identifier) is None
                        or identifier in seen
                    ):
                        raise LSGKMSafetyError("invalid or duplicate FASTA identifier")
                    seen.add(identifier)
                else:
                    if identifier is None:
                        raise LSGKMSafetyError("FASTA sequence precedes its identifier")
                    if not line:
                        raise LSGKMSafetyError("blank FASTA sequence line")
                    pieces.append(line)
    except (OSError, UnicodeError) as error:
        raise LSGKMSafetyError(f"FASTA cannot be read safely: {error}") from error
    finish()
    if not records:
        raise LSGKMSafetyError("FASTA contains no records")
    lengths = {len(sequence) for _, sequence in records}
    if len(lengths) != 1:
        raise LSGKMSafetyError("FASTA contains mixed sequence lengths")
    return records


def _read_bounded_model(path: Path) -> tuple[str, bool]:
    if path.is_symlink() or not path.is_file():
        raise LSGKMSafetyError("model must be a regular non-symlink file")
    size = path.stat().st_size
    if size <= 0 or size > MAX_MODEL_COMPRESSED_BYTES:
        raise LSGKMSafetyError("model compressed/plain byte bound failed")
    with path.open("rb") as handle:
        magic = handle.read(2)
    compressed = magic == b"\x1f\x8b"
    opener = gzip.open if compressed else Path.open
    chunks: list[bytes] = []
    total = 0
    try:
        if compressed:
            handle_context = gzip.open(path, "rb")
        else:
            handle_context = path.open("rb")
        with handle_context as handle:
            while True:
                block = handle.read(1024 * 1024)
                if not block:
                    break
                total += len(block)
                if total > MAX_MODEL_UNCOMPRESSED_BYTES:
                    raise LSGKMSafetyError("model decompressed byte bound failed")
                chunks.append(block)
    except (OSError, EOFError, gzip.BadGzipFile) as error:
        raise LSGKMSafetyError(f"model gzip/plain read failed: {error}") from error
    try:
        text = b"".join(chunks).decode("ascii", errors="strict")
    except UnicodeDecodeError as error:
        raise LSGKMSafetyError("model is not strict ASCII") from error
    if not text.endswith("\n") or "\r" in text or "\x00" in text:
        raise LSGKMSafetyError("model newline or control-byte contract failed")
    return text, compressed


def validate_project_model(
    path: str | Path,
    *,
    support_vector_length: int = 300,
    maximum_support_vectors: int = 10_000_000,
) -> dict[str, object]:
    """Validate the exact binary C-SVC model grammar before native loading."""

    source = Path(path)
    text, compressed = _read_bounded_model(source)
    lines = text.splitlines()
    if len(lines) < len(EXPECTED_HEADER_KEYS) + 1 or any(not line for line in lines):
        raise LSGKMSafetyError("model is empty, truncated, or contains blank lines")
    header_lines = lines[: len(EXPECTED_HEADER_KEYS)]
    observed_keys = tuple(line.split(" ", 1)[0] for line in header_lines)
    if observed_keys != EXPECTED_HEADER_KEYS:
        raise LSGKMSafetyError("model header order or key set differs")
    expected_singletons = {
        "svm_type": "c_svc",
        "kernel_type": "gkm_esttrunc",
        "L": "11",
        "k": "7",
        "d": "3",
        "norc": "0",
        "nr_class": "2",
    }
    parsed: dict[str, list[str]] = {}
    for line in header_lines:
        fields = line.split()
        parsed[fields[0]] = fields[1:]
    for key, value in expected_singletons.items():
        if parsed.get(key) != [value]:
            raise LSGKMSafetyError(f"model header differs for {key}")
    if parsed.get("SV"):
        raise LSGKMSafetyError("SV header must not carry values")
    try:
        total_sv = int(parsed["total_sv"][0])
        rho = float(parsed["rho"][0])
        labels = tuple(int(value) for value in parsed["label"])
        nr_sv = tuple(int(value) for value in parsed["nr_sv"])
    except (KeyError, IndexError, ValueError) as error:
        raise LSGKMSafetyError("model numeric header cannot be parsed") from error
    if (
        len(parsed["total_sv"]) != 1
        or len(parsed["rho"]) != 1
        or total_sv <= 0
        or total_sv > maximum_support_vectors
        or not math.isfinite(rho)
        or labels != (1, -1)
        or len(nr_sv) != 2
        or any(value <= 0 for value in nr_sv)
        or sum(nr_sv) != total_sv
    ):
        raise LSGKMSafetyError("model class, rho, or support-vector census differs")
    support_lines = lines[len(EXPECTED_HEADER_KEYS) :]
    if len(support_lines) != total_sv:
        raise LSGKMSafetyError("model support-vector row count differs")
    coefficient_min = math.inf
    coefficient_max = -math.inf
    for index, line in enumerate(support_lines):
        fields = line.split()
        if len(fields) != 2:
            raise LSGKMSafetyError("support-vector row grammar differs")
        try:
            coefficient = float(fields[0])
        except ValueError as error:
            raise LSGKMSafetyError("support-vector coefficient is invalid") from error
        if not math.isfinite(coefficient) or coefficient == 0.0:
            raise LSGKMSafetyError("support-vector coefficient is zero or non-finite")
        _validate_sequence(
            fields[1],
            label=f"support vector {index}",
            maximum_length=2047,
            exact_length=support_vector_length,
        )
        coefficient_min = min(coefficient_min, coefficient)
        coefficient_max = max(coefficient_max, coefficient)
    return {
        "path": source.as_posix(),
        "sha256": file_sha256(source),
        "size_bytes": source.stat().st_size,
        "gzip": compressed,
        "svm_type": "c_svc",
        "kernel_type": "gkm_esttrunc",
        "L": 11,
        "k": 7,
        "d": 3,
        "norc": 0,
        "nr_class": 2,
        "total_sv": total_sv,
        "rho": rho,
        "label": list(labels),
        "nr_sv": list(nr_sv),
        "support_vector_length": support_vector_length,
        "coefficient_min": coefficient_min,
        "coefficient_max": coefficient_max,
    }


def parse_predictions(
    path: str | Path,
    *,
    expected_identifiers: Sequence[str] | None = None,
) -> dict[str, float]:
    source = Path(path)
    if source.is_symlink() or not source.is_file():
        raise LSGKMSafetyError("prediction output must be a regular file")
    if source.stat().st_size > MAX_PREDICTION_BYTES:
        raise LSGKMSafetyError("prediction output exceeds byte bound")
    values: dict[str, float] = {}
    order: list[str] = []
    try:
        with source.open("r", encoding="ascii", errors="strict", newline="") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\r\n")
                fields = line.split("\t")
                if len(fields) != 2 or not fields[0] or fields[0] in values:
                    raise LSGKMSafetyError("prediction row grammar or identity differs")
                if _IDENTIFIER.fullmatch(fields[0]) is None:
                    raise LSGKMSafetyError("prediction identifier is invalid")
                try:
                    score = float(fields[1])
                except ValueError as error:
                    raise LSGKMSafetyError("prediction score is invalid") from error
                if not math.isfinite(score):
                    raise LSGKMSafetyError("prediction score is non-finite")
                values[fields[0]] = score
                order.append(fields[0])
    except (OSError, UnicodeError) as error:
        raise LSGKMSafetyError(f"prediction output cannot be read safely: {error}") from error
    if not values:
        raise LSGKMSafetyError("prediction output is empty")
    if expected_identifiers is not None and order != list(expected_identifiers):
        raise LSGKMSafetyError("prediction row order or identifier set differs")
    return values


def direct_gkmsvm_alt_minus_ref(
    *, reference_score: float, alternative_score: float
) -> float:
    if not math.isfinite(reference_score) or not math.isfinite(alternative_score):
        raise LSGKMSafetyError("direct gkm-SVM score is non-finite")
    return alternative_score - reference_score


def deltasvm_scores(
    reference_context: str,
    alternative_context: str,
    canonical_11mer_scores: Mapping[str, float],
) -> dict[str, float]:
    """Return native REF-minus-ALT and one-time-converted ALT-minus-REF scores."""

    _validate_sequence(
        reference_context,
        label="deltaSVM REF context",
        maximum_length=21,
        exact_length=21,
    )
    _validate_sequence(
        alternative_context,
        label="deltaSVM ALT context",
        maximum_length=21,
        exact_length=21,
    )
    differences = [
        index
        for index, (reference, alternative) in enumerate(
            zip(reference_context, alternative_context, strict=True)
        )
        if reference != alternative
    ]
    if differences != [10]:
        raise LSGKMSafetyError("deltaSVM v1 admits one center SNP only")

    def total(sequence: str) -> float:
        values: list[float] = []
        for start in range(11):
            key = canonical_kmer(sequence[start : start + 11])
            try:
                value = float(canonical_11mer_scores[key])
            except (KeyError, TypeError, ValueError) as error:
                raise LSGKMSafetyError(f"missing or invalid 11-mer score: {key}") from error
            if not math.isfinite(value):
                raise LSGKMSafetyError(f"non-finite 11-mer score: {key}")
            values.append(value)
        return math.fsum(values)

    reference_sum = total(reference_context)
    alternative_sum = total(alternative_context)
    native = reference_sum - alternative_sum
    return {
        "reference_11mer_sum": reference_sum,
        "alternative_11mer_sum": alternative_sum,
        "native_ref_minus_alt": native,
        "canonical_alt_minus_ref": -native,
    }


def canonical_score_table(rows: Iterable[tuple[str, float]]) -> dict[str, float]:
    """Normalize a scored nonredundant-k-mer table without averaging duplicates."""

    result: dict[str, float] = {}
    for sequence, raw_score in rows:
        key = canonical_kmer(sequence)
        score = float(raw_score)
        if not math.isfinite(score):
            raise LSGKMSafetyError("k-mer score is non-finite")
        if key in result and result[key] != score:
            raise LSGKMSafetyError("reverse-complement k-mer scores differ")
        result[key] = score
    if not result:
        raise LSGKMSafetyError("k-mer score table is empty")
    return result

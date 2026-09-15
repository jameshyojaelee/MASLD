"""Deterministic hashing primitives for benchmark requirements.

Only JSON-compatible values are hashable.  The conversion is intentionally
strict: an unordered set, a non-finite float, or an opaque Python object must
never acquire a process-dependent benchmark identity.
"""

from __future__ import annotations

from dataclasses import fields, is_dataclass
from enum import Enum
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence


class HashingError(ValueError):
    """Raised when a value cannot be represented by canonical JSON."""


def canonicalize(value: Any) -> Any:
    """Return a deterministic, JSON-compatible representation of ``value``."""

    if is_dataclass(value) and not isinstance(value, type):
        return {
            field.name: canonicalize(getattr(value, field.name))
            for field in fields(value)
            if not field.name.startswith("_")
        }
    if isinstance(value, Enum):
        return canonicalize(value.value)
    if isinstance(value, Path):
        return value.as_posix()
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise HashingError("non-finite floats are not valid canonical JSON")
        return value
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) for key in value):
            raise HashingError("canonical JSON object keys must be strings")
        return {
            key: canonicalize(value[key])
            for key in sorted(value)
        }
    if isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray, memoryview)
    ):
        return [canonicalize(item) for item in value]
    raise HashingError(
        f"unsupported canonical JSON value: {type(value).__name__}"
    )


def canonical_json(value: Any) -> str:
    """Serialize ``value`` as RFC-8259-compatible canonical JSON."""

    try:
        return json.dumps(
            canonicalize(value),
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
    except (TypeError, ValueError) as error:
        raise HashingError(f"value is not canonical-JSON serializable: {error}") from error


def canonical_json_bytes(value: Any) -> bytes:
    return canonical_json(value).encode("utf-8")


def sha256_bytes(payload: bytes) -> str:
    if not isinstance(payload, bytes):
        raise TypeError("sha256_bytes requires bytes")
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: str | Path, *, chunk_size: int = 8 * 1024 * 1024) -> str:
    source = Path(path)
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    if not source.is_file():
        raise FileNotFoundError(f"not a regular file: {source}")
    digest = hashlib.sha256()
    with source.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_sha256(value: Any) -> str:
    return sha256_bytes(canonical_json_bytes(value))


def is_sha256(value: object) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def require_sha256(value: object, *, field_name: str = "sha256") -> str:
    if not is_sha256(value):
        raise HashingError(
            f"{field_name} must be a 64-character lowercase hexadecimal SHA-256"
        )
    return value


__all__ = [
    "HashingError",
    "canonical_json",
    "canonical_json_bytes",
    "canonical_sha256",
    "canonicalize",
    "is_sha256",
    "require_sha256",
    "sha256_bytes",
    "sha256_file",
]

#!/usr/bin/env python3
"""Shared, dependency-free helpers for the outcome-isolated Myojin workflow."""

from __future__ import annotations

import csv
import hashlib
import os
import re
import tempfile
from pathlib import Path
from typing import Iterable, Mapping, Sequence


RELEASE_ID = "program-context-v2-candidate-2026-08-07"
HLF_MODEL_ID = "ACH-000393"
PROJECT_ROOT = Path(__file__).resolve().parents[4]
CANDIDATE_ROOT = (
    PROJECT_ROOT
    / "Analysis/Multimodal_Program_Projection/candidates"
    / RELEASE_ID
    / "myojin_hlf"
)


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def md5_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.md5()  # noqa: S324 - published source integrity checksum
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def require_within(path: Path, root: Path) -> Path:
    resolved = path.resolve()
    resolved.relative_to(root.resolve())
    return resolved


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_name, path)
    finally:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)


def write_tsv(path: Path, rows: Iterable[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(fields),
                delimiter="\t",
                lineterminator="\n",
                extrasaction="raise",
            )
            writer.writeheader()
            for row in rows:
                writer.writerow({field: row.get(field, "") for field in fields})
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_name, path)
    finally:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


GENE_WITH_ENTREZ = re.compile(r"^(.+?)\s+\((\d+)\)$")


def split_gene_header(value: str) -> tuple[str, str]:
    value = value.strip()
    match = GENE_WITH_ENTREZ.match(value)
    if match:
        return match.group(1).strip(), match.group(2)
    return value, ""


def normalized_model_name(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", value.upper())


def stable_bundle_sha256(paths: Sequence[Path], base: Path) -> tuple[str, list[dict[str, object]]]:
    records: list[dict[str, object]] = []
    for path in sorted(paths, key=lambda p: str(p.resolve())):
        resolved = path.resolve()
        try:
            relative = str(resolved.relative_to(base.resolve()))
        except ValueError:
            relative = str(resolved)
        records.append(
            {
                "relative_path": relative,
                "bytes": resolved.stat().st_size,
                "sha256": sha256_file(resolved),
            }
        )
    canonical = "".join(
        f"{row['relative_path']}\t{row['bytes']}\t{row['sha256']}\n" for row in records
    )
    return sha256_text(canonical), records

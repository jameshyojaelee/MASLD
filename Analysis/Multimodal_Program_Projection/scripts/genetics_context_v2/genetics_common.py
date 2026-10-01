#!/usr/bin/env python3
"""Shared, dependency-free contract helpers for the genetics v2 preflight."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]
# GEN_CONTEXT_ROOT_REL relocates the owned root for a rerun on new inputs.
OWNED_ROOT_REL = Path(os.environ.get(
    "GEN_CONTEXT_ROOT_REL",
    "Analysis/Multimodal_Program_Projection/candidates/"
    "program-context-v2-candidate-2026-08-07/genetics_context",
))
DEFAULT_CANDIDATE_ROOT = PROJECT_ROOT / OWNED_ROOT_REL

FORBIDDEN_ROOTS_REL = (
    Path("GWAS/finemapping/results"),
    Path("RNA-seq/results/causal_inference"),
    Path("RNA-seq/results/multi_evidence"),
    Path("figures/main"),
    Path("RNA-seq/results/manuscript_release"),
)

QUARANTINED_REL = (
    Path("RNA-seq/results/causal_inference/sceqtl/ieqtl_disease_genes.csv"),
    Path("RNA-seq/results/causal_inference/sceqtl/ieqtl_summary.csv"),
)

MISSING_STRINGS = {"", "NA", "NaN", "nan", "None", "null", "NULL", "."}


class ContractError(RuntimeError):
    """Raised when a fail-closed contract condition is not met."""


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def assert_candidate_root(project_root: Path, candidate_root: Path) -> Path:
    project = project_root.resolve()
    expected = (project / OWNED_ROOT_REL).resolve()
    candidate = candidate_root.resolve()
    if candidate != expected:
        raise ContractError(
            f"candidate root must be exactly {expected}; received {candidate}"
        )
    for rel in FORBIDDEN_ROOTS_REL:
        forbidden = (project / rel).resolve()
        if is_relative_to(candidate, forbidden) or is_relative_to(forbidden, candidate):
            raise ContractError(f"candidate root overlaps forbidden path: {forbidden}")
    return candidate


def assert_owned_output(candidate_root: Path, output_path: Path) -> Path:
    candidate = candidate_root.resolve()
    output = output_path.resolve()
    if not is_relative_to(output, candidate):
        raise ContractError(f"output escapes candidate root: {output}")
    return output


def resolve_source(project_root: Path, path_text: str) -> Path:
    candidate = Path(path_text)
    return candidate if candidate.is_absolute() else project_root / candidate


def sha256_file(path: Path, chunk_bytes: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(chunk_bytes)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def clean(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text in MISSING_STRINGS else text


def ensembl_base(value: object) -> str:
    return clean(value).split(".", 1)[0]


def parse_float(value: object) -> float | None:
    text = clean(value)
    if not text:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def delimiter_for(fmt: str) -> str:
    if fmt in {"tsv", "tsv_gz"}:
        return "\t"
    if fmt == "csv":
        return ","
    raise ContractError(f"no tabular delimiter for format {fmt!r}")


def iter_table(path: Path, fmt: str) -> Iterator[dict[str, str]]:
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter_for(fmt))
        if reader.fieldnames is None:
            raise ContractError(f"missing header: {path}")
        for row in reader:
            yield dict(row)


def read_table(path: Path, fmt: str) -> list[dict[str, str]]:
    return list(iter_table(path, fmt))


def inspect_table(path: Path, fmt: str) -> tuple[list[str], int]:
    with _open_text(path) as handle:
        reader = csv.reader(handle, delimiter=delimiter_for(fmt))
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ContractError(f"empty table: {path}") from exc
        n_rows = sum(1 for _ in reader)
    return header, n_rows


def atomic_write_tsv(
    path: Path, rows: Iterable[Mapping[str, object]], fieldnames: Sequence[str]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    materialized = list(rows)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        tmp = Path(handle.name)
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fieldnames),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(materialized)
    os.replace(tmp, path)


def atomic_write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        tmp = Path(handle.name)
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp, path)


def snapshot_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="wb",
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        tmp = Path(handle.name)
    try:
        shutil.copyfile(source, tmp)
        os.replace(tmp, destination)
    finally:
        if tmp.exists():
            tmp.unlink()


def write_stage_seal(
    candidate_root: Path,
    stage: str,
    outputs: Sequence[Path],
    upstream: Sequence[Path] = (),
) -> Path:
    candidate = candidate_root.resolve()
    payload = {
        "stage": stage,
        "created_at_utc": utc_now(),
        "candidate_root": str(candidate),
        "upstream": [
            {
                "path": str(path.resolve()),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for path in upstream
        ],
        "outputs": [
            {
                "relative_path": str(path.resolve().relative_to(candidate)),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for path in outputs
        ],
    }
    seal = candidate / "work" / "stage_seals" / f"{stage}.json"
    assert_owned_output(candidate, seal)
    atomic_write_json(seal, payload)
    return seal


def require_columns(row: Mapping[str, object], columns: Sequence[str], context: str) -> None:
    missing = [column for column in columns if column not in row]
    if missing:
        raise ContractError(f"{context} missing columns: {', '.join(missing)}")


def load_config_tsv(path: Path) -> list[dict[str, str]]:
    return read_table(path, "tsv")


def verify_no_quarantined_inputs(paths: Iterable[Path], project_root: Path) -> None:
    quarantined = {(project_root / rel).resolve() for rel in QUARANTINED_REL}
    offenders = sorted(str(path.resolve()) for path in paths if path.resolve() in quarantined)
    if offenders:
        raise ContractError("quarantined context files are prohibited: " + ", ".join(offenders))


#!/usr/bin/env python3
"""Shared contracts for the isolated ATAC Context v3 candidate."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Iterable, Mapping, Sequence


RELEASE_ID = "atac-context-v3-candidate-2026-08-11-r1"
PROGRAM_RELEASE_ID = "program-context-v2-candidate-2026-08-07"
SEED = 42
MIN_CELLS = 20
MIN_DONORS_PER_GROUP = 4
MIN_PROGRAM_GENES = 8
MIN_PROGRAM_WEIGHT = 0.20
MIN_MAPPED_POSTERIOR = 0.95
ACCESSIBLE_POSTERIOR = 0.50

STANDARD_CHROMS = tuple([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])
COHORTS = ("GSE244832", "GSE281367")
LINEAGES = ("hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk")

LABEL_MAP = {
    "GSE244832": {
        "Hepatocytes": "hepatocyte",
        "Fibroblasts": "stellate",
        "Macrophages": "macrophage",
        "Cholangiocytes": "cholangiocyte",
        "T_cells": "t_nk",
        "Resident_NK": "t_nk",
        "Circulating_NK_NKT": "t_nk",
    },
    "GSE281367": {
        "Hepatocyte": "hepatocyte",
        "Stellate_Cell": "stellate",
        "Macrophage": "macrophage",
        "Kupffer_Cell": "macrophage",
        "Cholangiocyte": "cholangiocyte",
        "NK_T_Cell": "t_nk",
    },
}

EXPECTED_PRIMARY_DONORS = {
    ("GSE244832", "hepatocyte"): (9, 5),
    ("GSE244832", "stellate"): (9, 5),
    ("GSE244832", "macrophage"): (8, 4),
    ("GSE244832", "cholangiocyte"): (5, 1),
    ("GSE244832", "t_nk"): (3, 1),
    ("GSE281367", "hepatocyte"): (6, 6),
    ("GSE281367", "stellate"): (6, 6),
    ("GSE281367", "macrophage"): (6, 6),
    ("GSE281367", "cholangiocyte"): (6, 6),
    ("GSE281367", "t_nk"): (6, 6),
}

EVIDENCE_STATES = {
    "supported",
    "discordant",
    "source_dependent",
    "indeterminate",
    "untestable",
}
GENETIC_STATES = {
    "replicated_accessible",
    "source_dependent",
    "partial",
    "indeterminate",
    "untestable",
}

INPUT_MANIFEST_COLUMNS = (
    "release_id",
    "role",
    "cohort",
    "donor_id",
    "relative_or_absolute_path",
    "bytes",
    "sha256",
)
DONOR_QC_COLUMNS = (
    "release_id",
    "cohort",
    "donor_id",
    "condition",
    "lineage",
    "n_cells",
    "min_cells",
    "contrast_eligible",
    "exclusion_reason",
)
PEAK_MANIFEST_COLUMNS = (
    "release_id",
    "lineage",
    "peak_id",
    "chrom",
    "start0",
    "end",
    "width",
    "supported_gse244832",
    "supported_gse281367",
    "supported_both",
    "blacklist_overlap",
)
GATE_COLUMNS = (
    "release_id",
    "gate",
    "status",
    "reason",
    "manifest_sha256",
    "created_utc",
)


class ContractError(RuntimeError):
    """Raised when a candidate contract would be violated."""


def project_root() -> Path:
    return Path(__file__).resolve().parents[4]


def default_candidate_root() -> Path:
    return (
        project_root()
        / "Analysis/Multimodal_Program_Projection/candidates"
        / RELEASE_ID
    )


def candidate_root(value: str | Path | None = None, *, fixture_mode: bool = False) -> Path:
    root = Path(value or os.environ.get("ATAC_V3_CANDIDATE_ROOT", default_candidate_root()))
    resolved = root.resolve(strict=False)
    if fixture_mode:
        return resolved
    expected_parent = (
        project_root() / "Analysis/Multimodal_Program_Projection/candidates"
    ).resolve()
    if resolved.parent != expected_parent or resolved.name != RELEASE_ID:
        raise ContractError(f"unsafe ATAC v3 candidate root: {resolved}")
    return resolved


def require_new_path(path: Path) -> None:
    if path.exists() or path.is_symlink():
        raise ContractError(f"refusing to overwrite existing path: {path}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def read_tsv(path: Path, required: Sequence[str] = ()) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise ContractError(f"missing required table: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames or []
        missing = [column for column in required if column not in columns]
        if missing:
            raise ContractError(f"{path} missing columns: {missing}")
        return columns, [dict(row) for row in reader]


def write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    require_new_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as handle:
        tmp = Path(handle.name)
        writer = csv.DictWriter(
            handle,
            fieldnames=list(columns),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def write_json(path: Path, payload: Mapping[str, object]) -> None:
    require_new_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(payload, sort_keys=True, indent=2) + "\n"
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        tmp = Path(handle.name)
        handle.write(encoded)
    os.replace(tmp, path)


def parse_bool(value: object, label: str = "value") -> bool:
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    raise ContractError(f"invalid boolean for {label}: {value!r}")


def bh_adjust(pvalues: Sequence[float]) -> list[float]:
    """Benjamini-Hochberg adjustment with stable deterministic tie handling."""
    n = len(pvalues)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: (pvalues[i], i))
    adjusted = [1.0] * n
    running = 1.0
    for rank_index in range(n - 1, -1, -1):
        index = order[rank_index]
        rank = rank_index + 1
        running = min(running, float(pvalues[index]) * n / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def classify_da(q244: float, q281: float, effect244: float, effect281: float) -> str:
    sig244 = q244 < 0.05
    sig281 = q281 < 0.05
    if sig244 and sig281:
        return "supported" if effect244 * effect281 > 0 else "discordant"
    if sig244 or sig281:
        return "source_dependent"
    return "indeterminate"


def classify_genetic(
    mapped_mass: float,
    shared_mass: float,
    cohort244_mass: float,
    cohort281_mass: float,
) -> str:
    if mapped_mass < MIN_MAPPED_POSTERIOR:
        return "untestable"
    if shared_mass >= ACCESSIBLE_POSTERIOR:
        return "replicated_accessible"
    if max(cohort244_mass, cohort281_mass) >= ACCESSIBLE_POSTERIOR:
        return "source_dependent"
    if max(shared_mass, cohort244_mass, cohort281_mass) > 0:
        return "partial"
    return "indeterminate"


def canonical_manifest_hash(rows: Sequence[Mapping[str, object]], columns: Sequence[str]) -> str:
    payload = "".join(
        "\t".join(str(row.get(column, "")) for column in columns) + "\n"
        for row in rows
    )
    return sha256_text(payload)

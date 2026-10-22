#!/usr/bin/env python3
"""Shared fail-closed utilities for the Plan 13 compatible-Visium rerun.

The candidate producer deliberately calls the byte-pinned v1 Moran engine
instead of maintaining a second scientific implementation.  The adapter owns
only registry translation, isolated output, provenance, and release gates.
This gives the v1 regression run the strongest possible semantic check: the
same matched-control, graph, residualization, donor-collapse, null, and
sensitivity functions are exercised with parameterized inputs.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from types import ModuleType
from typing import Iterable, Mapping, Sequence


RELEASE_ID = "program-context-v2-candidate-2026-08-07"
V1_ENGINE_SHA256 = "d198583734ebc02af62973e644963a4d50b87d6de37faa7ef8084b854b94fa94"
V1_REGISTRY_SHA256 = "79fe355dc36a764f67a04171947de9dc1df867b2260a66d6319fc66756ea2e69"
V1_MEMBERSHIP_SHA256 = "900b6e7bd514994fef2dea2c976f5cd0eb715c715a859cb86ced4bdcaf7f5a4c"
V2_REGISTRY_SHA256 = "b134b5f8e46a716e857af22271aa54e1294432b9f9ab735559d7f603cdd29fe7"
V2_MEMBERSHIP_SHA256 = "feec3fc9ccaaffaa9abc1ac5d593cc06845ed8140460c8873a51bc8652d7336b"
V2_READY_SHA256 = "cf109fcb2518470600764cbfb88dde59681e59c0dbb61dc3e0f3d7dfb2a45966"

N_NULL = 9_999
N_SENSITIVITY_NULL = 999
SEED = 42
GRAPH_DISTANCE_MULTIPLIER = 2.5
MIN_GRAPH_SPOTS = 8

CORE_OUTPUTS = (
    "spatial_program_results.tsv",
    "spatial_section_results.tsv",
    "spatial_matching_audit.tsv",
    "spatial_graph_audit.tsv",
    "spatial_null_summary.tsv",
)


@dataclass(frozen=True)
class Paths:
    project_root: Path
    candidate_root: Path
    native_root: Path
    hotspot_root: Path
    v1_engine: Path
    v1_registry: Path
    v1_membership: Path
    v1_results: Path
    gene_metadata: Path
    datasets: Mapping[str, Path]


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[4]


def build_paths(project_root: Path | None = None) -> Paths:
    root = (project_root or project_root_from_script()).resolve()
    candidate = (
        root
        / "Analysis/Multimodal_Program_Projection/candidates"
        / RELEASE_ID
        / "spatial_context"
    )
    mm = root / "Analysis/Multimodal_Program_Projection"
    return Paths(
        project_root=root,
        candidate_root=candidate,
        native_root=candidate / "native_spatial",
        hotspot_root=mm / "candidates" / RELEASE_ID / "hotspot",
        v1_engine=mm / "scripts/04_spatial_projection.py",
        v1_registry=mm / "results/frozen_programs.tsv",
        v1_membership=mm / "results/frozen_program_membership.tsv",
        v1_results=mm / "results/spatial",
        gene_metadata=root / "data/gencode_v49_gene_metadata.tsv.gz",
        datasets={
            "GSE192741": root
            / "Analysis/Spatial/results/cell2location/spatial_model/spatial_deconvolved.h5ad",
            "Vu_et_al_2025": root
            / "Analysis/Spatial/results/cell2location/spatial_model_vu/spatial_deconvolved_vu.h5ad",
        },
    )


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def aggregate_file_hash(paths: Iterable[Path], base: Path) -> str:
    records = []
    for path in sorted(paths, key=lambda item: item.relative_to(base).as_posix()):
        records.append(
            f"{path.relative_to(base).as_posix()}\0{path.stat().st_size}\0{sha256_file(path)}\n"
        )
    return hashlib.sha256("".join(records).encode("utf-8")).hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def canonical_rows_sha256(rows: Sequence[Mapping[str, object]], columns: Sequence[str]) -> str:
    payload = "\n".join(
        "\t".join(str(row.get(column, "")) for column in columns) for row in rows
    ) + "\n"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def require_file_hash(path: Path, expected: str, label: str) -> None:
    if not path.is_file():
        raise RuntimeError(f"missing {label}: {path}")
    observed = sha256_file(path)
    if observed != expected:
        raise RuntimeError(f"{label} hash drift: expected {expected}, observed {observed}")


def verify_v1_anchors(paths: Paths) -> None:
    require_file_hash(paths.v1_engine, V1_ENGINE_SHA256, "v1 spatial producer")
    require_file_hash(paths.v1_registry, V1_REGISTRY_SHA256, "v1 registry")
    require_file_hash(paths.v1_membership, V1_MEMBERSHIP_SHA256, "v1 membership")
    for name in CORE_OUTPUTS:
        if not (paths.v1_results / name).is_file():
            raise RuntimeError(f"missing canonical v1 regression target: {paths.v1_results / name}")
    summary = read_tsv(paths.candidate_root / "v1_preservation_summary.tsv")
    if len(summary) != 1 or summary[0].get("status") != "pass_byte_identical":
        raise RuntimeError("SP-INT-01 v1 byte-preservation gate is not passing")


def verify_hotspot_ready(paths: Paths) -> dict[str, str]:
    registry = paths.hotspot_root / "program_registry_v2.tsv"
    membership = paths.hotspot_root / "program_membership_v2.tsv"
    ready_path = paths.hotspot_root / "READY"
    require_file_hash(registry, V2_REGISTRY_SHA256, "Plan 20 v2 registry")
    require_file_hash(membership, V2_MEMBERSHIP_SHA256, "Plan 20 v2 membership")
    require_file_hash(ready_path, V2_READY_SHA256, "Plan 20 READY")
    rows = read_tsv(ready_path)
    if len(rows) != 1:
        raise RuntimeError("Plan 20 READY must contain exactly one row")
    row = rows[0]
    required = {
        "release_id": RELEASE_ID,
        "status": "ready_for_external_testing",
        "registry_sha256": V2_REGISTRY_SHA256,
        "membership_table_sha256": V2_MEMBERSHIP_SHA256,
        "external_outcomes_read": "FALSE",
    }
    for key, expected in required.items():
        if row.get(key) != expected:
            raise RuntimeError(f"Plan 20 READY {key} mismatch: expected {expected!r}, got {row.get(key)!r}")
    transitive = {
        "validation_status_sha256": paths.hotspot_root / "validation_status.tsv",
        "release_manifest_sha256": paths.hotspot_root / "release_manifest.tsv",
    }
    for key, path in transitive.items():
        require_file_hash(path, row[key], f"Plan 20 {path.name}")
    return row


def v1_universe(paths: Paths) -> list[dict[str, object]]:
    registry = read_tsv(paths.v1_registry)
    if len(registry) != 22:
        raise RuntimeError(f"expected 22 immutable v1 programs, found {len(registry)}")
    return [
        {
            "program_id": row["program_id"],
            "legacy_program_id": row["program_id"],
            "cell_type": row["cell_type"],
            "module": row["module"],
            "program_name": row["program_name"],
            "display_order": row["display_order"],
            "membership_sha256": "legacy_table_level_only",
        }
        for row in registry
    ]


def v2_universe(paths: Paths) -> list[dict[str, object]]:
    verify_hotspot_ready(paths)
    rows = [row for row in read_tsv(paths.hotspot_root / "program_registry_v2.tsv") if row["robust_display"] == "TRUE"]
    rows.sort(key=lambda row: (row["cell_type"], int(row["module"]), row["program_uid"]))
    if not rows:
        raise RuntimeError("Plan 20 robust_display universe is empty")
    allowed_lineages = {"hepatocytes", "fibroblasts", "macrophages", "cholangiocytes"}
    unsupported = sorted({row["cell_type"] for row in rows}.difference(allowed_lineages))
    if unsupported:
        raise RuntimeError(f"v2 Visium universe contains unsupported lineages: {unsupported}")
    out = []
    for order, row in enumerate(rows, start=1):
        out.append(
            {
                "program_id": row["program_uid"],
                "legacy_program_id": f"{row['cell_type']}::{row['module']}",
                "cell_type": row["cell_type"],
                "module": row["module"],
                "program_name": row["module_name"],
                "display_order": order,
                "membership_sha256": row["membership_sha256"],
                "primary_direction": row["primary_direction"],
                "external_test_eligible": row["external_test_eligible"],
            }
        )
    return out


def universe_hash(rows: Sequence[Mapping[str, object]]) -> str:
    columns = (
        "program_id",
        "legacy_program_id",
        "cell_type",
        "module",
        "program_name",
        "display_order",
        "membership_sha256",
    )
    return canonical_rows_sha256(rows, columns)


def load_legacy_engine(paths: Paths) -> ModuleType:
    """Import the pinned v1 engine without allowing a canonical write target.

    The v1 module creates its configured output directory at import time.  We
    require that directory to pre-exist and assert that importing does not
    change its metadata.  Candidate code never calls the v1 ``main`` function.
    """
    verify_v1_anchors(paths)
    canonical_out = paths.v1_results.resolve()
    if not canonical_out.is_dir():
        raise RuntimeError(f"canonical v1 output directory is missing: {canonical_out}")
    before = canonical_out.stat()
    old_smoke = os.environ.get("FIG4_SPATIAL_SMOKE")
    old_project_root = os.environ.get("MASLD_PROJECT_ROOT")
    os.environ["FIG4_SPATIAL_SMOKE"] = "FALSE"
    os.environ["MASLD_PROJECT_ROOT"] = str(paths.project_root)
    try:
        spec = importlib.util.spec_from_file_location("_pinned_v1_spatial_engine", paths.v1_engine)
        if spec is None or spec.loader is None:
            raise RuntimeError("could not construct import specification for the v1 engine")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        if old_smoke is None:
            os.environ.pop("FIG4_SPATIAL_SMOKE", None)
        else:
            os.environ["FIG4_SPATIAL_SMOKE"] = old_smoke
        if old_project_root is None:
            os.environ.pop("MASLD_PROJECT_ROOT", None)
        else:
            os.environ["MASLD_PROJECT_ROOT"] = old_project_root
    after = canonical_out.stat()
    if (before.st_mtime_ns, before.st_size, before.st_ino) != (after.st_mtime_ns, after.st_size, after.st_ino):
        raise RuntimeError("importing the pinned engine changed canonical v1 output-directory metadata")
    if Path(module.OUT).resolve() != canonical_out:
        raise RuntimeError(f"pinned v1 engine resolved an unexpected output root: {module.OUT}")
    if Path(module.BASE).resolve() != paths.project_root:
        raise RuntimeError(f"pinned v1 engine resolved an unexpected project root: {module.BASE}")
    module.SMOKE = False
    module.N_NULL = N_NULL
    module.N_SENSITIVITY_NULL = N_SENSITIVITY_NULL
    module.SEED = SEED
    module.GRAPH_DISTANCE_MULTIPLIER = GRAPH_DISTANCE_MULTIPLIER
    module.MIN_GRAPH_SPOTS = MIN_GRAPH_SPOTS
    return module


def selected_v2_membership_rows(paths: Paths) -> list[dict[str, str]]:
    selected = {str(row["program_id"]): row for row in v2_universe(paths)}
    rows = []
    seen = set()
    for row in read_tsv(paths.hotspot_root / "program_membership_v2.tsv"):
        uid = row["program_uid"]
        if uid not in selected:
            continue
        if row["membership_sha256"] != selected[uid]["membership_sha256"]:
            raise RuntimeError(f"per-row membership hash mismatch for {uid}")
        if row["mapped_symbol_status"] != "gencode_v49_unique_symbol_confirmed":
            continue
        symbol = row["mapped_symbol"].strip()
        if not symbol:
            raise RuntimeError(f"confirmed membership row lacks mapped symbol for {uid}")
        key = (uid, row["source_gene"], symbol)
        if key in seen:
            raise RuntimeError(f"duplicate selected membership row: {key}")
        seen.add(key)
        rows.append(row)
    missing = sorted(set(selected).difference({row["program_uid"] for row in rows}))
    if missing:
        raise RuntimeError(f"selected programs have no confirmed mapped genes: {missing}")
    return rows


def prepare_engine_inputs(paths: Paths, version: str):
    """Return pandas registry/membership frames in the exact v1 engine schema."""
    import pandas as pd

    if version == "v1":
        registry = pd.read_csv(paths.v1_registry, sep="\t")
        membership = pd.read_csv(paths.v1_membership, sep="\t")
        membership = membership[
            membership["mapped_symbol"].fillna(False) & membership["gene_symbol"].notna()
        ].copy()
        return registry, membership, v1_universe(paths)
    if version != "v2":
        raise ValueError(f"unknown registry version: {version}")
    universe = v2_universe(paths)
    registry = pd.DataFrame(
        [
            {
                "program_id": row["program_id"],
                "display_order": row["display_order"],
                "cell_type": row["cell_type"],
                "module": int(row["module"]),
                "program_name": row["program_name"],
            }
            for row in universe
        ]
    )
    membership = pd.DataFrame(
        [
            {
                "program_id": row["program_uid"],
                "cell_type": row["cell_type"],
                "module": int(row["module"]),
                "source_gene": row["source_gene"],
                "gene_symbol": row["mapped_symbol"],
                "mapped_symbol": True,
                "original_l1_weight": float(row["original_l1_weight"]),
            }
            for row in selected_v2_membership_rows(paths)
        ]
    )
    return registry, membership, universe


def score_variants(z, weights):
    """Compute primary/equal/leave-highest-weight scores exactly as v1."""
    import numpy as np

    weights = np.asarray(weights, dtype=float)
    if z.ndim != 2 or z.shape[1] != len(weights) or len(weights) < 2:
        raise ValueError("score fixture requires a 2D matrix and at least two aligned weights")
    weights = weights / weights.sum()
    primary = z @ weights
    equal = np.mean(z, axis=1)
    top = int(np.argmax(weights))
    leave_weights = np.delete(weights, top)
    leave_weights = leave_weights / leave_weights.sum()
    leave = np.delete(z, top, axis=1) @ leave_weights
    return primary, equal, leave, top, leave_weights


def json_dump(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")

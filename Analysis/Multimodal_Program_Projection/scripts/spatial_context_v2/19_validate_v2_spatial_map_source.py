#!/usr/bin/env python3
"""Independently validate the candidate v2 per-spot spatial map source."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import math
import sys
from collections import Counter
from pathlib import Path

from visium_rerun_lib import (
    RELEASE_ID,
    build_paths,
    read_tsv,
    sha256_file,
    utc_now,
    verify_hotspot_ready,
    verify_v1_anchors,
    write_tsv,
)


class Checks:
    def __init__(self) -> None:
        self.rows: list[dict[str, str]] = []

    def require(self, condition: bool, check_id: str, detail: str) -> None:
        self.rows.append(
            {"check_id": check_id, "status": "PASS" if condition else "FAIL", "detail": detail}
        )

    @property
    def passed(self) -> bool:
        return bool(self.rows) and all(row["status"] == "PASS" for row in self.rows)


def read_gzip_tsv(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError("map source lacks header")
        return list(reader)


def independent_selection(samples: list[dict[str, str]]) -> dict[str, str]:
    selected = {}
    for dataset in sorted({row["dataset"] for row in samples}):
        part = [row for row in samples if row["dataset"] == dataset]
        counts = sorted(int(row["n_spots"]) for row in part)
        n = len(counts)
        median = (
            float(counts[n // 2])
            if n % 2
            else (counts[n // 2 - 1] + counts[n // 2]) / 2.0
        )
        chosen = min(
            part,
            key=lambda row: (abs(int(row["n_spots"]) - median), row["technical_id"]),
        )
        selected[dataset] = chosen["technical_id"]
    return selected


def canonical_row_hash(rows: list[dict[str, str]], fields: list[str]) -> str:
    payload = "".join(
        "\t".join(row[field] for field in fields) + "\n"
        for row in sorted(
            rows,
            key=lambda row: (row["dataset"], row["program_uid"], row["spot_id"]),
        )
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def validate(project_root: Path | None) -> tuple[Path, Checks, list[dict[str, str]]]:
    paths = build_paths(project_root)
    verify_v1_anchors(paths)
    verify_hotspot_ready(paths)
    root = paths.candidate_root / "map_source_v2"
    if not root.is_dir():
        raise RuntimeError(f"missing map source: {root}")
    for name in ("validation_report.tsv", "release_manifest.tsv", "MAP_SOURCE_READY"):
        if (root / name).exists():
            raise RuntimeError(f"refusing to overwrite map validation artifact: {root / name}")
    checks = Checks()

    build = read_tsv(root / "BUILD_COMPLETE")
    map_path = root / "per_spot_program_map.tsv.gz"
    checks.require(
        len(build) == 1
        and build[0].get("status") == "built_pending_independent_validation"
        and build[0].get("map_sha256") == sha256_file(map_path),
        "build_seal",
        "BUILD_COMPLETE links the immutable gzip source",
    )
    for row in read_tsv(root / "source_manifest.tsv"):
        source = paths.project_root / row["relative_path"]
        ok = (
            source.is_file()
            and source.stat().st_size == int(row["bytes"])
            and sha256_file(source) == row["sha256"]
        )
        checks.require(ok, f"source::{row['source_role']}", "size and SHA256 rederived")

    native_samples = read_tsv(paths.native_root / "v2_candidate/native_sample_manifest.tsv")
    expected_selection = independent_selection(native_samples)
    selection_rows = read_tsv(root / "map_selection.tsv")
    observed_selection = {row["dataset"]: row["reporting_unit_id"] for row in selection_rows}
    checks.require(
        observed_selection == expected_selection,
        "outcome_independent_selection_rederived",
        f"observed={observed_selection};expected={expected_selection}",
    )
    checks.require(
        len(selection_rows) == 2
        and all(row["program_outcomes_used_for_selection"] == "FALSE" for row in selection_rows)
        and all(row["same_unit_for_all_programs"] == "TRUE" for row in selection_rows),
        "selection_scope",
        "two datasets; same unit for both programs; no program outcome used",
    )
    resolution = {row["dataset"]: row["biological_unit_resolution"] for row in selection_rows}
    checks.require(
        resolution
        == {
            "GSE192741": "resolved_human_donor",
            "Vu_et_al_2025": "unresolved_physical_array",
        },
        "biological_unit_boundary",
        str(resolution),
    )

    rows = read_gzip_tsv(map_path)
    checks.require(
        len(rows) == int(build[0]["row_count"]) and len(rows) > 0,
        "map_row_count",
        f"rows={len(rows)}",
    )
    keys = [(row["dataset"], row["program_uid"], row["spot_id"]) for row in rows]
    checks.require(len(keys) == len(set(keys)), "unique_map_key", "dataset/program/spot unique")
    programs = sorted({row["program_uid"] for row in rows})
    datasets = sorted({row["dataset"] for row in rows})
    checks.require(len(programs) == 2 and len(datasets) == 2, "complete_map_family", f"2 programs x 2 datasets; {len(rows)} rows")
    checks.require(
        all(row["reporting_unit_id"] == expected_selection[row["dataset"]] for row in rows),
        "single_selected_unit_per_dataset",
        "all rows resolve to the independently selected unit",
    )
    checks.require(
        all(row["graph_eligible"] == "TRUE" for row in rows),
        "graph_eligible_only",
        "small disconnected tissue islands excluded",
    )

    numeric_ok = True
    max_abs_z = 0.0
    for row in rows:
        values = [float(row[field]) for field in ("x", "y", "residual_program_score_z")]
        numeric_ok &= all(math.isfinite(value) for value in values)
        max_abs_z = max(max_abs_z, abs(values[2]))
    checks.require(numeric_ok and max_abs_z < 20, "finite_map_values", f"max_abs_z={max_abs_z:.4g}")

    registry = {
        row["program_uid"]: row
        for row in read_tsv(paths.hotspot_root / "program_registry_v2.tsv")
        if row["robust_display"] == "TRUE"
    }
    membership_ok = all(
        row["program_uid"] in registry
        and row["membership_sha256"] == registry[row["program_uid"]]["membership_sha256"]
        for row in rows
    )
    checks.require(membership_ok, "membership_identity", "all rows link to frozen v2 membership")

    audit = read_tsv(root / "map_scoring_audit.tsv")
    audit_by_key = {(row["dataset"], row["program_uid"]): row for row in audit}
    native = {
        (row["dataset"], row["program_id"]): row
        for row in read_tsv(paths.native_root / "v2_candidate/spatial_program_results.tsv")
    }
    counts = Counter((row["dataset"], row["program_uid"]) for row in rows)
    audit_ok = len(audit_by_key) == 4
    for key, row in audit_by_key.items():
        source = native.get(key)
        audit_ok &= source is not None
        if source is None:
            continue
        audit_ok &= int(row["n_graph_eligible_spots"]) == counts[key]
        audit_ok &= int(row["n_genes_measured"]) == int(source["n_measured"])
        audit_ok &= math.isclose(
            float(row["retained_l1_weight"]),
            float(source["retained_l1_weight"]),
            rel_tol=1e-12,
            abs_tol=1e-14,
        )
        audit_ok &= row["inferential_use"] == "FALSE_illustrative_map_only"
    checks.require(audit_ok, "native_scoring_link", "coverage/L1/spot counts rederive")

    execution = {row["parameter"]: row["value"] for row in read_tsv(root / "execution_manifest.tsv")}
    checks.require(
        execution.get("program_outcomes_used_for_selection") == "FALSE"
        and execution.get("program_count") == "2"
        and execution.get("map_row_count") == str(len(rows)),
        "execution_contract",
        "selection, family, and row count linked",
    )
    return root, checks, rows


def write_outputs(root: Path, checks: Checks, rows: list[dict[str, str]]) -> None:
    report = root / "validation_report.tsv"
    write_tsv(report, ("check_id", "status", "detail"), checks.rows)
    if not checks.passed:
        raise RuntimeError(f"map-source validation failed: {report}")
    manifest = root / "release_manifest.tsv"
    manifest_rows = []
    for path in sorted(root.iterdir()):
        if not path.is_file() or path.name in {manifest.name, "MAP_SOURCE_READY"}:
            continue
        manifest_rows.append(
            {
                "relative_path": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    write_tsv(manifest, ("relative_path", "bytes", "sha256"), manifest_rows)
    ready = root / "MAP_SOURCE_READY"
    write_tsv(
        ready,
        (
            "release_id",
            "status",
            "row_count",
            "map_source_sha256",
            "validation_report_sha256",
            "release_manifest_sha256",
            "canonical_row_sha256",
            "validated_utc",
        ),
        [
            {
                "release_id": RELEASE_ID,
                "status": "ready_candidate_illustrative_maps",
                "row_count": len(rows),
                "map_source_sha256": sha256_file(root / "per_spot_program_map.tsv.gz"),
                "validation_report_sha256": sha256_file(report),
                "release_manifest_sha256": sha256_file(manifest),
                "canonical_row_sha256": canonical_row_hash(rows, list(rows[0])),
                "validated_utc": utc_now(),
            }
        ],
    )
    print(f"validated v2 spatial map source: checks={len(checks.rows)} rows={len(rows)}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=None)
    args = parser.parse_args()
    try:
        root, checks, rows = validate(args.project_root)
        write_outputs(root, checks, rows)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

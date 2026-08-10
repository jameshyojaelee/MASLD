#!/usr/bin/env python3
"""Freeze/check the validation-hotfix v2 Visium execution addendum."""

from __future__ import annotations

import argparse
import hashlib
import sys
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


FREEZE_FILES = ("semantic_runner_specification.tsv", "source_manifest.tsv")


def aggregate_manifest(path: Path, base: Path) -> str:
    rows = read_tsv(path)
    payload = []
    for row in rows:
        source = base / row["relative_path"]
        if (
            not source.is_file()
            or source.stat().st_size != int(row["bytes"])
            or sha256_file(source) != row["sha256"]
        ):
            raise RuntimeError(f"base freeze manifest drift: {source}")
        payload.append(
            f"{row['relative_path']}\0{row['bytes']}\0{row['sha256']}\n"
        )
    return hashlib.sha256("".join(payload).encode()).hexdigest()


def resolve_ready_report(v1_root: Path) -> tuple[Path, dict[str, str]]:
    ready = v1_root / "READY"
    rows = read_tsv(ready)
    if len(rows) != 1 or rows[0].get("status") != "pass_v1_regression":
        raise RuntimeError("semantic v1 READY is absent or invalid")
    row = rows[0]
    report = (v1_root / row.get("validation_report_path", "validation_report.tsv")).resolve()
    if report.parent != v1_root.resolve() or not report.is_file():
        raise RuntimeError("semantic v1 READY report path is unsafe or missing")
    if row.get("validation_report_sha256") != sha256_file(report):
        raise RuntimeError("semantic v1 READY report hash mismatch")
    return report, row


def write(paths) -> None:
    verify_v1_anchors(paths)
    verify_hotspot_ready(paths)
    if (paths.native_root / "v2_candidate").exists():
        raise RuntimeError("semantic runner must be frozen before v2 candidate outcomes exist")
    root = paths.native_root / "semantic_runner_freeze"
    owned = [root / name for name in FREEZE_FILES] + [root / "freeze_manifest.tsv"]
    if any(path.exists() for path in owned):
        raise RuntimeError("refusing to overwrite semantic-runner freeze")
    root.mkdir(parents=True, exist_ok=True)

    base_manifest = paths.native_root / "freeze/freeze_manifest.tsv"
    base_aggregate = aggregate_manifest(base_manifest, paths.candidate_root)
    v1_root = paths.native_root / "v1_regression"
    ready_report, ready_row = resolve_ready_report(v1_root)
    script_root = Path(__file__).resolve().parent

    sources = [
        (base_manifest, "outcome_blind_native_freeze_manifest"),
        (v1_root / "READY", "v1_semantic_ready"),
        (ready_report, "v1_semantic_validation_report"),
        (v1_root / "validation_report.tsv", "v1_original_strict_failed_report"),
        (v1_root / "validator_hotfix_manifest.tsv", "v1_validator_hotfix_manifest"),
        (script_root / "04b_freeze_visium_semantic_runner.py", "semantic_runner_freezer"),
        (script_root / "05b_run_visium_v2_semantic_ready.py", "semantic_ready_producer"),
        (script_root / "06c_validate_visium_v2_semantic.py", "semantic_v2_validator"),
        (script_root / "run_visium_v2_semantic.sbatch", "semantic_v2_wrapper"),
        (script_root / "05_run_visium_rerun.py", "original_frozen_producer"),
        (script_root / "06_validate_visium_rerun.py", "original_frozen_validator"),
    ]
    source_rows = []
    for path, role in sources:
        if not path.is_file():
            raise RuntimeError(f"missing semantic-runner source {role}: {path}")
        source_rows.append(
            {
                "relative_path": path.relative_to(paths.project_root).as_posix(),
                "source_role": role,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    write_tsv(
        root / "source_manifest.tsv",
        ("relative_path", "source_role", "bytes", "sha256"),
        source_rows,
    )
    write_tsv(
        root / "semantic_runner_specification.tsv",
        ("parameter", "value"),
        [
            {"parameter": "release_id", "value": RELEASE_ID},
            {"parameter": "task_id", "value": "SP-INT-03-semantic-runner"},
            {
                "parameter": "scientific_engine_change",
                "value": "FALSE_original_byte_pinned_engine_and_model_unchanged",
            },
            {
                "parameter": "validation_change",
                "value": "READY_report_resolution_and_finite_null_semantic_bound_only",
            },
            {
                "parameter": "minimum_finite_null_fraction",
                "value": "0.95",
            },
            {
                "parameter": "v1_ready_sha256",
                "value": sha256_file(v1_root / "READY"),
            },
            {
                "parameter": "v1_validation_report_path",
                "value": ready_report.name,
            },
            {
                "parameter": "v1_validation_report_sha256",
                "value": sha256_file(ready_report),
            },
            {
                "parameter": "base_freeze_manifest_aggregate_sha256",
                "value": base_aggregate,
            },
            {"parameter": "v2_outcomes_read", "value": "FALSE"},
            {"parameter": "frozen_at_utc", "value": utc_now()},
        ],
    )
    manifest_rows = []
    for name in FREEZE_FILES:
        path = root / name
        manifest_rows.append(
            {
                "relative_path": path.relative_to(paths.candidate_root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    write_tsv(
        root / "freeze_manifest.tsv",
        ("relative_path", "bytes", "sha256"),
        manifest_rows,
    )
    print(
        "froze SP-INT-03 semantic runner before v2 outcomes: "
        f"v1_READY={sha256_file(v1_root / 'READY')}"
    )


def check(paths) -> None:
    verify_v1_anchors(paths)
    verify_hotspot_ready(paths)
    root = paths.native_root / "semantic_runner_freeze"
    manifest = root / "freeze_manifest.tsv"
    rows = read_tsv(manifest)
    if not rows:
        raise RuntimeError("semantic-runner freeze manifest is empty")
    for row in rows:
        path = paths.candidate_root / row["relative_path"]
        if (
            not path.is_file()
            or path.stat().st_size != int(row["bytes"])
            or sha256_file(path) != row["sha256"]
        ):
            raise RuntimeError(f"semantic-runner freeze artifact drift: {path}")
    for row in read_tsv(root / "source_manifest.tsv"):
        path = paths.project_root / row["relative_path"]
        if (
            not path.is_file()
            or path.stat().st_size != int(row["bytes"])
            or sha256_file(path) != row["sha256"]
        ):
            raise RuntimeError(f"semantic-runner source drift: {path}")
    spec = {
        row["parameter"]: row["value"]
        for row in read_tsv(root / "semantic_runner_specification.tsv")
    }
    report, _ = resolve_ready_report(paths.native_root / "v1_regression")
    if (
        spec.get("scientific_engine_change")
        != "FALSE_original_byte_pinned_engine_and_model_unchanged"
        or spec.get("v2_outcomes_read") != "FALSE"
        or spec.get("v1_ready_sha256")
        != sha256_file(paths.native_root / "v1_regression/READY")
        or spec.get("v1_validation_report_sha256") != sha256_file(report)
    ):
        raise RuntimeError("semantic-runner specification drift")
    print("SP-INT-03 semantic-runner freeze check passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("write", "check"))
    parser.add_argument("--project-root", type=Path, default=None)
    args = parser.parse_args()
    try:
        paths = build_paths(args.project_root)
        write(paths) if args.mode == "write" else check(paths)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

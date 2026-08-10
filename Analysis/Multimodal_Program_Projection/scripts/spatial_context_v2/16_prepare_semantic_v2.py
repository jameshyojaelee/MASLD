#!/usr/bin/env python3
"""Prepare an isolated evidence-state-v2 rebuild without mutating historical bundles."""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

from contract_lib import (
    ContractError,
    RELEASE_ID,
    SEMANTIC_CONTRACT_ID,
    read_tsv,
    sha256_file,
    write_tsv,
)


BASELINE_COLUMNS = ("source_scope", "relative_path", "bytes", "sha256")
COPY_COLUMNS = ("source_relative_path", "target_relative_path", "bytes", "sha256", "copy_role")
PREREQUISITE_FILES = (
    "map_selection_manifest.tsv",
    "v1_preservation_anchor.tsv",
    "v1_preservation_baseline.tsv",
    "v1_preservation_comparison.tsv",
    "v1_preservation_final.tsv",
    "v1_preservation_scope.tsv",
    "v1_preservation_summary.tsv",
)
PREREQUISITE_TREES = ("protein_atac_native", "native_spatial")
EXPECTED_NATIVE_V2_READY_SHA256 = (
    "fcff53888cd8a29adc817752be4f3798e7fb3bcd7903f02622a90af729089898"
)


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[4]


def _relative_files(root: Path) -> list[Path]:
    return sorted(path for path in root.rglob("*") if path.is_file())


def _baseline_rows(root: Path, scope: str) -> list[dict[str, object]]:
    return [
        {
            "source_scope": scope,
            "relative_path": path.relative_to(root).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in _relative_files(root)
    ]


def _copy_file(source: Path, target: Path) -> None:
    if not source.is_file():
        raise ContractError(f"missing semantic-v2 prerequisite: {source}")
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    if source.stat().st_size != target.stat().st_size or sha256_file(source) != sha256_file(target):
        raise ContractError(f"semantic-v2 prerequisite copy drift: {source}")


def prepare(source_root: Path, target_root: Path, yak_native_root: Path) -> None:
    if target_root.exists():
        raise ContractError(f"refusing to overwrite semantic-v2 candidate root: {target_root}")
    if not source_root.is_dir() or not yak_native_root.is_dir():
        raise ContractError("historical Plan 13 and Yakubovsky roots must exist")
    native_ready = source_root / "native_spatial/v2_candidate/READY"
    _, ready_rows = read_tsv(
        native_ready,
        ("registry_version", "status", "validation_report_sha256"),
    )
    if len(ready_rows) != 1 or ready_rows[0]["status"] != "pass_v2_candidate":
        raise ContractError("provenance-complete native spatial v2 READY is not available")
    if sha256_file(native_ready) != EXPECTED_NATIVE_V2_READY_SHA256:
        raise ContractError("native spatial v2 READY does not match the reviewed semantic-runner seal")
    report_relative = ready_rows[0].get("validation_report_path", "validation_report.tsv")
    report = (native_ready.parent / report_relative).resolve()
    if (
        report.parent != native_ready.parent.resolve()
        or not report.is_file()
        or sha256_file(report) != ready_rows[0]["validation_report_sha256"]
    ):
        raise ContractError("native spatial v2 READY-linked validation report drift")

    target_root.mkdir(parents=True)
    copy_rows: list[dict[str, object]] = []
    for relative in PREREQUISITE_FILES:
        source = source_root / relative
        target = target_root / relative
        _copy_file(source, target)
        copy_rows.append(
            {
                "source_relative_path": source.relative_to(source_root.parent).as_posix(),
                "target_relative_path": target.relative_to(target_root).as_posix(),
                "bytes": target.stat().st_size,
                "sha256": sha256_file(target),
                "copy_role": "outcome_blind_historical_prerequisite",
            }
        )
    for relative in PREREQUISITE_TREES:
        source_tree = source_root / relative
        if not source_tree.is_dir():
            raise ContractError(f"missing semantic-v2 prerequisite tree: {source_tree}")
        for source in _relative_files(source_tree):
            target = target_root / relative / source.relative_to(source_tree)
            _copy_file(source, target)
            copy_rows.append(
                {
                    "source_relative_path": source.relative_to(source_root.parent).as_posix(),
                    "target_relative_path": target.relative_to(target_root).as_posix(),
                    "bytes": target.stat().st_size,
                    "sha256": sha256_file(target),
                    "copy_role": "sealed_native_prerequisite",
                }
            )

    baseline_rows = _baseline_rows(source_root, "historical_plan13_root")
    baseline_rows.extend(_baseline_rows(yak_native_root, "historical_yakubovsky_root"))
    write_tsv(target_root / "semantic_v2_source_baseline.tsv", BASELINE_COLUMNS, baseline_rows)
    write_tsv(target_root / "semantic_v2_prerequisite_copies.tsv", COPY_COLUMNS, copy_rows)
    write_tsv(
        target_root / "semantic_v2_preparation_status.tsv",
        (
            "release_id",
            "semantic_contract_id",
            "status",
            "source_root",
            "yak_native_root",
            "n_baseline_files",
            "n_copied_files",
            "canonical_write_authorized",
        ),
        [
            {
                "release_id": RELEASE_ID,
                "semantic_contract_id": SEMANTIC_CONTRACT_ID,
                "status": "isolated_semantic_v2_prerequisites_prepared",
                "source_root": str(source_root),
                "yak_native_root": str(yak_native_root),
                "n_baseline_files": len(baseline_rows),
                "n_copied_files": len(copy_rows),
                "canonical_write_authorized": False,
            }
        ],
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=project_root_from_script())
    parser.add_argument("--source-root", type=Path, default=None)
    parser.add_argument("--target-root", type=Path, default=None)
    parser.add_argument("--yak-native-root", type=Path, default=None)
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    release_root = (
        project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID
    )
    source_root = (args.source_root or (release_root / "spatial_context")).resolve()
    target_root = (
        args.target_root or (release_root / "spatial_context_semantic_v2_2026-08-08")
    ).resolve()
    yak_native_root = (
        args.yak_native_root
        or (project_root / "Analysis/Spatial/candidates" / RELEASE_ID / "yakubovsky2026")
    ).resolve()
    try:
        prepare(source_root, target_root, yak_native_root)
        print(f"PASS: prepared isolated semantic-v2 root {target_root}")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

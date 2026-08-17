#!/usr/bin/env python3
"""Seal the fixture-tested genetics implementation without biological outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path


RELEASE_ID = "atac-context-v3-candidate-2026-08-11-r1"
ROOT = Path(__file__).resolve().parents[4]
PACKAGE = Path(__file__).resolve().parent
CANDIDATE = (
    ROOT / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID
).resolve()
SOURCES = (
    "07_prepare_genetic_replay.py",
    "08_liftover_variants.R",
    "09_aggregate_genetic_context.py",
    "11_validate_genetics.py",
    "13_render_genetic.R",
    "20_validate_genetics_implementation.py",
    "21_verify_replay_batch_inputs.py",
    "22_validate_replay_batch.py",
    "genetics_context.py",
    "genetics_export_helpers.R",
    "README.md",
    "schemas/output_schemas.json",
    "slurm/08_genetics_batch.sbatch",
    "slurm/09_genetics_validate.sbatch",
    "slurm/14_validate_genetics_implementation.sbatch",
    "tests/smoke_coloc_susie_export.R",
    "tests/test_contracts.py",
    "tests/test_post_gate_contracts.py",
)


class ValidationError(RuntimeError):
    pass


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pending-root", type=Path, required=True)
    parser.add_argument("--non-genetic-ready-sha256", required=True)
    parser.add_argument("--candidate-audit-ready-sha256", required=True)
    parser.add_argument("--integration-ready-sha256", required=True)
    parser.add_argument("--fig4f-ready-sha256", required=True)
    parser.add_argument("--runner-before-sha256", required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise ValidationError(f"missing or unsafe table: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def resolve_artifact(value: str) -> Path:
    candidate_path = CANDIDATE / value
    return candidate_path if candidate_path.is_file() else ROOT / value


def verify_ready_gate(path: Path, gate: str) -> None:
    rows = read(path)
    if not rows or any(
        row.get("release_id") != RELEASE_ID
        or row.get("gate") != gate
        or row.get("status") != "READY"
        for row in rows
    ):
        raise ValidationError(f"invalid {gate} seal")
    for row in rows:
        artifact = resolve_artifact(row["artifact"])
        if not artifact.is_file() or sha256(artifact) != row["sha256"]:
            raise ValidationError(f"{gate} artifact hash mismatch: {row['artifact']}")


def verify_candidate_audit() -> None:
    rows = read(CANDIDATE / "CANDIDATE_AUDIT_READY")
    if not rows or any(
        row.get("release_id") != RELEASE_ID
        or row.get("gate") != "CANDIDATE_AUDIT_READY"
        or row.get("status") != "NON_GENETIC_COMPLETE_GENETICS_BLOCKED_UPSTREAM"
        for row in rows
    ):
        raise ValidationError("candidate audit seal is invalid")
    for row in rows:
        artifact = ROOT / row["artifact"]
        if not artifact.is_file() or sha256(artifact) != row["sha256"]:
            raise ValidationError(f"candidate audit artifact changed: {artifact}")


def atomic_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    if path.exists() or path.is_symlink():
        raise ValidationError(f"refusing to overwrite validation artifact: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def require_log(path: Path, markers: tuple[str, ...]) -> None:
    text = path.read_text(encoding="utf-8")
    if not all(marker in text for marker in markers):
        raise ValidationError(f"validation log lacks a required marker: {path}")


def main() -> None:
    args = arguments()
    pending = args.pending_root.resolve()
    genetics = (CANDIDATE / "genetics").resolve()
    if pending.parent != genetics or not pending.name.startswith(
        ".implementation_validation.pending."
    ):
        raise ValidationError(f"unsafe implementation validation path: {pending}")
    final = genetics / "implementation_validation"
    seal_path = CANDIDATE / "GENETICS_IMPLEMENTATION_READY"
    if final.exists() or final.is_symlink() or seal_path.exists() or seal_path.is_symlink():
        raise ValidationError("refusing to overwrite genetics implementation validation")

    expected_seals = {
        CANDIDATE / "NON_GENETIC_READY": args.non_genetic_ready_sha256,
        CANDIDATE / "CANDIDATE_AUDIT_READY": args.candidate_audit_ready_sha256,
        CANDIDATE / "integration/INTEGRATION_READY": args.integration_ready_sha256,
        CANDIDATE / "integration/FIG4F_READY": args.fig4f_ready_sha256,
    }
    for path, expected in expected_seals.items():
        if not path.is_file() or sha256(path) != expected:
            raise ValidationError(f"sealed non-genetic gate changed during validation: {path}")
    verify_ready_gate(CANDIDATE / "NON_GENETIC_READY", "NON_GENETIC_READY")
    verify_ready_gate(CANDIDATE / "integration/INTEGRATION_READY", "INTEGRATION_READY")
    verify_ready_gate(CANDIDATE / "integration/FIG4F_READY", "FIG4F_READY")
    verify_candidate_audit()

    runner = ROOT / "GWAS/finemapping/src/06_susie_coloc.R"
    if sha256(runner) != args.runner_before_sha256:
        raise ValidationError("active COLOC runner changed during genetics fixture validation")
    runner_text = runner.read_text(encoding="utf-8")
    if "ATAC_V3_" in runner_text:
        raise ValidationError("active COLOC runner contains candidate replay code")
    patch_targets = (
        "library(susieR)",
        'registry <- read.delim("config/gwas_registry.tsv", stringsAsFactors = FALSE)',
        'OUT_DIR <- file.path(FM_DIR, paste0("results/susie_coloc", OUT_SUFFIX), gwas_name)',
        "egenes <- unique(eqtl_all$ENSG)",
        'susie_method <- "susie"',
    )
    if any(runner_text.count(target) != 1 for target in patch_targets):
        raise ValidationError("active COLOC runner no longer supports an exact isolated patch")

    blocked = read(CANDIDATE / "genetics/BLOCKED_UPSTREAM_RELEASE.tsv")
    if len(blocked) != 1 or blocked[0].get("reason") != "promoted_corrected_coloc_manifest_unavailable":
        raise ValidationError("genetics upstream blocker is absent or changed")
    for forbidden in (
        CANDIDATE / "GENETIC_READY",
        CANDIDATE / "FULL_READY",
        genetics / "prepared",
        genetics / "replay_execution",
        genetics / "liftover",
        genetics / "context",
        genetics / "figures",
    ):
        if forbidden.exists() or forbidden.is_symlink():
            raise ValidationError(f"biological genetics stage exists before promotion: {forbidden}")

    logs = pending / "logs"
    require_log(logs / "static_checks.log", ("PY_COMPILE_PASS", "R_PARSE_PASS", "DIFF_CHECK_PASS"))
    require_log(logs / "fixture_tests.log", ("Ran ", "OK"))
    require_log(logs / "coloc_export_fixture.log", ("COLOC_SUSIE_EXPORT_FIXTURE\tPASS",))
    require_log(logs / "resource_scope.log", ("PASS",))

    source_rows = []
    for relative in SOURCES:
        path = PACKAGE / relative
        if not path.is_file() or path.is_symlink():
            raise ValidationError(f"missing or unsafe genetics implementation source: {path}")
        source_rows.append({
            "release_id": RELEASE_ID,
            "role": "genetics_implementation_source",
            "path": path.relative_to(ROOT).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    source_rows.append({
        "release_id": RELEASE_ID,
        "role": "immutable_active_coloc_runner",
        "path": runner.relative_to(ROOT).as_posix(),
        "bytes": runner.stat().st_size,
        "sha256": sha256(runner),
    })
    atomic_tsv(
        pending / "source_manifest.tsv",
        ("release_id", "role", "path", "bytes", "sha256"),
        source_rows,
    )
    completed = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    validation_rows = [
        {"release_id": RELEASE_ID, "check": check, "status": "PASS", "completed_utc": completed}
        for check in (
            "sealed_non_genetic_candidate_unchanged",
            "active_coloc_runner_unchanged",
            "isolated_runner_patch_fixture",
            "replay_batch_and_posterior_fixtures",
            "unique_liftover_and_allele_normalization_fixtures",
            "candidate_path_and_no_overwrite_contracts",
            "resource_scope",
            "upstream_promotion_fail_closed",
        )
    ]
    atomic_tsv(
        pending / "validation_results.tsv",
        ("release_id", "check", "status", "completed_utc"),
        validation_rows,
    )
    os.replace(pending, final)

    artifacts = [
        final / "source_manifest.tsv",
        final / "validation_results.tsv",
        final / "logs/static_checks.log",
        final / "logs/fixture_tests.log",
        final / "logs/coloc_export_fixture.log",
        final / "logs/resource_scope.log",
        CANDIDATE / "genetics/BLOCKED_UPSTREAM_RELEASE.tsv",
        *expected_seals.keys(),
    ]
    seal_rows = [{
        "release_id": RELEASE_ID,
        "gate": "GENETICS_IMPLEMENTATION_READY",
        "status": "FIXTURE_VALIDATED_UPSTREAM_BLOCKED",
        "artifact": path.relative_to(CANDIDATE).as_posix(),
        "sha256": sha256(path),
        "completed_utc": completed,
    } for path in artifacts]
    atomic_tsv(
        seal_path,
        ("release_id", "gate", "status", "artifact", "sha256", "completed_utc"),
        seal_rows,
    )
    print("GENETICS_IMPLEMENTATION_READY\tFIXTURE_VALIDATED_UPSTREAM_BLOCKED")


if __name__ == "__main__":
    main()

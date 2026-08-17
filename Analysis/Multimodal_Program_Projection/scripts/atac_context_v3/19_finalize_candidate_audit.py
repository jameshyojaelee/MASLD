#!/usr/bin/env python3
"""Independent terminal audit for the sealed, candidate-only ATAC v3 package."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import subprocess
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence


RELEASE_ID = "atac-context-v3-candidate-2026-08-11-r1"
ROOT = Path(__file__).resolve().parents[4]
PACKAGE = Path(__file__).resolve().parent
EXPECTED = (
    ROOT / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID
).resolve()
STANDARD_CHROMS = {f"chr{number}" for number in range(1, 23)} | {"chrX", "chrY"}
LINEAGES = {"hepatocyte", "stellate", "macrophage", "cholangiocyte", "t_nk"}
EXPECTED_PEAKS = {
    "hepatocyte": 206874,
    "stellate": 173191,
    "macrophage": 142456,
    "cholangiocyte": 113584,
    "t_nk": 40553,
}
EXPECTED_SHARED = {
    "hepatocyte": 131663,
    "stellate": 110048,
    "macrophage": 55300,
    "cholangiocyte": 45901,
    "t_nk": 4597,
}
PRODUCTION_JOBS = {
    "19750902": ("compute_smoke", "COMPLETED", "ACCEPTED"),
    "19750997": ("input_freeze", "COMPLETED", "ACCEPTED"),
    "19756191": ("chromvar_environment", "COMPLETED", "ACCEPTED"),
    "19774885": ("peak_support_smoke", "COMPLETED", "ACCEPTED"),
    "19776566": ("native_peak_calls_and_support", "FAILED", "RECOVERED_VALIDATED"),
    "19797475": ("consensus_finalization", "COMPLETED", "ACCEPTED"),
    "19797649": ("independent_consensus_validation", "COMPLETED", "ACCEPTED"),
    "19797819": ("gse244832_exact_counts", "COMPLETED", "ACCEPTED"),
    "19797820": ("gse281367_exact_counts", "COMPLETED", "ACCEPTED"),
    "19813420": ("da_and_program_inference", "COMPLETED", "ACCEPTED"),
    "19816710": ("official_chromvar", "COMPLETED", "ACCEPTED"),
    "19824258": ("independent_non_genetic_validation", "COMPLETED", "ACCEPTED"),
    "19827530": ("candidate_non_genetic_integration", "COMPLETED", "ACCEPTED"),
    "19830687": ("integration_compatibility_validation", "COMPLETED", "ACCEPTED"),
    "19832428": ("figure4f_render_and_validation", "COMPLETED", "ACCEPTED"),
}
DOCUMENTS = (
    "docs/README.md",
    "docs/STATUS.md",
    "docs/RESULTS.md",
    "docs/ROADMAP.md",
    "docs/technical/DATASETS_AND_PIPELINES.md",
)
POST_GATE_SOURCES = (
    "Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/15_integrate_non_genetic.py",
    "Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/16_validate_non_genetic_integration.py",
    "Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/17_render_integration_figure4f.R",
    "Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/18_validate_integration_figure4f.py",
    "Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/19_finalize_candidate_audit.py",
    "Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/tests/test_post_gate_contracts.py",
    "Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/slurm/10_integrate_non_genetic.sbatch",
    "Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/slurm/11_validate_non_genetic_integration.sbatch",
    "Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/slurm/12_render_validate_integration_figure4f.sbatch",
    "Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/slurm/13_finalize_candidate_audit.sbatch",
    "Analysis/ATAC/Integration/scripts/35_atac_integration.py",
)


class AuditError(RuntimeError):
    pass


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-root", type=Path, default=EXPECTED)
    parser.add_argument("--scope-log", type=Path, required=True)
    parser.add_argument("--tests-log", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    if not path.is_file():
        raise AuditError(f"missing required table: {path}")
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def atomic_tsv(
    path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    if path.exists() or path.is_symlink():
        raise AuditError(f"refusing to overwrite audit output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
        writer = csv.DictWriter(
            handle,
            fieldnames=list(columns),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def bh_adjust(pvalues: list[float]) -> list[float]:
    count = len(pvalues)
    order = sorted(range(count), key=lambda index: (pvalues[index], index))
    adjusted = [1.0] * count
    running = 1.0
    for order_index in range(count - 1, -1, -1):
        index = order[order_index]
        rank = order_index + 1
        running = min(running, pvalues[index] * count / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def require_close(observed: float, expected: float, label: str) -> None:
    if not math.isclose(observed, expected, rel_tol=1e-7, abs_tol=1e-10):
        raise AuditError(f"{label}: observed {observed}, expected {expected}")


def verify_gate(root: Path, relative: str, expected_gate: str) -> int:
    rows = read_rows(root / relative)
    if not rows or any(
        row.get("release_id") != RELEASE_ID
        or row.get("gate") != expected_gate
        or row.get("status") != "READY"
        for row in rows
    ):
        raise AuditError(f"invalid {expected_gate} seal")
    for row in rows:
        artifact_text = row["artifact"]
        artifact = root / artifact_text
        if not artifact.is_file():
            artifact = ROOT / artifact_text
        if not artifact.is_file() or sha256(artifact) != row["sha256"]:
            raise AuditError(f"{expected_gate} artifact hash mismatch: {artifact_text}")
    return len(rows)


def audit_manifest(root: Path) -> tuple[int, int]:
    rows = read_rows(root / "input_manifest.tsv")
    changed_sources = []
    for row in rows:
        path = Path(row["relative_or_absolute_path"])
        artifact = path if path.is_absolute() else ROOT / path
        if not artifact.is_file():
            raise AuditError(f"frozen input missing: {artifact}")
        observed = sha256(artifact)
        source_changed = observed != row["sha256"]
        if source_changed:
            if row["role"] == "producer_source" and path.as_posix().endswith(
                "/15_integrate_non_genetic.py"
            ):
                changed_sources.append(path.as_posix())
            else:
                raise AuditError(f"frozen input changed: {artifact}")
        if int(row["bytes"]) != artifact.stat().st_size and not source_changed:
            raise AuditError(f"frozen input size changed: {artifact}")
    if changed_sources != [
        "Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/15_integrate_non_genetic.py"
    ]:
        raise AuditError(f"unexpected changed sealed producer set: {changed_sources}")
    return len(rows), len(changed_sources)


def audit_peaks_and_counts(root: Path) -> tuple[int, int, int]:
    rows = read_rows(root / "consensus_peak_manifest.tsv")
    by_lineage = Counter(row["lineage"] for row in rows)
    shared = Counter(row["lineage"] for row in rows if row["supported_both"] == "TRUE")
    if dict(by_lineage) != EXPECTED_PEAKS or dict(shared) != EXPECTED_SHARED:
        raise AuditError(f"peak census mismatch: {dict(by_lineage)}; shared={dict(shared)}")
    coordinates = set()
    for row in rows:
        if row["release_id"] != RELEASE_ID or row["chrom"] not in STANDARD_CHROMS:
            raise AuditError("peak release or chromosome contract failed")
        if int(row["width"]) != 500 or int(row["end"]) - int(row["start0"]) != 500:
            raise AuditError("peak width contract failed")
        if row["blacklist_overlap"] != "FALSE":
            raise AuditError("blacklist overlap in consensus atlas")
        key = (row["lineage"], row["chrom"], row["start0"], row["end"])
        if key in coordinates:
            raise AuditError(f"duplicate consensus coordinate: {key}")
        coordinates.add(key)

    matrices = sorted((root / "counts").glob("*/*.mtx.gz"))
    if len(matrices) != 10:
        raise AuditError(f"primary sparse matrix family size is {len(matrices)}, expected 10")
    for matrix in matrices:
        cohort = matrix.parent.name
        lineage = matrix.name.removesuffix(".mtx.gz")
        donor_path = matrix.with_name(f"{lineage}.donors.tsv")
        peak_path = matrix.with_name(f"{lineage}.peaks.tsv")
        with gzip.open(matrix, "rt", encoding="utf-8") as handle:
            dimensions = next(line for line in handle if not line.startswith("%")).split()
        n_donors, n_peaks, _ = (int(value) for value in dimensions)
        if n_donors != len(read_rows(donor_path)) or n_peaks != len(read_rows(peak_path)):
            raise AuditError(f"matrix dimension mismatch: {cohort}/{lineage}")
        if n_peaks != EXPECTED_PEAKS[lineage]:
            raise AuditError(f"matrix peak family mismatch: {cohort}/{lineage}")
    recount = read_rows(root / "recount/exact_fragment_recount.tsv")
    if len(recount) != 10 or any(
        row["record_count_match"] != "TRUE"
        or row["counting_unit"] != "deduplicated_fragment_record"
        for row in recount
    ):
        raise AuditError("exact fragment recount contract failed")
    return len(rows), len(matrices), len(recount)


def audit_condition_blindness(root: Path) -> int:
    rows = read_rows(root / "condition_blindness/condition_blindness_audit.tsv")
    if len(rows) != 10 or any(
        row["pass"] != "TRUE"
        or row["condition_permutation_seed"] != "42"
        or row["peak_coordinate_hash_before"] != row["peak_coordinate_hash_after"]
        or row["filter_hash_before"] != row["filter_hash_after"]
        or row["promoter_matrix_hash_before"] != row["promoter_matrix_hash_after"]
        or row["program_testability_hash_before"] != row["program_testability_hash_after"]
        for row in rows
    ):
        raise AuditError("condition-blindness audit failed")
    return len(rows)


def audit_da(root: Path) -> tuple[int, Counter[str], int]:
    summary = read_rows(root / "da/da_lineage_summary.tsv")
    if {row["lineage"] for row in summary} != {"hepatocyte", "stellate", "macrophage"}:
        raise AuditError("DA lineage family mismatch")
    expected_joint = {"hepatocyte": 30650, "stellate": 9264, "macrophage": 0}
    for row in summary:
        lineage = row["lineage"]
        if int(row["n_jointly_testable"]) != expected_joint[lineage]:
            raise AuditError(f"DA denominator mismatch: {lineage}")
        if lineage == "macrophage" and row["lineage_state"] != "untestable":
            raise AuditError("macrophage DA should remain untestable")
        if lineage == "stellate" and (
            int(row["n_supported"]) != 4
            or row["supported_sensitivity_stable"] != "TRUE"
            or row["dynamic_main_eligible"] != "TRUE"
        ):
            raise AuditError("stellate DA support gate mismatch")

    q_families = 0
    for cohort in ("GSE244832", "GSE281367"):
        for lineage in ("hepatocyte", "stellate", "macrophage"):
            rows = read_rows(root / f"da/cohort_primary_{cohort}_{lineage}.tsv")
            if not rows:
                continue
            recalculated = bh_adjust([float(row["pvalue"]) for row in rows])
            for index, (row, expected) in enumerate(zip(rows, recalculated)):
                require_close(float(row["qvalue"]), expected, f"DA BH {cohort}/{lineage}/{index}")
            q_families += 1

    joined = read_rows(root / "da/da_peak_results.tsv.gz")
    states = Counter(row["evidence_state"] for row in joined)
    if len(joined) != 39914 or states != {
        "indeterminate": 38066,
        "source_dependent": 1844,
        "supported": 4,
    }:
        raise AuditError(f"joined DA state family mismatch: {states}")
    if any(row["sensitivity_no_reversal"] != "TRUE" for row in joined if row["evidence_state"] == "supported"):
        raise AuditError("supported DA peak reverses in sensitivity")

    pvalues = [float(row["directional_binomial_pvalue"]) for row in summary]
    qvalues = bh_adjust(pvalues)
    for row, expected in zip(summary, qvalues):
        require_close(float(row["directional_binomial_qvalue"]), expected, "DA directional BH")
    return len(joined), states, q_families + 1


def audit_programs(root: Path) -> tuple[int, Counter[str], int]:
    rows = read_rows(root / "programs/program_atac_results.tsv")
    keys = [(row["cohort"], row["program_uid"]) for row in rows]
    if len(rows) != 234 or len(set(keys)) != 234:
        raise AuditError("program family is not exactly 117 programs by two cohorts")
    cohort_counts = Counter(row["cohort"] for row in rows)
    if cohort_counts != {"GSE244832": 117, "GSE281367": 117}:
        raise AuditError(f"program cohort family mismatch: {cohort_counts}")
    input_rows = read_rows(root / "input_manifest.tsv")
    registry_path = next(
        Path(row["relative_or_absolute_path"])
        for row in input_rows if row["role"] == "program_registry"
    )
    registry = read_rows(ROOT / registry_path)
    registry_hashes = {row["program_uid"]: row["membership_sha256"] for row in registry}
    if len(registry_hashes) != 117:
        raise AuditError("frozen registry does not contain 117 unique programs")
    for row in rows:
        if registry_hashes.get(row["program_uid"]) != row["membership_sha256"]:
            raise AuditError(f"program membership hash drift: {row['program_uid']}")
    for cohort in ("GSE244832", "GSE281367"):
        tested = [
            row for row in rows
            if row["cohort"] == cohort and row["contrast_testable"] == "TRUE"
        ]
        recalculated = bh_adjust([float(row["pvalue"]) for row in tested])
        for index, (row, expected) in enumerate(zip(tested, recalculated)):
            require_close(float(row["qvalue"]), expected, f"program BH {cohort}/{index}")
    states = Counter(row["cross_cohort_state"] for row in rows)
    if states != {"indeterminate": 156, "untestable": 78}:
        raise AuditError(f"program cross-cohort state mismatch: {states}")
    if any(row["replicated"] == "TRUE" for row in rows):
        raise AuditError("unexpected replicated program")
    return len(rows), states, 2


def audit_chromvar(root: Path) -> tuple[int, Counter[str], int]:
    rows = read_rows(root / "chromvar/chromvar_results.tsv.gz")
    if len(rows) != 8790:
        raise AuditError(f"chromVAR result family is {len(rows)}, expected 8790")
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["contrast_testable"] == "TRUE":
            grouped[(row["cohort"], row["lineage"])].append(row)
    if len(grouped) != 8 or sum(map(len, grouped.values())) != 7032:
        raise AuditError("chromVAR tested-family denominator mismatch")
    for key, family in grouped.items():
        recalculated = bh_adjust([float(row["pvalue"]) for row in family])
        for index, (row, expected) in enumerate(zip(family, recalculated)):
            require_close(float(row["qvalue"]), expected, f"chromVAR BH {key}/{index}")
    session = (root / "chromvar/sessionInfo.txt").read_text(encoding="utf-8")
    for package in ("chromVAR_1.32.0", "motifmatchr_1.32.0", "JASPAR2024_0.99.7"):
        if package not in session:
            raise AuditError(f"official chromVAR session missing {package}")
    states = Counter(row["cross_cohort_state"] for row in rows)
    if states.get("supported", 0) != 0 or states.get("replicated", 0) != 0:
        raise AuditError("unexpected replicated chromVAR motif")
    return len(rows), states, len(grouped)


def audit_integration(root: Path) -> tuple[int, int, int]:
    states = read_rows(root / "integration/fig4f_atac_v3_state_contract.tsv")
    if len(states) != 234:
        raise AuditError("Figure 4F state contract row count mismatch")
    source = read_rows(root / "integration/figures/fig4f_atac_v3_contract_source.tsv")
    if len(source) != 30:
        raise AuditError("Figure 4F source metric grid mismatch")
    before = root / "compatibility/default_before/l8_atac_columns.csv"
    after = root / "compatibility/default_after/l8_atac_columns.csv"
    if before.read_bytes() != after.read_bytes():
        raise AuditError("default integration compatibility output changed")
    v3 = root / "compatibility/with_v3/l8_atac_columns.csv"
    with v3.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 27187:
        raise AuditError("candidate Gene Catalog adapter row count mismatch")
    mapped = sum(row.get("atac_v3_release_id") == RELEASE_ID for row in rows)
    if mapped != 14510:
        raise AuditError(f"candidate Gene Catalog mapping denominator changed: {mapped}")
    return len(states), len(source), mapped


def audit_genetics_gate(root: Path) -> str:
    blocked = read_rows(root / "genetics/BLOCKED_UPSTREAM_RELEASE.tsv")
    if len(blocked) != 1 or blocked[0].get("status") != "blocked_upstream_release":
        raise AuditError("genetics fail-closed gate is invalid")
    if (root / "GENETIC_READY").exists() or (root / "FULL_READY").exists():
        raise AuditError("genetics or full seal exists without promoted upstream release")
    if (root / "genetics/replay_plan.tsv").exists():
        raise AuditError("biological genetics replay rows exist while upstream is blocked")
    return blocked[0]["status"]


def audit_jobs(root: Path) -> list[dict[str, str]]:
    command = [
        "sacct", "-n", "-P", "-j", ",".join(PRODUCTION_JOBS),
        "--format=JobIDRaw,JobName,Partition,State,Elapsed,ReqMem,ExitCode",
    ]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    observed = {}
    for line in result.stdout.splitlines():
        fields = line.split("|")
        if len(fields) < 7 or "." in fields[0] or fields[0] not in PRODUCTION_JOBS:
            continue
        observed[fields[0]] = fields
    if set(observed) != set(PRODUCTION_JOBS):
        raise AuditError(f"production job provenance incomplete: {set(PRODUCTION_JOBS) - set(observed)}")
    rows = []
    recovery = read_rows(root / "failed_attempts/consensus_merge_recovery.tsv")
    recovered = next((row for row in recovery if row["job_id"] == "19776566"), None)
    if recovered is None or recovered["biological_outputs_used"] != "cohort-native peak tables and support audit only":
        raise AuditError("partial peak-call job recovery provenance is missing")
    for job_id, (stage, expected_state, disposition) in PRODUCTION_JOBS.items():
        values = observed[job_id]
        if values[3] != expected_state:
            raise AuditError(f"production job state drift: {job_id} {values[3]}")
        if expected_state == "COMPLETED" and values[6] != "0:0":
            raise AuditError(f"completed production job has nonzero exit: {job_id} {values[6]}")
        if expected_state == "FAILED" and values[6] == "0:0":
            raise AuditError(f"failed recovery source has zero exit: {job_id}")
        rows.append({
            "release_id": RELEASE_ID,
            "stage": stage,
            "job_id": job_id,
            "job_name": values[1],
            "partition": values[2],
            "state": values[3],
            "elapsed": values[4],
            "requested_memory": values[5],
            "exit_code": values[6],
            "disposition": disposition,
        })
    return rows


def audit_docs_and_sources(root: Path) -> list[dict[str, object]]:
    results = (ROOT / "docs/RESULTS.md").read_text(encoding="utf-8")
    status = (ROOT / "docs/STATUS.md").read_text(encoding="utf-8")
    roadmap = (ROOT / "docs/ROADMAP.md").read_text(encoding="utf-8")
    required_results = (
        "676,658",
        "39,914",
        "All four supported",
        "exactly 234",
        "8,790",
        "candidate",
        "not promoted",
    )
    if any(value not in results for value in required_results):
        raise AuditError("RESULTS.md lacks sealed ATAC candidate facts")
    if "genetics blocked upstream" not in status.lower():
        raise AuditError("STATUS.md lacks genetics blocker")
    if "ATAC v3 variant-posterior replay" not in roadmap:
        raise AuditError("ROADMAP.md lacks remaining variant replay work")

    rows = []
    for relative in DOCUMENTS + POST_GATE_SOURCES:
        path = ROOT / relative
        if not path.is_file():
            raise AuditError(f"post-gate source missing: {relative}")
        rows.append({
            "release_id": RELEASE_ID,
            "role": "human_authority" if relative.startswith("docs/") else "post_gate_source",
            "path": relative,
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    diff = root / "integration/35_atac_integration.v3.diff"
    rows.append({
        "release_id": RELEASE_ID,
        "role": "integration_diff",
        "path": str(diff.relative_to(ROOT)),
        "bytes": diff.stat().st_size,
        "sha256": sha256(diff),
    })
    return rows


def main() -> None:
    args = arguments()
    supplied = args.candidate_root
    root = supplied.resolve()
    if root != EXPECTED or supplied.is_symlink():
        raise AuditError(f"unsafe candidate root: {root}")
    audit_dir = root / "audit"
    seal = root / "CANDIDATE_AUDIT_READY"
    if audit_dir.exists() or seal.exists():
        raise AuditError("refusing to overwrite terminal audit")
    if not args.scope_log.is_file() or not args.tests_log.is_file():
        raise AuditError("compute validation logs are required")
    scope_text = args.scope_log.read_text(encoding="utf-8")
    tests_text = args.tests_log.read_text(encoding="utf-8")
    if "PASS" not in scope_text or "OK" not in tests_text:
        raise AuditError("scope or fixture test log does not contain a passing terminal marker")

    checks: list[dict[str, object]] = []

    def passed(check: str, details: Mapping[str, object]) -> None:
        checks.append({
            "release_id": RELEASE_ID,
            "check": check,
            "status": "PASS",
            "details_json": json.dumps(details, sort_keys=True),
        })

    non_genetic_seal_rows = verify_gate(root, "NON_GENETIC_READY", "NON_GENETIC_READY")
    integration_seal_rows = verify_gate(root, "integration/INTEGRATION_READY", "INTEGRATION_READY")
    figure_seal_rows = verify_gate(root, "integration/FIG4F_READY", "FIG4F_READY")
    passed("sealed_gates", {
        "non_genetic_rows": non_genetic_seal_rows,
        "integration_rows": integration_seal_rows,
        "figure4f_rows": figure_seal_rows,
    })

    manifest_rows, changed_sources = audit_manifest(root)
    passed("frozen_inputs", {"rows": manifest_rows, "documented_post_gate_changes": changed_sources})
    peaks, matrices, recounts = audit_peaks_and_counts(root)
    passed("peaks_counts_recount", {"peaks": peaks, "primary_matrices": matrices, "exact_recounts": recounts})
    passed("condition_blindness", {"families": audit_condition_blindness(root), "seed": 42})
    da_rows, da_states, da_bh = audit_da(root)
    passed("differential_accessibility", {"rows": da_rows, "states": dict(da_states), "bh_families": da_bh})
    program_rows, program_states, program_bh = audit_programs(root)
    passed("program_projection", {"rows": program_rows, "states": dict(program_states), "bh_families": program_bh})
    motif_rows, motif_states, motif_bh = audit_chromvar(root)
    passed("official_chromvar", {"rows": motif_rows, "states": dict(motif_states), "bh_families": motif_bh})
    fig_states, fig_cells, mapped_genes = audit_integration(root)
    passed("candidate_only_integration", {"program_states": fig_states, "figure_cells": fig_cells, "mapped_catalog_genes": mapped_genes})
    genetics_status = audit_genetics_gate(root)
    passed("genetics_fail_closed", {"status": genetics_status, "genetic_ready": False, "full_ready": False})
    post_gate_rows = audit_docs_and_sources(root)
    passed("documentation_and_sources", {"post_gate_rows": len(post_gate_rows)})
    job_rows = audit_jobs(root)
    passed("slurm_provenance", {
        "terminal_jobs": len(job_rows),
        "accepted_jobs": sum(row["disposition"] == "ACCEPTED" for row in job_rows),
        "recovered_validated_jobs": sum(row["disposition"] == "RECOVERED_VALIDATED" for row in job_rows),
    })
    passed("compute_checks", {
        "scope_log": str(args.scope_log.relative_to(ROOT)),
        "scope_sha256": sha256(args.scope_log),
        "tests_log": str(args.tests_log.relative_to(ROOT)),
        "tests_sha256": sha256(args.tests_log),
    })

    atomic_tsv(
        audit_dir / "completion_audit.tsv",
        ("release_id", "check", "status", "details_json"),
        checks,
    )
    atomic_tsv(
        audit_dir / "post_gate_source_manifest.tsv",
        ("release_id", "role", "path", "bytes", "sha256"),
        post_gate_rows,
    )
    atomic_tsv(
        audit_dir / "job_audit.tsv",
        ("release_id", "stage", "job_id", "job_name", "partition", "state", "elapsed", "requested_memory", "exit_code", "disposition"),
        job_rows,
    )
    completed = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    seal_rows = []
    for path in (
        audit_dir / "completion_audit.tsv",
        audit_dir / "post_gate_source_manifest.tsv",
        audit_dir / "job_audit.tsv",
        args.scope_log,
        args.tests_log,
    ):
        seal_rows.append({
            "release_id": RELEASE_ID,
            "gate": "CANDIDATE_AUDIT_READY",
            "status": "NON_GENETIC_COMPLETE_GENETICS_BLOCKED_UPSTREAM",
            "artifact": str(path.relative_to(ROOT)),
            "sha256": sha256(path),
            "completed_utc": completed,
        })
    atomic_tsv(
        seal,
        ("release_id", "gate", "status", "artifact", "sha256", "completed_utc"),
        seal_rows,
    )
    print(
        f"CANDIDATE_AUDIT_READY: {len(checks)} independent checks passed; "
        "non-genetic candidate complete; genetics remains blocked upstream"
    )


if __name__ == "__main__":
    main()

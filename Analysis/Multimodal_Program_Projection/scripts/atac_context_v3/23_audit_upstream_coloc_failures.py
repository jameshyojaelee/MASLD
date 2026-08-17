#!/usr/bin/env python3
"""Read-only structural audit of failed corrected-COLOC array task outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
RELEASE_ID = "atac-context-v3-candidate-2026-08-11-r1"
CANDIDATE = (
    ROOT / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID
).resolve()
ARRAY_JOB_ID = "19648525"
EXPECTED_TASKS = 1100
EXPECTED_STUDIES = 50
REQUIRED_COLUMNS = {
    "gwas_name", "gene", "ensembl", "chr", "PP.H4.abf", "PP.H4.susie",
    "method", "n_snps",
}
ALLOWED_METHODS = {"susie", "abf_fallback", "abf_only"}


class AuditError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, columns: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_probability(value: str) -> bool:
    if value in {"", "NA", "NaN"}:
        return True
    try:
        number = float(value)
    except ValueError:
        return False
    return math.isfinite(number) and 0 <= number <= 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    monitor_root = (CANDIDATE / "upstream_monitor").resolve()
    if output.parent != monitor_root or not output.name.startswith("coloc-failure-audit-"):
        raise AuditError(f"unsafe upstream audit output: {output}")
    if output.exists() or output.is_symlink():
        raise AuditError(f"refusing to overwrite upstream audit: {output}")
    stale = sorted(monitor_root.glob(f".{output.name}.pending.*")) if monitor_root.exists() else []
    if stale:
        raise AuditError(f"unresolved prior upstream audit: {stale}")
    pending = monitor_root / f".{output.name}.pending.{os.getpid()}"
    pending.mkdir(parents=True)

    registry = ROOT / "GWAS/finemapping/config/gwas_registry.tsv"
    with registry.open("r", encoding="utf-8", newline="") as handle:
        studies = sorted({row["study_name"] for row in csv.DictReader(handle, delimiter="\t")})
    if len(studies) != EXPECTED_STUDIES:
        raise AuditError(f"expected {EXPECTED_STUDIES} registry studies, observed {len(studies)}")

    accounting = subprocess.run(
        [
            "sacct", "-X", "-j", ARRAY_JOB_ID, "--starttime", "2026-08-08",
            "-n", "-P", "--format=JobID,State,Start,End,ExitCode",
        ],
        check=True, text=True, capture_output=True,
    ).stdout.splitlines()
    failed = []
    for line in accounting:
        fields = line.split("|")
        if len(fields) < 5 or not fields[0].startswith(f"{ARRAY_JOB_ID}_"):
            continue
        task_text = fields[0].split("_", 1)[1]
        if not task_text.isdigit():
            continue
        task = int(task_text)
        if fields[1].split()[0] == "FAILED":
            failed.append((task, fields[2], fields[3], fields[4]))
    if not failed:
        raise AuditError("no failed corrected-COLOC tasks were reported by Slurm")

    rows = []
    for task, start, end, exit_code in sorted(failed):
        if not 0 <= task < EXPECTED_TASKS:
            raise AuditError(f"array task outside 0..1099: {task}")
        study = studies[task // 22]
        chromosome = task % 22 + 1
        result = (
            ROOT / "GWAS/finemapping/results/susie_coloc_rerun" / study
            / f"susie_coloc_chr{chromosome}.csv"
        )
        stderr = (
            ROOT / "GWAS/finemapping/logs/coloc_rerun"
            / f"coloc_{ARRAY_JOB_ID}_{task}.err"
        )
        stdout = stderr.with_suffix(".out")
        reasons = []
        n_rows = 0
        methods: set[str] = set()
        genes: set[str] = set()
        anonymous_abf_only_rows = 0
        header: list[str] = []
        if not result.is_file() or result.stat().st_size == 0:
            reasons.append("missing_or_empty_result")
        else:
            with result.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle)
                header = list(reader.fieldnames or [])
                if not REQUIRED_COLUMNS.issubset(header):
                    reasons.append("schema_incomplete")
                for row in reader:
                    n_rows += 1
                    if row.get("gwas_name") != study:
                        reasons.append("gwas_identity_mismatch")
                    try:
                        observed_chr = int(float(row.get("chr", "")))
                    except ValueError:
                        observed_chr = -1
                    if observed_chr != chromosome:
                        reasons.append("chromosome_identity_mismatch")
                    gene_id = row.get("ensembl", "")
                    method = row.get("method", "")
                    methods.add(method)
                    if method not in ALLOWED_METHODS:
                        reasons.append("invalid_method")
                    if not gene_id:
                        susie_value = row.get("PP.H4.susie", "")
                        try:
                            abf_value = float(row.get("PP.H4.abf", "nan"))
                        except ValueError:
                            abf_value = float("nan")
                        if (
                            method == "abf_only"
                            and susie_value in {"", "NA", "NaN"}
                            and math.isfinite(abf_value)
                            and abf_value <= 0.5
                        ):
                            anonymous_abf_only_rows += 1
                        else:
                            reasons.append("anonymous_row_not_excludable")
                    elif gene_id in genes:
                        reasons.append("duplicate_ensembl")
                    else:
                        genes.add(gene_id)
                    if not parse_probability(row.get("PP.H4.abf", "")):
                        reasons.append("invalid_abf_probability")
                    if not parse_probability(row.get("PP.H4.susie", "")):
                        reasons.append("invalid_susie_probability")
            if n_rows == 0:
                reasons.append("zero_result_rows")
        error_text = stderr.read_text(encoding="utf-8", errors="replace") if stderr.is_file() else ""
        output_text = stdout.read_text(encoding="utf-8", errors="replace") if stdout.is_file() else ""
        expected_postwrite_error = "Error in ve(chk) : could not find function \"ve\"" in error_text
        if not expected_postwrite_error:
            reasons.append("failure_signature_differs")
        if "Results:" not in output_text or "Genes tested:" not in output_text:
            reasons.append("completion_summary_absent")
        reasons = sorted(set(reasons))
        rows.append({
            "release_id": RELEASE_ID,
            "array_job_id": ARRAY_JOB_ID,
            "array_task_id": task,
            "study": study,
            "chromosome": chromosome,
            "slurm_start": start,
            "slurm_end": end,
            "slurm_exit_code": exit_code,
            "result_path": str(result),
            "result_bytes": result.stat().st_size if result.is_file() else 0,
            "result_sha256": sha256(result) if result.is_file() else "",
            "n_rows": n_rows,
            "n_unique_ensembl": len(genes),
            "n_anonymous_excluded_abf_only_rows": anonymous_abf_only_rows,
            "methods": ",".join(sorted(methods)),
            "expected_postwrite_error": str(expected_postwrite_error).upper(),
            "structural_status": (
                (
                    "complete_output_with_anonymous_excluded_abf_only_row"
                    if anonymous_abf_only_rows > 0
                    else "complete_output_after_postwrite_failure"
                )
                if not reasons else "requires_upstream_review"
            ),
            "review_reasons": ",".join(reasons),
        })

    status_counts: dict[str, int] = {}
    for row in rows:
        status = str(row["structural_status"])
        status_counts[status] = status_counts.get(status, 0) + 1
    n_structurally_complete = sum(
        count for status, count in status_counts.items() if status.startswith("complete_output_")
    )
    completed = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    write_tsv(
        pending / "failed_task_output_audit.tsv",
        tuple(rows[0]), rows,
    )
    write_tsv(
        pending / "audit_summary.tsv",
        ("release_id", "status", "n_failed_tasks", "n_structurally_complete", "n_review", "completed_utc"),
        [{
            "release_id": RELEASE_ID,
            "status": "READ_ONLY_RECOVERY_EVIDENCE_NOT_PROMOTION",
            "n_failed_tasks": len(rows),
            "n_structurally_complete": n_structurally_complete,
            "n_review": status_counts.get("requires_upstream_review", 0),
            "completed_utc": completed,
        }],
    )
    write_tsv(
        pending / "input_manifest.tsv",
        ("release_id", "role", "path", "sha256"),
        [
            {"release_id": RELEASE_ID, "role": "gwas_registry", "path": str(registry), "sha256": sha256(registry)},
            {
                "release_id": RELEASE_ID, "role": "active_coloc_runner",
                "path": str(ROOT / "GWAS/finemapping/src/06_susie_coloc.R"),
                "sha256": sha256(ROOT / "GWAS/finemapping/src/06_susie_coloc.R"),
            },
        ],
    )
    os.replace(pending, output)
    print(
        "COLOC_FAILURE_OUTPUT_AUDIT\t"
        f"failed={len(rows)}\tstructurally_complete="
        f"{n_structurally_complete}\t"
        f"review={status_counts.get('requires_upstream_review', 0)}"
    )


if __name__ == "__main__":
    main()

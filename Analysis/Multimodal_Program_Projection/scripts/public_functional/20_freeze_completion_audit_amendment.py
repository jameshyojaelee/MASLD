#!/usr/bin/env python3
"""Freeze the post-validation method-compliance audit before remediation.

This is an append-only correction contract. It does not alter the outcome-blind
scientific hypotheses, frozen programs, evidence classes, contrasts, or
promotion thresholds. It preserves the previous terminal tables and code so
the method-compliant rerun cannot erase the flawed-but-historical bundle.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import shutil
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[4]
CANDIDATE_ROOT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09"
SCRIPT_ROOT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/scripts/public_functional"
ARCHIVE_ROOT = CANDIDATE_ROOT / "superseded_bundles/2026-08-09_pre_completion_audit_method_fix"
EXPECTED_SPEC = "50c204c596a967b7565e88c562a418476af6a4c056e9d3ac1c9273acd7369169"
EXPECTED_RELEASE = "3e92aaf179868f143a02c5034c30dcb1b687d4b5bf99feb1cbaa44edd8ca1fc9"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(text)
    temporary.replace(path)


def main() -> None:
    validated_path = CANDIDATE_ROOT / "VALIDATED.json"
    if not validated_path.is_file():
        raise RuntimeError("The pre-remediation terminal bundle is missing VALIDATED.json")
    validated = json.loads(validated_path.read_text())
    if validated.get("specification_sha256") != EXPECTED_SPEC:
        raise RuntimeError("Specification drift before completion-audit freeze")
    if validated.get("release_manifest_sha256") != EXPECTED_RELEASE:
        raise RuntimeError("Unexpected pre-remediation release identity")
    if sha256(CANDIDATE_ROOT / "release_manifest.tsv") != EXPECTED_RELEASE:
        raise RuntimeError("Pre-remediation release_manifest.tsv hash mismatch")
    if ARCHIVE_ROOT.exists():
        raise RuntimeError(f"Archive already exists: {ARCHIVE_ROOT}")

    candidate_files = [
        "VALIDATED.json",
        "release_manifest.tsv",
        "validation_report.tsv",
        "figure5_verdict.tsv",
        "integration_status.json",
        "program_effects.tsv",
        "program_testability.tsv",
        "sensitivity.tsv",
        "evidence_class_effects.tsv",
        "evidence_class_sensitivity.tsv",
        "cross_assay_support.tsv",
        "null_permutation_audit.tsv",
        "source_gate_status.tsv",
        "source_schema_audit.tsv",
    ]
    analysis_files = []
    for dataset in ["GSE200418", "GSE207889", "GSE253380", "GSE106737"]:
        for filename in [
            "program_effects.tsv",
            "program_testability.tsv",
            "sensitivity.tsv",
            "evidence_class_effects.tsv",
            "evidence_class_sensitivity.tsv",
            "gene_level_effects.tsv",
        ]:
            relative = Path("analyses") / dataset / filename
            if (CANDIDATE_ROOT / relative).is_file():
                analysis_files.append(str(relative))
    code_files = [
        "03_run_source_gates.py",
        "04_run_bulk_assays.R",
        "05_build_schlo_pseudobulk.py",
        "06_run_schlo_programs.R",
        "07_run_evidence_classes.R",
        "11_run_null_permutations.py",
        "12_build_final_integration.py",
        "13_validate_final.py",
        "functional_analysis_common.R",
        "test_public_functional.R",
    ]

    manifest = []
    for relative_text in candidate_files + analysis_files:
        relative = Path(relative_text)
        source = CANDIDATE_ROOT / relative
        if not source.is_file():
            raise RuntimeError(f"Missing required pre-remediation artifact: {source}")
        destination = ARCHIVE_ROOT / "candidate" / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        manifest.append(
            {
                "kind": "candidate_artifact",
                "source": str(source.relative_to(PROJECT_ROOT)),
                "archive": str(destination.relative_to(PROJECT_ROOT)),
                "size_bytes": source.stat().st_size,
                "sha256": sha256(source),
            }
        )
    for filename in code_files:
        source = SCRIPT_ROOT / filename
        destination = ARCHIVE_ROOT / "code" / filename
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        manifest.append(
            {
                "kind": "producer_code",
                "source": str(source.relative_to(PROJECT_ROOT)),
                "archive": str(destination.relative_to(PROJECT_ROOT)),
                "size_bytes": source.stat().st_size,
                "sha256": sha256(source),
            }
        )

    manifest_path = ARCHIVE_ROOT / "archive_manifest.tsv"
    lines = ["kind\tsource\tarchive\tsize_bytes\tsha256"]
    for row in manifest:
        lines.append("\t".join(str(row[key]) for key in ["kind", "source", "archive", "size_bytes", "sha256"]))
    atomic_text(manifest_path, "\n".join(lines) + "\n")

    amendment = {
        "amendment_id": "completion-audit-method-compliance-v1",
        "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "candidate_id": "public-functional-map-2026-08-09",
        "specification_sha256": EXPECTED_SPEC,
        "superseded_release_manifest_sha256": EXPECTED_RELEASE,
        "archive_manifest_sha256": sha256(manifest_path),
        "audit_findings": [
            {
                "requirement": "GSE207889 replicate-by-lineage pseudobulks use edgeR/limma-voom",
                "observed_violation": "06_run_schlo_programs.R used logCPM plus a two-value t-test rather than voom and a paired limma model",
                "required_correction": "use voom on pseudobulk counts and fit condition contrasts with deposited replicate as a blocking fixed effect; retain replicate directions as sensitivity",
            },
            {
                "requirement": "GSE106737 highest-median-expression probe primary with median-across-probes sensitivity",
                "observed_violation": "the primary collapse was implemented but no median-across-probes sensitivity was emitted",
                "required_correction": "emit program and evidence-class median-across-probes sensitivities without changing the primary mapping rule",
            },
            {
                "requirement": "source and exclusion accounting is validated, not merely narrated",
                "observed_violation": "the PCLS gate did not assert the source-stated 112 attempted / 13 failed / 99 deposited accounting, and the final validator did not test it",
                "required_correction": "parse and validate the local GEO SOFT statement and add fail-closed terminal checks",
            },
            {
                "requirement": "final validation covers the preregistered assay-native methods and mandatory sensitivities",
                "observed_violation": "149 passing checks covered table presence and multiplicity but not the two method requirements above",
                "required_correction": "add explicit model/sensitivity/source-accounting checks and rerun terminal validation",
            },
        ],
        "unchanged_scientific_contract": [
            "117 frozen programs and memberships",
            "two primary program UIDs and weights",
            "evidence-class definitions",
            "all contrasts and expected directions",
            "multiplicity families",
            "main-figure promotion gates",
            "no alternative dataset, contrast, or outcome-selected subset",
        ],
        "allowed_code_changes": [
            "03_run_source_gates.py",
            "04_run_bulk_assays.R",
            "06_run_schlo_programs.R",
            "07_run_evidence_classes.R",
            "12_build_final_integration.py",
            "13_validate_final.py",
            "test_public_functional.R",
            "new remediation freeze/finalization wrappers and documentation",
        ],
        "prohibited_changes": [
            "program discovery, membership, weight, label, or expected-direction changes",
            "new datasets or contrasts",
            "relaxed q-value, donor-agreement, robustness, or cross-assay gates",
            "deletion or mutation of the archived pre-remediation bundle",
            "canonical output or main-figure writes",
        ],
    }
    amendment_path = CANDIDATE_ROOT / "COMPLETION_AUDIT_AMENDMENT_01.json"
    atomic_text(amendment_path, json.dumps(amendment, indent=2, sort_keys=True) + "\n")
    atomic_text(CANDIDATE_ROOT / "COMPLETION_AUDIT_AMENDMENT_01_SHA256", sha256(amendment_path) + "\n")
    print(json.dumps({"status": "frozen", "amendment_sha256": sha256(amendment_path), "archived_files": len(manifest)}, indent=2))


if __name__ == "__main__":
    main()

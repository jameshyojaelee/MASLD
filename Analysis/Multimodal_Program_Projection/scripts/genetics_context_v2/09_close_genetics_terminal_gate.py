#!/usr/bin/env python3
"""Create a hashed terminal GEN closure without rewriting historical gates."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from genetics_common import (
    ContractError,
    DEFAULT_CANDIDATE_ROOT,
    PROJECT_ROOT,
    assert_candidate_root,
    atomic_write_tsv,
    clean,
    ensembl_base,
    read_table,
    sha256_file,
)


RELEASE_ID = "program-context-v2-candidate-2026-08-07"
EXPECTED_PREFLIGHT_SHA256 = (
    "b639664c3fc9199a54719a75b94bd60a4b11e57a303d91e65134894f98cf8aeb"
)
EXPECTED_HONG_GATE_SHA256 = (
    "1e9a2a48b0a6488008267c06a3cbf47861eb1b0e51c28a571f6a406fc60c2d25"
)
EXPECTED_HONG_VALIDATION_SHA256 = (
    "fe3856c9fa73166cde55a0320b8f8e1cb5d5d83897fbf0be8f98b01925cf49e8"
)
EXPECTED_SOURCE_EGENES = 6_564
EXPECTED_MAPPED_SYMBOL_ROWS = 6_583
EXPECTED_REPRESENTED_SOURCE_ENSGS = 6_562


def fail(condition: bool, message: str) -> None:
    if condition:
        raise ContractError(message)


def read_source_primary_ensgs(path: Path) -> set[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fail(
            reader.fieldnames is None
            or not {"Ensembl", "Signal"}.issubset(reader.fieldnames),
            "Broadaway snapshot schema drift",
        )
        result = {
            ensembl_base(row["Ensembl"])
            for row in reader
            if clean(row["Signal"]).startswith("1:")
        }
    fail("" in result, "blank source Ensembl identifier")
    return result


def identifier_boundary(candidate: Path) -> dict[str, object]:
    source_snapshot = (
        candidate
        / "work/input_snapshots/Liver_eQTL_Meta_Leads_ST3_20240530.tsv"
    )
    source_ensgs = read_source_primary_ensgs(source_snapshot)
    observability = read_table(candidate / "gene_observability.tsv", "tsv")
    mapped_rows = [
        row for row in observability if clean(row["source_defined_egene"]) == "true"
    ]
    assignments = [
        ensembl_base(identifier)
        for row in mapped_rows
        for identifier in clean(row["source_eqtl_ensembl"]).split(";")
        if clean(identifier)
    ]
    represented = set(assignments)
    fail(not represented.issubset(source_ensgs), "mapped identifier outside source eGenes")
    missing = sorted(source_ensgs - represented)
    return {
        "source_wide_unique_source_defined_egenes": len(source_ensgs),
        "mapped_symbol_annotation_rows": len(mapped_rows),
        "mapped_source_ensg_assignments": len(assignments),
        "represented_unique_source_ensgs": len(represented),
        "source_ensgs_not_represented_in_symbol_rows": len(missing),
        "unrepresented_source_ensgs": ";".join(missing),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--candidate-root", type=Path, default=DEFAULT_CANDIDATE_ROOT)
    args = parser.parse_args()

    fail(args.project_root.resolve() != PROJECT_ROOT.resolve(), "project-root drift")
    project_root = args.project_root.resolve()
    candidate = assert_candidate_root(project_root, args.candidate_root)

    preflight_path = candidate / "gate_status.tsv"
    hong_gate_path = candidate / "hong_source_gate_status.tsv"
    hong_validation_path = candidate / "hong_validation_report.tsv"
    for path in (preflight_path, hong_gate_path, hong_validation_path):
        fail(not path.is_file(), f"required terminal source missing: {path}")
    fail(
        sha256_file(preflight_path) != EXPECTED_PREFLIGHT_SHA256,
        "historical preflight gate drift",
    )
    fail(
        sha256_file(hong_gate_path) != EXPECTED_HONG_GATE_SHA256,
        "terminal Hong gate drift",
    )
    fail(
        sha256_file(hong_validation_path) != EXPECTED_HONG_VALIDATION_SHA256,
        "terminal Hong validation drift",
    )

    preflight = {
        row["gate_id"]: row for row in read_table(preflight_path, "tsv")
    }
    fail(
        preflight.get("GEN03_HONG_SOURCE_AUDIT", {}).get("status")
        != "not_started_outside_preflight_scope",
        "historical GEN03 preflight row no longer matches the superseded state",
    )
    fail(
        preflight.get("PREFLIGHT_OVERALL", {}).get("status")
        != "preflight_complete_claim_integration_blocked",
        "historical PREFLIGHT_OVERALL row no longer matches the superseded state",
    )

    hong_gate = {
        row["gate_id"]: row for row in read_table(hong_gate_path, "tsv")
    }
    fail(
        hong_gate.get("GEN03_OVERALL", {}).get("status")
        != "coverage_limited_terminal",
        "GEN03 terminal status is not coverage_limited_terminal",
    )
    fail(
        hong_gate.get("GEN03_CONTEXT_RESCUE_AUTHORIZATION", {}).get("status")
        != "prohibited",
        "context-rescue prohibition drift",
    )
    fail(
        hong_gate.get("GEN03_NEGATIVE_CLAIM_AUTHORIZATION", {}).get("status")
        != "prohibited",
        "negative-claim prohibition drift",
    )

    validation = {
        row["check_id"]: row for row in read_table(hong_validation_path, "tsv")
    }
    fail(
        validation.get("terminal_authorization", {}).get("status") != "pass"
        or validation.get("terminal_authorization", {}).get("observed")
        != "coverage_limited_terminal",
        "independent Hong terminal validation is absent or failed",
    )

    boundary = identifier_boundary(candidate)
    fail(
        boundary["source_wide_unique_source_defined_egenes"]
        != EXPECTED_SOURCE_EGENES,
        "source-wide eGene count drift",
    )
    fail(
        boundary["mapped_symbol_annotation_rows"] != EXPECTED_MAPPED_SYMBOL_ROWS,
        "mapped symbol-row count drift",
    )
    fail(
        boundary["represented_unique_source_ensgs"]
        != EXPECTED_REPRESENTED_SOURCE_ENSGS,
        "represented source-Ensembl count drift",
    )

    closure_row = {
        "release_id": RELEASE_ID,
        "closure_id": "GEN_TERMINAL_CLOSURE_V1",
        "terminal_status": "coverage_limited_terminal",
        "downstream_selection_authority": "terminal_closure",
        "authoritative_gate_path": "hong_source_gate_status.tsv",
        "authoritative_gate_sha256": sha256_file(hong_gate_path),
        "independent_validation_path": "hong_validation_report.tsv",
        "independent_validation_sha256": sha256_file(hong_validation_path),
        "historical_preflight_path": "gate_status.tsv",
        "historical_preflight_sha256": sha256_file(preflight_path),
        "superseded_preflight_gate_ids": "GEN03_HONG_SOURCE_AUDIT;PREFLIGHT_OVERALL",
        "supersession_scope": "downstream_gate_selection_only",
        "historical_preflight_preserved": "true",
        "context_rescue_authorized": "false",
        "negative_claim_authorized": "false",
        "allowed_downstream_use": "phenotype_provenance;positive_only_broadaway_summary;all_joint_34_of_447_descriptive_interface",
        "prohibited_downstream_use": "context_rescue;rescued_fraction;context_enrichment;source_negative;powered_genetic_null",
        **boundary,
        "caption_identifier_rule": "6564=source-wide unique source-defined eGenes;6583=mapped symbol annotation rows;6562=represented unique source ENSGs;never substitute one unit for another",
        "canonical_promotion_status": "not_promoted",
    }
    closure_path = candidate / "terminal_closure.tsv"
    atomic_write_tsv(closure_path, [closure_row], list(closure_row))

    script_dir = Path(__file__).resolve().parent
    code_paths = [
        script_dir / "09_close_genetics_terminal_gate.py",
        script_dir / "10_validate_genetics_terminal_closure.py",
        script_dir / "run_genetics_terminal_closure.sbatch",
    ]
    for path in code_paths:
        fail(not path.is_file(), f"required reviewed closure code missing: {path}")
    sources = [
        (preflight_path, "historical_preflight_preserved"),
        (hong_gate_path, "terminal_gate_authority"),
        (hong_validation_path, "terminal_independent_validation"),
        (candidate / "source_eqtl_reproduction.tsv", "source_count_contract"),
        (candidate / "gene_observability.tsv", "mapped_symbol_contract"),
        (candidate / "gene_mapping_audit.tsv", "mapping_count_audit"),
        (
            candidate
            / "work/input_snapshots/Liver_eQTL_Meta_Leads_ST3_20240530.tsv",
            "source_ensg_identity_contract",
        ),
        (closure_path, "terminal_closure_output"),
        (code_paths[0], "terminal_closure_producer"),
        (code_paths[1], "independent_closure_validator"),
        (code_paths[2], "review_only_slurm_wrapper"),
    ]
    manifest_rows = []
    for path, role in sources:
        fail(not path.is_file(), f"closure manifest source missing: {path}")
        manifest_rows.append(
            {
                "relative_path": str(path.resolve().relative_to(project_root)),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "role": role,
                "canonical_promotion_status": "not_promoted",
            }
        )
    manifest_path = candidate / "terminal_closure_manifest.tsv"
    atomic_write_tsv(manifest_path, manifest_rows, list(manifest_rows[0]))

    fail(
        sha256_file(preflight_path) != EXPECTED_PREFLIGHT_SHA256,
        "historical preflight changed while closure was created",
    )
    print(
        "PASS: GEN terminal closure created; stale preflight preserved and "
        "superseded only for downstream gate selection"
    )


if __name__ == "__main__":
    main()

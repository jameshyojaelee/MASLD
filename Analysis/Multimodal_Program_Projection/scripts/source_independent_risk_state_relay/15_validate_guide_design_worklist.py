#!/usr/bin/env python3
"""Validate separated promoter-perturbation and exact-edit worklists."""

from __future__ import annotations

import json

from relay_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def main() -> None:
    seal_path = CANDIDATE_ROOT / "GUIDE_DESIGN_WORKLIST_SEALED.json"
    if not seal_path.is_file():
        raise RuntimeError("Missing guide-design seal")
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("status") != "separate_stage_a_promoter_and_stage_b_exact_designs_pending":
        raise RuntimeError("Invalid guide-design release status")
    if seal.get("experimental_targets_frozen") is not False:
        raise RuntimeError("Guide-design worklist froze targets")
    for name, expected in seal["output_sha256"].items():
        path = CANDIDATE_ROOT / name
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Guide-design output hash mismatch: {name}")

    stage_a = read_tsv(CANDIDATE_ROOT / "stage_a_promoter_design_worklist.tsv")
    stage_b = read_tsv(CANDIDATE_ROOT / "stage_b_exact_edit_worklist.tsv")
    contracts = read_tsv(CANDIDATE_ROOT / "external_design_contract.tsv")
    manifests = read_tsv(CANDIDATE_ROOT / "guide_design_input_manifest.tsv")
    if len(stage_a) != int(seal["n_stage_a_promoter_requests"]):
        raise RuntimeError("Stage-A promoter-request count mismatch")
    if len(stage_b) != int(seal["n_stage_b_accessible_variant_designs"]):
        raise RuntimeError("Stage-B exact-edit row count mismatch")
    if len(stage_a) != int(seal["n_routed_pairs"]):
        raise RuntimeError("Stage A must contain one promoter request per routed pair")
    if len({row["promoter_request_uid"] for row in stage_a}) != len(stage_a):
        raise RuntimeError("Duplicate Stage-A promoter request UID")
    if len({row["orientation_uid"] for row in stage_a}) != len(stage_a):
        raise RuntimeError("Duplicate Stage-A routed orientation")
    if len({row["pair_family_uid"] for row in stage_a}) != int(
        seal["n_routed_pair_families"]
    ):
        raise RuntimeError("Stage-A routed pair-family count mismatch")
    if len({row["coarse_locus_uid"] for row in stage_a}) != int(
        seal["n_routed_physical_loci"]
    ):
        raise RuntimeError("Stage-A routed physical-locus count mismatch")
    if len({row["design_uid"] for row in stage_b}) != len(stage_b):
        raise RuntimeError("Duplicate Stage-B exact-edit design UID")
    if len(contracts) != len(stage_a) + len(stage_b):
        raise RuntimeError("External design-contract cardinality mismatch")
    if {row["design_arm"] for row in contracts} != {
        "stage_a_promoter_gene_perturbation",
        "stage_b_exact_edit",
    }:
        raise RuntimeError("External design-arm universe drift")
    for row in stage_a:
        expected_mode = {
            "risk_increases_expression": "CRISPRa",
            "risk_decreases_expression": "CRISPRi",
        }.get(row["oriented_risk_effect"])
        if row["stage_a_perturbation_mode"] != expected_mode:
            raise RuntimeError(
                f"Stage-A risk-direction/mode mismatch: {row['orientation_uid']}"
            )
        expected_architecture = (
            row["dominant_celltype"] == "Hepatocytes"
            and row["accessibility_evidence_class"]
            == "two_cohort_hepatocyte_accessibility"
        )
        if yes(row["primary_architecture_eligible_before_guide_gate"]) != expected_architecture:
            raise RuntimeError(
                f"Stage-A architecture gate mismatch: {row['orientation_uid']}"
            )
        if row["target_freeze_status"] != (
            "prohibited_until_external_design_and_protocol_gate"
        ):
            raise RuntimeError("A Stage-A promoter request opened target selection")

    if any(len(row["sequence_401bp"]) != 401 for row in stage_b):
        raise RuntimeError("Sequence-window length drift")
    if any(
        row["reference_allele_forward"] != row["sequence_401bp"][200]
        for row in stage_b
    ):
        raise RuntimeError("Reference allele is not centered in sequence window")
    if any(
        row["target_freeze_status"]
        != "prohibited_until_external_design_and_protocol_gate"
        for row in stage_b
    ):
        raise RuntimeError("A Stage-B exact-edit row opened target selection")
    exact_count = sum(yes(row["exact_edit_exportable"]) for row in stage_b)
    if exact_count != int(seal["n_exact_edit_exportable"]):
        raise RuntimeError("Exact-edit export count mismatch")
    for row in stage_b:
        if yes(row["exact_edit_exportable"]):
            if not row["primedesign_input"] or not row["alternate_allele_forward"]:
                raise RuntimeError(f"Incomplete exact-edit export: {row['design_uid']}")
        elif row["primedesign_input"]:
            raise RuntimeError(f"Unresolved edit has PrimeDesign input: {row['design_uid']}")

    for row in manifests:
        source = row["source_path"]
        path = (
            PROJECT_ROOT / source
            if not source.startswith("/")
            else __import__("pathlib").Path(source)
        )
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]):
            raise RuntimeError(f"Guide-design source drift: {path}")
        if sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Guide-design source hash drift: {path}")

    print(
        "GUIDE_DESIGN_WORKLIST_VALIDATION_PASS "
        f"stage_a_pairs={len(stage_a)} stage_b_variants={len(stage_b)} "
        f"exact_edit_exportable={exact_count} "
        "external_designs_pending=true targets_frozen=false"
    )


if __name__ == "__main__":
    main()

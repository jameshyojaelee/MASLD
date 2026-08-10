#!/usr/bin/env python3
"""Independently validate the full-goal completion audit and missing-data contract."""

from __future__ import annotations

import json

from relay_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


def main() -> None:
    seal = json.loads((CANDIDATE_ROOT / "GOAL_COMPLETION_AUDIT.json").read_text(encoding="utf-8"))
    if seal.get("status") != "goal_not_achieved_new_data_required":
        raise RuntimeError("Invalid full-goal audit status")
    false_fields = [
        "paired_human_reversal_achieved", "orthogonal_translation_achieved",
        "regulatory_perturbation_achieved", "full_risk_to_state_architecture_achieved",
        "paper_promotion_authorized", "new_scientific_outcomes_inspected",
    ]
    if any(seal.get(field) is not False for field in false_fields):
        raise RuntimeError("Full-goal audit overstates completion or outcome access")
    if seal.get("public_data_escalation_closed") is not True or seal.get("retrospective_terminal_evidence_audited") is not True:
        raise RuntimeError("Full-goal audit omits terminal public-data boundary")
    for relative, expected in seal["output_sha256"].items():
        path = CANDIDATE_ROOT / relative
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Full-goal audit output drift: {relative}")

    source_rows = read_tsv(CANDIDATE_ROOT / "goal_audit_source_manifest.tsv")
    expected_source_roles = {
        "plan43_release", "plan43_promotion", "plan43_human", "plan43_protein",
        "plan43_spatial", "plan43_cropseq", "plan44_terminal_release",
        "plan44_promotion", "plan45_handoff_seal", "plan45_response_marker",
        "producer", "validator",
    }
    if {row["role"] for row in source_rows} != expected_source_roles:
        raise RuntimeError("Full-goal source universe drift")
    sources = {}
    for row in source_rows:
        path = PROJECT_ROOT / row["source_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Full-goal source drift: {path}")
        sources[row["role"]] = path

    promotion = {row["component_id"]: row for row in read_tsv(sources["plan43_promotion"])}
    required_failed = {
        "human_histologic_reversal", "regulatory_risk_cropseq_bridge",
        "protein_translation", "spatial_lineage_remodeling",
        "complete_risk_to_state_architecture",
    }
    if not required_failed <= set(promotion) or any(promotion[key]["gate_pass"] != "False" for key in required_failed):
        raise RuntimeError("Plan 43 terminal failures were not preserved")
    human = read_tsv(sources["plan43_human"])
    primary = {
        row["dataset_id"]: row
        for row in human
        if row["score_scheme"] == "weighted__primary" and row["inferential_status"] == "primary"
    }
    if not {"GSE106737", "GSE83452", "GSE48452"} <= set(primary):
        raise RuntimeError("Plan 43 paired-human primary universe drift")
    if float(primary["GSE83452"]["p"]) < 0.95 or primary["GSE48452"]["direction_agrees"] != "FALSE":
        raise RuntimeError("Plan 43 paired-human failure semantics drift")
    protein = {row["dataset_id"]: row for row in read_tsv(sources["plan43_protein"])}
    if float(protein["PXD052787"]["q"]) < 0.70:
        raise RuntimeError("Plan 43 secretome null semantics drift")
    crop = read_tsv(sources["plan43_cropseq"])
    if len(crop) != 1 or crop[0]["status"] != "skipped_insufficient_oriented_loci" or crop[0]["inference_authorized"] != "false":
        raise RuntimeError("Plan 43 CROP-seq skip semantics drift")
    spatial = read_tsv(sources["plan43_spatial"])
    if not any(row["status"] == "skipped_missing_prespecified_gene_level_spatial_variability_and_coverage_covariates" for row in spatial):
        raise RuntimeError("Plan 43 spatial source-gate failure drift")
    handoff = json.loads(sources["plan45_handoff_seal"].read_text(encoding="utf-8"))
    if handoff.get("experimental_targets_frozen") is not False or handoff.get("scientific_outcomes_inspected") is not False:
        raise RuntimeError("Plan 45 handoff improperly treated as biological evidence")
    response = json.loads(sources["plan45_response_marker"].read_text(encoding="utf-8"))
    if response.get("target_identity_present") is not False or response.get("scientific_outcomes_present") is not False:
        raise RuntimeError("Plan 45 response scaffold improperly treated as biological evidence")

    matrix = read_tsv(CANDIDATE_ROOT / "ultimate_goal_requirement_matrix.tsv")
    if [row["requirement_id"] for row in matrix] != [f"UG{i:02d}" for i in range(1, 12)]:
        raise RuntimeError("Ultimate-goal requirement universe drift")
    allowed = {"partial", "failed", "missing", "not_achieved", "achieved_design_only"}
    if not {row["current_verdict"] for row in matrix} <= allowed:
        raise RuntimeError("Ultimate-goal verdict vocabulary drift")
    if any(not row["decisive_missing_evidence"] or not row["completion_test"] for row in matrix):
        raise RuntimeError("Ultimate-goal audit lacks a decisive completion test")
    if any(row["current_verdict"] == "achieved" for row in matrix):
        raise RuntimeError("Ultimate-goal audit improperly claims biological completion")
    if seal.get("n_requirements") != len(matrix) or seal.get("n_achieved_biological_requirements") != 0:
        raise RuntimeError("Ultimate-goal completion count drift")

    package = read_tsv(CANDIDATE_ROOT / "minimal_new_data_package.tsv")
    if [row["package_id"] for row in package] != [f"NDP{i:02d}" for i in range(1, 9)]:
        raise RuntimeError("New-data package universe drift")
    if any(not row["freeze_before_outcome"] or not row["pass_gate"] or not row["failure_interpretation"] for row in package):
        raise RuntimeError("New-data package lacks a prospective gate")
    prohibited_units = {"cell", "spot", "section", "guide", "well", "organoid"}
    if any(row["biological_unit"].strip().lower() in prohibited_units for row in package):
        raise RuntimeError("New-data package uses a technical unit as the biological unit")

    claims = read_tsv(CANDIDATE_ROOT / "claim_promotion_logic.tsv")
    if [row["claim_id"] for row in claims] != [f"CL{i:02d}" for i in range(1, 7)]:
        raise RuntimeError("Claim-promotion universe drift")
    architecture = next(row for row in claims if row["claim_id"] == "CL04")
    required_architecture = {f"UG{i:02d}" for i in range(2, 10)}
    if set(architecture["required_components"].split(";")) != required_architecture:
        raise RuntimeError("Full architecture gate omits a load-bearing requirement")
    if seal.get("n_new_data_components") != len(package) or seal.get("n_claim_rules") != len(claims):
        raise RuntimeError("Full-goal package count drift")
    print(
        "GOAL_COMPLETION_AUDIT_VALIDATION_PASS requirements=11 achieved_biological=0 "
        "new_data_components=8 full_claim=false promotion=false"
    )


if __name__ == "__main__":
    main()

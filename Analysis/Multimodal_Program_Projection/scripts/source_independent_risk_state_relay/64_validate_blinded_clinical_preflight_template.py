#!/usr/bin/env python3
"""Independently rederive the frozen outcome-free Plan 46 clinical template."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[4]
ROOT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/source-independent-risk-state-relay-clinical-preflight-template-v1-2026-08-10"
EXPECTED_AXIS_SHA256 = "03242574b4a4aee810000cd2304abb0ba9f3c960ee7b2e0f584c3e3c8a180e33"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    seal_path = ROOT / "CLINICAL_PREFLIGHT_TEMPLATE_SEALED.json"
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    expected = {
        "status": "sealed_outcome_free_two_cohort_clinical_preflight_template",
        "n_frozen_programs": 117,
        "n_cohort_roles": 2,
        "n_responder_endpoints": 3,
        "n_estimands": 3,
        "n_orthogonal_assays": 3,
        "bundle_read_only": True,
        "real_partner_response_present": False,
        "participant_rows_present": False,
        "new_clinical_outcomes_present": False,
        "new_molecular_outcomes_present": False,
        "cohort_selected": False,
        "outcome_key_release_authorized": False,
        "analysis_authorized": False,
        "paper_promotion_authorized": False,
        "state_axis_sha256": EXPECTED_AXIS_SHA256,
    }
    for field, value in expected.items():
        if seal.get(field) != value:
            raise RuntimeError(f"Clinical preflight seal drift: {field}")
    if (ROOT / "CLINICAL_PREFLIGHT_TEMPLATE_SEAL_SHA256.txt").read_text(encoding="utf-8").strip() != sha256_file(seal_path):
        raise RuntimeError("Clinical preflight seal sidecar mismatch")
    for relative, digest in seal["output_sha256"].items():
        path = ROOT / relative
        if not path.is_file() or sha256_file(path) != digest:
            raise RuntimeError(f"Clinical preflight output drift: {relative}")
    source = read_tsv(ROOT / "source_manifest.tsv")
    if {row["role"] for row in source} != {"plan43", "plan46", "plan46b", "plan46c", "state_axis", "state_spec", "producer", "validator"}:
        raise RuntimeError("Clinical preflight source universe drift")
    for row in source:
        path = PROJECT_ROOT / row["path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Clinical preflight source drift: {path}")
    binding = read_tsv(ROOT / "state_axis_binding.tsv")
    if len(binding) != 1 or binding[0]["state_axis_sha256"] != EXPECTED_AXIS_SHA256 or binding[0]["n_programs"] != "117":
        raise RuntimeError("Clinical preflight state-axis binding drift")
    axis = read_tsv(PROJECT_ROOT / binding[0]["state_axis_path"])
    if len(axis) != 117 or abs(sum(abs(float(row["primary_loading"])) for row in axis) - 1.0) > 1e-10:
        raise RuntimeError("Clinical preflight state-axis rederivation failed")
    endpoints = read_tsv(ROOT / "histologic_responder_definitions.tsv")
    if [row["endpoint_id"] for row in endpoints] != [
        "primary_histologic_improvement", "mash_resolution_no_fibrosis_worsening",
        "fibrosis_improvement",
    ]:
        raise RuntimeError("Clinical responder endpoint drift")
    primary = endpoints[0]
    for phrase in (
        "NAS_decrease_ge_2", "ballooning_or_lobular_inflammation_improvement",
        "no_fibrosis_worsening",
    ):
        if phrase not in primary["improver_definition"]:
            raise RuntimeError(f"Clinical primary responder rule lost: {phrase}")
    if primary["substitution_allowed"] != "false":
        raise RuntimeError("Clinical responder substitution became permitted")
    estimands = read_tsv(ROOT / "prospective_human_estimands.tsv")
    if [row["estimand_id"] for row in estimands] != ["HUM_PRIMARY", "HUM_COVARIATE_SENSITIVITY", "HUM_TWO_COHORT_META"]:
        raise RuntimeError("Clinical estimand universe drift")
    if estimands[0]["expected_direction"] != "negative" or "cannot_rescue" not in estimands[2]["promotion_role"]:
        raise RuntimeError("Clinical direction or replication rule drift")
    protocols = read_tsv(ROOT / "cohort_protocol_template.tsv")
    if [row["portfolio_role"] for row in protocols] != ["cohort_A", "cohort_B"]:
        raise RuntimeError("Clinical cohort protocol role drift")
    if any(row["clinical_response_candidate_id"] or row["cohort_uid"] or row["protocol_sha256"] or row["protocol_locked_utc"] or row["outcome_key_released_utc"] for row in protocols):
        raise RuntimeError("Clinical cohort template was prefilled")
    if any(row["molecular_outcome_accessed"] != "false" for row in protocols):
        raise RuntimeError("Clinical cohort template crossed the outcome firewall")
    power = read_tsv(ROOT / "blinded_power_audit_template.tsv")
    if len(power) != 2 or any(row["power_gate"] != "pending" or row["state_axis_or_responder_linked_molecular_outcome_inspected"] != "false" for row in power):
        raise RuntimeError("Clinical blinded power template drift")
    orthogonal = read_tsv(ROOT / "orthogonal_translation_protocol.tsv")
    if {row["assay_id"] for row in orthogonal} != {"tissue_proteomics", "liver_linked_secreted_proteomics", "physical_niche"}:
        raise RuntimeError("Clinical orthogonal assay universe drift")
    if any(row["substitution_allowed"] != "false" for row in orthogonal):
        raise RuntimeError("Clinical orthogonal substitution became permitted")
    manifest = read_tsv(ROOT / "file_manifest.tsv")
    for row in manifest:
        path = ROOT / row["relative_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Clinical preflight manifest drift: {path}")
    for path in ROOT.rglob("*"):
        if path.is_file() and path.stat().st_mode & 0o222:
            raise RuntimeError(f"Clinical preflight payload remains writable: {path}")
    if ROOT.stat().st_mode & 0o222:
        raise RuntimeError("Clinical preflight root remains writable")
    print(
        "CLINICAL_PREFLIGHT_TEMPLATE_VALIDATION_PASS "
        "programs=117 endpoints=3 estimands=3 orthogonal_assays=3 "
        "cohorts_selected=false outcomes=false"
    )


if __name__ == "__main__":
    main()

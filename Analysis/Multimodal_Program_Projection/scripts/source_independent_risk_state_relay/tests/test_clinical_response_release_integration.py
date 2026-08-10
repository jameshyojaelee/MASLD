from __future__ import annotations

import csv
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import uuid
import json
from pathlib import Path


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = SCRIPT_ROOT.parents[3]
CANDIDATES = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates"


def read_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def write_rows(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_writable(root: Path) -> None:
    if not root.exists():
        return
    root.chmod(0o755)
    for path in root.rglob("*"):
        path.chmod(0o755 if path.is_dir() else 0o644)


class ClinicalResponseReleaseIntegrationTests(unittest.TestCase):
    def test_complete_metadata_only_response_seals_and_revalidates(self) -> None:
        token = uuid.uuid4().hex[:12]
        draft_id = f"source-independent-risk-state-relay-clinical-test-draft-{token}"
        release_id = f"source-independent-risk-state-relay-clinical-test-release-{token}"
        draft = CANDIDATES / draft_id
        release = CANDIDATES / release_id
        try:
            subprocess.run([
                sys.executable, str(SCRIPT_ROOT / "57_build_clinical_response_scaffold.py"),
                "--route-id", "nash_crn_as116", "--candidate-id", draft_id,
            ], check=True, capture_output=True, text=True)
            evidence = draft / "evidence" / "custodian_attestation.txt"
            evidence.parent.mkdir()
            evidence.write_text("Synthetic metadata-only test attestation. No participant rows or outcomes.\n", encoding="utf-8")
            relative = "evidence/custodian_attestation.txt"
            digest = hashlib.sha256(evidence.read_bytes()).hexdigest()
            write_rows(draft / "evidence_file_register.tsv", [{
                "evidence_id": "E01", "relative_path": relative,
                "description": "Synthetic metadata-only integration fixture",
                "sha256": digest, "size_bytes": str(evidence.stat().st_size),
            }], ["evidence_id", "relative_path", "description", "sha256", "size_bytes"])

            rows, fields = read_rows(draft / "clinical_cohort_feasibility_questions.tsv")
            for row in rows:
                row.update({
                    "cohort_role_candidate": "candidate_cohort_A", "answer": "yes",
                    "evidence_path": relative, "owner": "synthetic_custodian",
                    "signed_utc": "2026-08-10T12:00:00+00:00",
                })
            write_rows(draft / "clinical_cohort_feasibility_questions.tsv", rows, fields)

            rows, fields = read_rows(draft / "assay_availability_matrix.tsv")
            for row in rows:
                row.update({
                    "available": "yes", "n_baseline": "20", "n_followup": "20",
                    "n_authoritative_pairs": "20", "evidence_path": relative,
                })
            write_rows(draft / "assay_availability_matrix.tsv", rows, fields)

            rows, fields = read_rows(draft / "clinical_role_firewall.tsv")
            for row in rows:
                if row["role"] == "clinical_custodian":
                    row.update({
                        "assigned_person": "synthetic_custodian",
                        "accepted_utc": "2026-08-10T12:00:00+00:00",
                    })
            write_rows(draft / "clinical_role_firewall.tsv", rows, fields)

            rows, fields = read_rows(draft / "route_specific_questions.tsv")
            for row in rows:
                row.update({
                    "answer": "yes", "evidence_path": relative,
                    "owner": "synthetic_custodian",
                    "signed_utc": "2026-08-10T12:00:00+00:00",
                })
            write_rows(draft / "route_specific_questions.tsv", rows, fields)

            rows, fields = read_rows(draft / "public_inventory_audit.tsv")
            for row in rows:
                row["custodian_confirmation"] = "synthetic_confirmed_for_test_only"
            write_rows(draft / "public_inventory_audit.tsv", rows, fields)

            rows, fields = read_rows(draft / "blinded_cohort_census.tsv")
            values = {
                "unique_participants": "30",
                "authoritative_baseline_followup_pairs": "20",
                "complete_histology_pairs": "20",
                "paired_genomewide_rna_pairs": "20",
                "histologic_improvers": "8",
                "histologic_non_improvers": "12",
                "five_assay_complete_pairs": "20",
            }
            for row in rows:
                row.update({"value": values[row["metric"]], "evidence_path": relative})
            write_rows(draft / "blinded_cohort_census.tsv", rows, fields)

            rows, fields = read_rows(draft / "clinical_response_identity.tsv")
            rows[0].update({
                "cohort_uid": "synthetic_cohort", "institution_uid": "synthetic_inst_a",
                "trial_uid": "synthetic_trial_a", "named_custodian": "synthetic_custodian",
                "governance_route": "synthetic_governance", "credible_completion_utc": "2026-12-01T00:00:00+00:00",
                "prepared_by": "synthetic_custodian", "submitted_utc": "2026-08-10T12:00:00+00:00",
                "response_status": "submitted",
            })
            write_rows(draft / "clinical_response_identity.tsv", rows, fields)

            subprocess.run([
                sys.executable, str(SCRIPT_ROOT / "59_seal_and_adjudicate_clinical_response.py"),
                "--draft-root", str(draft), "--candidate-id", release_id,
            ], check=True, capture_output=True, text=True)
            result = subprocess.run([
                sys.executable, str(SCRIPT_ROOT / "60_validate_clinical_response_release.py"),
                "--root", str(release),
            ], check=False, capture_output=True, text=True)
            if result.returncode:
                self.fail(f"Clinical release validator failed:\n{result.stdout}\n{result.stderr}")
            self.assertIn("verdict=pass_to_blinded_census", result.stdout)
            self.assertIn("full_nested=true", result.stdout)
            self.assertIn("outcomes=false", result.stdout)
        finally:
            for root in (release, draft):
                make_writable(root)
                shutil.rmtree(root, ignore_errors=True)

    def test_two_sealed_metadata_responses_form_outcome_locked_portfolio(self) -> None:
        token = uuid.uuid4().hex[:12]
        root_a = CANDIDATES / f"source-independent-risk-state-relay-clinical-test-a-{token}"
        root_b = CANDIDATES / f"source-independent-risk-state-relay-clinical-test-b-{token}"
        portfolio_id = f"source-independent-risk-state-relay-clinical-test-portfolio-{token}"
        portfolio = CANDIDATES / portfolio_id
        try:
            base = {
                "status": "sealed_metadata_only_clinical_response_adjudication",
                "output_sha256": {},
                "participant_rows_present": False,
                "clinical_outcomes_present": False,
                "molecular_outcomes_present": False,
                "participant_data_request_authorized": False,
                "molecular_outcome_access_authorized": False,
                "analysis_authorized": False,
                "paper_promotion_authorized": False,
                "pass_to_blinded_census": True,
                "feasibility_verdict": "pass_to_blinded_census",
            }
            fixtures = (
                (root_a, {
                    **base, "route_id": "nash_crn_as116",
                    "candidate_cohort_role": "candidate_cohort_A",
                    "independence_group": "NASH_CRN", "institution_uid": "inst_a",
                    "trial_uid": "trial_a", "full_goal_nested_candidate": False,
                }),
                (root_b, {
                    **base, "route_id": "maestro_nash",
                    "candidate_cohort_role": "candidate_cohort_B",
                    "independence_group": "MADRIGAL_MAESTRO", "institution_uid": "inst_b",
                    "trial_uid": "trial_b", "full_goal_nested_candidate": True,
                }),
            )
            for root, seal in fixtures:
                root.mkdir()
                (root / "CLINICAL_RESPONSE_SEALED.json").write_text(
                    json.dumps(seal, indent=2, sort_keys=True) + "\n", encoding="utf-8"
                )
            with tempfile.TemporaryDirectory() as temporary:
                audit_root = Path(temporary)
                evidence = audit_root / "evidence" / "independence_attestation.txt"
                evidence.parent.mkdir()
                evidence.write_text("Synthetic independence evidence; no participant rows.\n", encoding="utf-8")
                relative = "evidence/independence_attestation.txt"
                write_rows(audit_root / "evidence_file_register.tsv", [{
                    "evidence_id": "E01", "relative_path": relative,
                    "description": "Synthetic independence fixture",
                    "sha256": hashlib.sha256(evidence.read_bytes()).hexdigest(),
                    "size_bytes": str(evidence.stat().st_size),
                }], ["evidence_id", "relative_path", "description", "sha256", "size_bytes"])
                fields = [
                    "cohort_a_candidate_id", "cohort_b_candidate_id",
                    "institutionally_independent", "trial_independent",
                    "participant_independent", "overlap_resolved",
                    "cohort_b_locked_before_cohort_a_outcomes", "evidence_path",
                    "owner", "signed_utc",
                ]
                write_rows(audit_root / "pairwise_independence_audit.tsv", [{
                    "cohort_a_candidate_id": root_a.name,
                    "cohort_b_candidate_id": root_b.name,
                    "institutionally_independent": "yes", "trial_independent": "yes",
                    "participant_independent": "yes", "overlap_resolved": "yes",
                    "cohort_b_locked_before_cohort_a_outcomes": "yes",
                    "evidence_path": relative, "owner": "synthetic_validator",
                    "signed_utc": "2026-08-10T12:00:00+00:00",
                }], fields)
                build = subprocess.run([
                    sys.executable, str(SCRIPT_ROOT / "61_adjudicate_clinical_response_portfolio.py"),
                    "--cohort-a-root", str(root_a), "--cohort-b-root", str(root_b),
                    "--audit-root", str(audit_root), "--candidate-id", portfolio_id,
                ], check=False, capture_output=True, text=True)
                if build.returncode:
                    self.fail(f"Clinical portfolio builder failed:\n{build.stdout}\n{build.stderr}")
                result = subprocess.run([
                    sys.executable, str(SCRIPT_ROOT / "62_validate_clinical_response_portfolio.py"),
                    "--root", str(portfolio),
                ], check=False, capture_output=True, text=True)
                if result.returncode:
                    self.fail(f"Clinical portfolio validator failed:\n{result.stdout}\n{result.stderr}")
                self.assertIn("verdict=pass_two_cohort_blinded_preflight", result.stdout)
                self.assertIn("orthogonal_nested=true", result.stdout)
                self.assertIn("outcomes=false", result.stdout)
        finally:
            for root in (portfolio, root_b, root_a):
                make_writable(root)
                shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()

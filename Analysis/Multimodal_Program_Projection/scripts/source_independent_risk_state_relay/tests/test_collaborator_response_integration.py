from __future__ import annotations

import csv
import hashlib
import os
import shutil
import subprocess
import sys
import unittest
import uuid
from pathlib import Path


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = SCRIPT_ROOT.parents[3]
CANDIDATES = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates"
HANDOFF_ROOT = CANDIDATES / "source-independent-risk-state-relay-stage-a-collaborator-handoff-v3-2026-08-10"


def update_tsv(path: Path, transform) -> None:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = [transform(dict(row)) for row in reader]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


@unittest.skipUnless(HANDOFF_ROOT.is_dir(), "sealed collaborator handoff unavailable")
class CollaboratorResponseIntegrationTests(unittest.TestCase):
    def test_all_yes_response_seals_and_rederives(self) -> None:
        token = uuid.uuid4().hex
        draft_id = f"source-independent-risk-state-relay-fixture-response-draft-{token}"
        release_id = f"source-independent-risk-state-relay-fixture-response-release-{token}"
        draft_root = CANDIDATES / draft_id
        release_root = CANDIDATES / release_id
        common = os.environ.copy()
        common["PLAN45_COLLABORATOR_HANDOFF_ROOT"] = str(HANDOFF_ROOT)
        try:
            build_env = {**common, "PLAN45_CANDIDATE_ID": draft_id}
            subprocess.run([sys.executable, str(SCRIPT_ROOT / "45_build_collaborator_response_scaffold.py")], env=build_env, check=True, capture_output=True, text=True)
            subprocess.run([sys.executable, str(SCRIPT_ROOT / "46_validate_collaborator_response_scaffold.py")], env=build_env, check=True, capture_output=True, text=True)

            evidence = draft_root / "evidence/platform_attestation.txt"
            evidence.parent.mkdir()
            evidence.write_text("synthetic integration fixture; not scientific evidence\n", encoding="utf-8")
            digest = hashlib.sha256(evidence.read_bytes()).hexdigest()
            evidence_size = evidence.stat().st_size
            update_tsv(
                draft_root / "platform_feasibility_questions.tsv",
                lambda row: {
                    **row, "answer": "yes", "evidence_path": "evidence/platform_attestation.txt",
                    "owner": "fixture_owner", "signed_utc": "2026-08-10T12:00:00+00:00",
                },
            )
            role_people = {
                "platform_lead": "fixture_platform", "wetlab_lead": "fixture_wetlab",
                "data_manager": "fixture_data", "qc_reviewer_1": "fixture_qc1",
                "qc_reviewer_2": "fixture_qc2", "analysis_lead": "fixture_analysis",
                "independent_validator": "fixture_validator",
            }
            update_tsv(
                draft_root / "role_and_blinding_firewall.tsv",
                lambda row: {
                    **row, "assigned_person": role_people[row["role"]],
                    "accepted_utc": "2026-08-10T12:00:00+00:00",
                },
            )
            update_tsv(
                draft_root / "collaborator_deliverable_register.tsv",
                lambda row: {**row, "status": "complete" if row["deliverable_id"] == "D01" else "not_started"},
            )
            update_tsv(
                draft_root / "response_identity.tsv",
                lambda row: {
                    **row, "organization": "fixture_org", "prepared_by": "fixture_preparer",
                    "created_utc": "2026-08-10T12:00:00+00:00", "response_status": "submitted",
                },
            )
            register = draft_root / "evidence_file_register.tsv"
            with register.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["evidence_id", "relative_path", "description", "sha256", "size_bytes"],
                    delimiter="\t", lineterminator="\n",
                )
                writer.writeheader()
                writer.writerow({
                    "evidence_id": "EV01", "relative_path": "evidence/platform_attestation.txt",
                    "description": "synthetic integration fixture", "sha256": digest,
                    "size_bytes": evidence_size,
                })

            seal_env = {
                **common,
                "PLAN45_CANDIDATE_ID": release_id,
                "PLAN45_COLLABORATOR_RESPONSE_DRAFT_ROOT": str(draft_root),
            }
            subprocess.run([sys.executable, str(SCRIPT_ROOT / "47_seal_and_adjudicate_collaborator_response.py")], env=seal_env, check=True, capture_output=True, text=True)
            result = subprocess.run([sys.executable, str(SCRIPT_ROOT / "48_validate_collaborator_response_release.py")], env=seal_env, check=True, capture_output=True, text=True)
            self.assertIn("target_disclosure=true", result.stdout)
            self.assertIn("target_frozen=false", result.stdout)
        finally:
            if draft_root.is_dir():
                shutil.rmtree(draft_root)
            if release_root.is_dir():
                shutil.rmtree(release_root)


if __name__ == "__main__":
    unittest.main()

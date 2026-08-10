from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_ROOT))


def load(name: str, filename: str):
    specification = importlib.util.spec_from_file_location(name, SCRIPT_ROOT / filename)
    assert specification and specification.loader
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


AUDIT = load("plan45_external_audit", "26_seal_external_guide_audits.py")
FREEZE = load("plan45_target_freeze", "27_adjudicate_and_freeze_stage_a_targets.py")


class ExternalGuideAuditTests(unittest.TestCase):
    def test_dna_and_dual_review_gates(self) -> None:
        self.assertTrue(AUDIT.valid_dna("ACGT" * 5, 18, 30))
        self.assertFalse(AUDIT.valid_dna("ACGTN" * 4, 18, 30))
        row = {"reviewer_1": "A", "reviewer_2": "B", "review_concordant": "true"}
        self.assertTrue(AUDIT.dual_review_pass(row))
        row["reviewer_2"] = "A"
        self.assertFalse(AUDIT.dual_review_pass(row))

    def test_stage_a_pass_requires_every_design_gate(self) -> None:
        row = {
            "audit_row_type": "guide",
            "design_status": "design_available",
            "source_recommended": "true",
            "sequence_qc_pass": "true",
            "off_target_review_pass": "true",
            "review_concordant": "true",
        }
        self.assertTrue(FREEZE.stage_a_guide_pass(row))
        for field in [
            "source_recommended", "sequence_qc_pass", "off_target_review_pass",
            "review_concordant",
        ]:
            failed = dict(row)
            failed[field] = "false"
            self.assertFalse(FREEZE.stage_a_guide_pass(failed), field)

    def test_stage_b_pass_requires_exact_edit_match(self) -> None:
        row = {
            "audit_row_type": "design",
            "design_status": "design_available",
            "source_recommended": "true",
            "complete_design": "true",
            "intended_edit_matches_worklist": "true",
            "off_target_review_pass": "true",
            "review_concordant": "true",
        }
        self.assertTrue(FREEZE.stage_b_design_pass(row))
        row["intended_edit_matches_worklist"] = "false"
        self.assertFalse(FREEZE.stage_b_design_pass(row))

    def test_distinct_design_selection_rejects_duplicate_sequences(self) -> None:
        rows = [
            {"source_design_rank": "1", "guide_id": "g1", "guide_sequence": "A" * 20},
            {"source_design_rank": "2", "guide_id": "g2", "guide_sequence": "A" * 20},
            {"source_design_rank": "3", "guide_id": "g3", "guide_sequence": "C" * 20},
        ]
        selected = FREEZE.select_distinct(rows, "guide_id", ["guide_sequence"], 2)
        self.assertEqual([row["guide_id"] for row in selected], ["g1", "g3"])


if __name__ == "__main__":
    unittest.main()

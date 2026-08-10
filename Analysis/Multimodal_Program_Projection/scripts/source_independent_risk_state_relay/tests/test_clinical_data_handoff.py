from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_ROOT))


def load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPT_ROOT / filename)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


HANDOFF = load("plan46_clinical_handoff", "51_build_clinical_data_handoff.py")


class ClinicalDataHandoffTests(unittest.TestCase):
    def test_question_and_role_universes_are_fixed(self) -> None:
        self.assertEqual([row[0] for row in HANDOFF.QUESTIONS], [f"CF{i:02d}" for i in range(1, 17)])
        self.assertEqual(len(HANDOFF.ROLES), 6)

    def test_required_human_and_orthogonal_assays_are_present(self) -> None:
        assays = {row[0] for row in HANDOFF.ASSAYS}
        self.assertEqual(
            assays,
            {"bulk_rna", "histology", "tissue_proteomics", "secreted_proteomics", "spatial_transcriptomics_or_proteomics"},
        )

    def test_dictionary_requests_no_direct_identifier(self) -> None:
        fields = {row[0] for row in HANDOFF.DATA_DICTIONARY}
        self.assertFalse(fields & {"name", "date_of_birth", "medical_record_number", "address"})
        self.assertIn("participant_id", fields)

    def test_brief_preserves_two_cohort_and_outcome_blind_design(self) -> None:
        text = HANDOFF.BRIEF.lower()
        self.assertIn("two genuinely independent", text)
        self.assertIn("do not send participant-level molecular outcomes initially", text)


if __name__ == "__main__":
    unittest.main()

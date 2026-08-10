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


HANDOFF = load("plan45_collaborator_handoff", "43_build_stage_a_collaborator_handoff.py")


class CollaboratorHandoffTests(unittest.TestCase):
    def test_feasibility_ids_are_complete_and_unique(self) -> None:
        self.assertEqual([row[0] for row in HANDOFF.QUESTIONS], [f"FQ{i:02d}" for i in range(1, 13)])

    def test_deliverable_ids_are_complete_and_unique(self) -> None:
        self.assertEqual([row[0] for row in HANDOFF.DELIVERABLES], [f"D{i:02d}" for i in range(1, 15)])

    def test_handoff_contains_no_named_target(self) -> None:
        legacy_targets = {"ACSL5", "ERCC2", "ZBTB41", "IL18R1", "SHMT1", "PNPLA6", "CENPQ", "NECAB2", "RANBP17"}
        upper = HANDOFF.HANDOFF.upper()
        self.assertFalse(any(target in upper for target in legacy_targets))

    def test_stage_b_questions_are_explicitly_separate(self) -> None:
        impacts = {row[0]: row[3] for row in HANDOFF.QUESTIONS}
        self.assertEqual(impacts["FQ11"], "stage_b_blocking")
        self.assertEqual(impacts["FQ12"], "stage_b_blocking")

    def test_sealed_bundle_and_response_are_separated(self) -> None:
        self.assertIn("read-only", HANDOFF.HANDOFF.lower())
        response = HANDOFF.RESPONSE_INSTRUCTIONS.lower()
        self.assertIn("outside this candidate root", response)
        self.assertIn("never write the response back", response)


if __name__ == "__main__":
    unittest.main()

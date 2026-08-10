from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


SCRIPT_ROOT = Path(__file__).resolve().parents[1]


def load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPT_ROOT / filename)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RESPONSE = load("plan45_collaborator_response", "47_seal_and_adjudicate_collaborator_response.py")


class CollaboratorResponseTests(unittest.TestCase):
    @staticmethod
    def answers(value: str = "yes") -> dict[str, str]:
        return {f"FQ{i:02d}": value for i in range(1, 13)}

    def test_all_yes_permits_target_disclosure_not_target_freeze(self) -> None:
        result = RESPONSE.classify(self.answers(), True, True)
        self.assertEqual(result["feasibility_verdict"], "complete_full_mechanism_feasible")
        self.assertTrue(result["target_disclosure_permitted"])
        self.assertNotIn("target_frozen", result)

    def test_stage_b_no_blocks_full_mechanism_and_disclosure(self) -> None:
        answers = self.answers()
        answers["FQ11"] = "no"
        result = RESPONSE.classify(answers, True, True)
        self.assertEqual(result["feasibility_verdict"], "complete_stage_a_only_not_full_mechanism")
        self.assertTrue(result["stage_a_planning_permitted"])
        self.assertFalse(result["target_disclosure_permitted"])

    def test_unresolved_stage_a_is_pending_and_blocks_planning(self) -> None:
        answers = self.answers()
        answers["FQ04"] = "unresolved"
        result = RESPONSE.classify(answers, True, True)
        self.assertEqual(result["feasibility_verdict"], "pending_stage_a_feasibility")
        self.assertFalse(result["stage_a_planning_permitted"])

    def test_stage_a_no_is_complete_infeasible(self) -> None:
        answers = self.answers()
        answers["FQ06"] = "no"
        result = RESPONSE.classify(answers, True, True)
        self.assertEqual(result["feasibility_verdict"], "complete_stage_a_infeasible")
        self.assertFalse(result["target_disclosure_permitted"])

    def test_role_or_d01_failure_blocks_disclosure(self) -> None:
        self.assertFalse(RESPONSE.classify(self.answers(), False, True)["target_disclosure_permitted"])
        self.assertFalse(RESPONSE.classify(self.answers(), True, False)["target_disclosure_permitted"])

    def test_response_code_contains_no_named_target(self) -> None:
        text = "\n".join(
            (SCRIPT_ROOT / name).read_text(encoding="utf-8")
            for name in [
                "45_build_collaborator_response_scaffold.py",
                "46_validate_collaborator_response_scaffold.py",
                "47_seal_and_adjudicate_collaborator_response.py",
                "48_validate_collaborator_response_release.py",
            ]
        ).upper()
        legacy = {"ACSL5", "ERCC2", "ZBTB41", "IL18R1", "SHMT1", "PNPLA6", "CENPQ", "NECAB2", "RANBP17"}
        self.assertFalse(any(gene in text for gene in legacy))


if __name__ == "__main__":
    unittest.main()

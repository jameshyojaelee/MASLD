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


AUDIT = load("plan45_goal_completion_audit", "49_build_goal_completion_audit.py")


class GoalCompletionAuditTests(unittest.TestCase):
    def test_full_goal_requirement_universe_is_complete(self) -> None:
        self.assertEqual([row[0] for row in AUDIT.REQUIREMENTS], [f"UG{i:02d}" for i in range(1, 12)])
        self.assertFalse(any(row[2] == "achieved" for row in AUDIT.REQUIREMENTS))

    def test_new_data_package_contains_two_human_cohorts_and_exact_edit(self) -> None:
        components = {row[0]: row[1] for row in AUDIT.NEW_DATA}
        self.assertIn("cohort A", components["NDP01"])
        self.assertIn("cohort B", components["NDP02"])
        self.assertIn("Exact endogenous", components["NDP07"])

    def test_full_architecture_requires_every_load_bearing_component(self) -> None:
        architecture = next(row for row in AUDIT.CLAIMS if row[0] == "CL04")
        self.assertEqual(set(architecture[2].split(";")), {f"UG{i:02d}" for i in range(2, 10)})

    def test_no_technical_unit_is_declared_as_biological_unit(self) -> None:
        prohibited = {"cell", "spot", "section", "guide", "well", "organoid"}
        self.assertFalse(any(row[2].lower() in prohibited for row in AUDIT.NEW_DATA))


if __name__ == "__main__":
    unittest.main()

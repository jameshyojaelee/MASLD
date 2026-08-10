#!/usr/bin/env python3
"""Fail-closed fixtures for Plan 46A clinical outreach."""

from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "55_build_clinical_outreach_packets.py"
SPEC = importlib.util.spec_from_file_location("clinical_outreach", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class ClinicalOutreachTests(unittest.TestCase):
    def test_route_and_independence_universe(self) -> None:
        self.assertEqual(set(MODULE.ROUTES), {"nash_crn_as116", "maestro_nash", "essence"})
        self.assertEqual(len({route["independence_group"] for route in MODULE.ROUTES.values()}), 3)
        self.assertEqual(sum(route["cohort_role"] == "candidate_cohort_A" for route in MODULE.ROUTES.values()), 1)
        self.assertEqual(sum(route["cohort_role"] == "candidate_cohort_B" for route in MODULE.ROUTES.values()), 2)

    def test_public_inventory_preserves_unverified_gates(self) -> None:
        for route in MODULE.ROUTES.values():
            self.assertEqual({row[2] for row in route["inventory"]}, {"public_aggregate", "unverified_route_gate"})

    def test_route_questions_are_metadata_only(self) -> None:
        for route in MODULE.ROUTES.values():
            self.assertGreaterEqual(len(route["questions"]), 6)
            text = " ".join(question[2].lower() for question in route["questions"])
            self.assertNotIn("show the expression result", text)
            self.assertNotIn("positive molecular result", text)

    def test_response_boundary_prohibits_outcomes_and_permission_assumption(self) -> None:
        text = " ".join(MODULE.RESPONSE.lower().split())
        for phrase in ("do not send participant-level", "does not grant permission", "does not", "promote the paper"):
            self.assertIn(phrase, text)

    def test_no_route_is_marked_as_accessible(self) -> None:
        for route in MODULE.ROUTES.values():
            self.assertNotIn("access_granted", route)
            self.assertNotIn("permission_granted", route)


if __name__ == "__main__":
    unittest.main()

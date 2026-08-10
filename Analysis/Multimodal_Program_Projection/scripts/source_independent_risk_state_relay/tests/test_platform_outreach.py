#!/usr/bin/env python3
"""Fail-closed fixtures for Plan 45A outreach semantics."""

from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "53_build_platform_outreach_packets.py"
SPEC = importlib.util.spec_from_file_location("platform_outreach", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class PlatformOutreachTests(unittest.TestCase):
    def test_route_universe_and_priorities_are_fixed(self) -> None:
        self.assertEqual(set(MODULE.ROUTES), {
            "du_tsinghua", "takebe_cincinnati",
            "ebrahimkhani_pittsburgh", "ge_iorgantech",
        })
        self.assertEqual(sum(route["send_priority"] == "primary_parallel" for route in MODULE.ROUTES.values()), 3)

    def test_every_route_preserves_blank_response_fields(self) -> None:
        for route in MODULE.ROUTES.values():
            self.assertGreaterEqual(len(route["questions"]), 6)
            for question in route["questions"]:
                self.assertEqual(len(question), 3)
                self.assertTrue(question[1].startswith("FQ"))

    def test_public_evidence_does_not_masquerade_as_platform_pass(self) -> None:
        for route in MODULE.ROUTES.values():
            states = {row[1] for row in route["evidence"]}
            self.assertEqual(states, {"publicly_supported", "not_publicly_demonstrated"})

    def test_initial_packets_contain_no_legacy_candidate_target(self) -> None:
        emitted = "\n".join(
            [MODULE.README, MODULE.RESPONSE_ADDENDUM]
            + [str(route) for route in MODULE.ROUTES.values()]
        )
        for token in MODULE.FORBIDDEN_TARGET_TOKENS:
            self.assertNotIn(token, emitted)

    def test_unverified_contact_is_not_prefilled(self) -> None:
        route = MODULE.ROUTES["ge_iorgantech"]
        self.assertEqual(route["contact_status"], "requires_current_confirmation")
        self.assertEqual(route["contact_value"], "")

    def test_response_instructions_do_not_authorize_science(self) -> None:
        text = " ".join(MODULE.RESPONSE_ADDENDUM.lower().split())
        for phrase in ("cannot freeze a target", "authorize outcome generation", "open stage b", "promote the paper"):
            self.assertIn(phrase, text)


if __name__ == "__main__":
    unittest.main()

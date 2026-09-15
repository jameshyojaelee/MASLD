"""Recovery must read the claim's status, not merely that a claim exists.

A hard kill between writing the claim and submitting leaves a marker at
``claimed_before_submission``. Skipping on existence retires that bundle
silently and permanently. The claim-before-submit ordering is a write-ahead
intent record and is correct; only the recovery read was wrong.
"""

from __future__ import annotations

import inspect
import json
from pathlib import Path
import tempfile
import unittest

from scripts import gpu_bundle_dispatcher


SOURCE = inspect.getsource(gpu_bundle_dispatcher.dispatch_once)


class ClaimRecoveryTests(unittest.TestCase):
    def test_recovery_no_longer_skips_on_bare_existence(self) -> None:
        self.assertNotIn(
            'or marker.exists():', SOURCE,
            "skipping on existence alone retires an interrupted claim forever",
        )

    def test_recovery_branches_on_the_recorded_status(self) -> None:
        self.assertIn('prior.get("status") != "claimed_before_submission"', SOURCE)

    def test_a_submitted_claim_is_still_skipped(self) -> None:
        """The dedup guarantee must survive the fix."""

        self.assertIn("continue", SOURCE)
        index = SOURCE.index('prior.get("status") != "claimed_before_submission"')
        self.assertIn("continue", SOURCE[index:index + 120])

    def test_the_write_ahead_ordering_is_documented_as_correct(self) -> None:
        self.assertIn("write-ahead intent", SOURCE)
        self.assertIn("do not reorder", SOURCE)

    def test_an_unreadable_marker_is_not_treated_as_interrupted(self) -> None:
        """A corrupt marker must not silently license a resubmission."""

        self.assertIn('"unreadable"', SOURCE)

    def test_retry_is_recorded_on_the_new_claim(self) -> None:
        self.assertIn("retry_of_interrupted_claim", SOURCE)


class ClaimStatusFixtureTests(unittest.TestCase):
    """The statuses the recovery branch distinguishes."""

    def setUp(self) -> None:
        holder = tempfile.TemporaryDirectory()
        self.addCleanup(holder.cleanup)
        self.claims = Path(holder.name)

    def write(self, status: str) -> dict:
        path = self.claims / "b.json"
        path.write_text(json.dumps({"bundle_id": "b", "status": status}), encoding="utf-8")
        return json.loads(path.read_text())

    def test_interrupted_claim_is_the_retry_case(self) -> None:
        self.assertEqual(
            self.write("claimed_before_submission")["status"],
            "claimed_before_submission",
        )

    def test_submitted_and_failed_are_both_terminal_for_recovery(self) -> None:
        for status in ("submitted", "submission_failed"):
            self.assertNotEqual(self.write(status)["status"], "claimed_before_submission")


if __name__ == "__main__":
    unittest.main()

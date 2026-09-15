"""The census must not advertise MPRALegNet as scoreable.

OVERALL_PLAN records MPRALegNet as terminally blocked from scoring, but the
terminal receipt was bound nowhere under config/, so the model census read as
conditionally open.  These tests pin the correction.
"""
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DISPOSITION = ROOT / "config/artifacts/models/mpralegnet/terminal_disposition_20260825.json"
RECEIPT = ROOT / "executions/model-check-230-21092108"


class TestMpralegnetTerminalDisposition(unittest.TestCase):
    def setUp(self):
        self.disp = json.loads(DISPOSITION.read_text())

    def test_disposition_exists_and_is_record_only(self):
        self.assertTrue(self.disp["record_only"])
        self.assertFalse(self.disp["consumed_by_runtime"])
        self.assertTrue(self.disp["base_registry_unchanged"])

    def test_scoring_ranking_and_champion_are_all_blocked(self):
        c = self.disp["binding_consequences"]
        for key in ("scoring_authorized", "may_produce_benchmark_metrics",
                    "may_enter_a_shortlist_or_ranking", "open_champion_eligible",
                    "standalone_champion_eligible"):
            self.assertFalse(c[key], f"{key} must be false")

    def test_binds_the_exact_terminal_receipt(self):
        r = self.disp["terminal_receipt"]
        self.assertEqual(r["path"], "executions/model-check-230-21092108")
        self.assertEqual(
            r["artifacts_sha256"],
            "cdf356b81aa19127f2eff45a8d213a4a9d360b9fed5f81edb4c7972d0cb9874d",
        )

    def test_disposition_agrees_with_the_receipt_it_binds(self):
        meta = json.loads((RECEIPT / "ARTIFACTS.json").read_text())["metadata"]
        r = self.disp["terminal_receipt"]
        self.assertEqual(meta["status"], "terminal_blocked")
        self.assertEqual(r["status"], meta["status"])
        self.assertEqual(r["terminal_reason"], meta["terminal_reason"])
        self.assertFalse(meta["provenance_equivalence_pass"])
        self.assertFalse(meta["standalone_champion_eligible"])

    def test_records_pending_versus_terminally_failed_distinction(self):
        d = self.disp["defect_corrected"]
        self.assertIn("PENDING", d["precise_discrepancy"])
        self.assertIn("FAILED", d["precise_discrepancy"])

    def test_ridge_seed_caveat_is_retained(self):
        self.assertIn(
            "not five stability replicates",
            self.disp["binding_consequences"]["ridge_seed_caveat"],
        )

    def test_no_outcome_or_metric_was_touched(self):
        for key in ("outcomes_read", "metrics_calculated", "sealed_data_read"):
            self.assertFalse(self.disp[key])


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

from pathlib import Path
import unittest

from scripts.audit_corgi_regular_smoke_plan import audit_smoke_plan


ROOT = Path(__file__).resolve().parents[2]
AUTHORITY = (
    ROOT
    / "config/artifacts/models/corgi/bounded_prediction_smoke_disposition_20260824.json"
)


class CorgiRegularSmokePlanAuditTests(unittest.TestCase):
    def test_bounded_prediction_is_complete_without_native_admission(self) -> None:
        receipt = audit_smoke_plan(ROOT, AUTHORITY)
        self.assertEqual(
            receipt["status"], "pass_bounded_outcome_free_prediction_smoke"
        )
        self.assertEqual(
            receipt["admission_class"], "restricted_development_comparator"
        )
        self.assertIs(receipt["native_lane_admitted"], False)
        self.assertIs(receipt["open_champion_eligible"], False)
        self.assertIs(receipt["resubmission_needed"], False)
        self.assertIs(receipt["prediction_arrays_opened"], False)
        self.assertIs(receipt["evaluator_outputs_opened"], False)
        self.assertIs(receipt["outcomes_opened"], False)
        self.assertIs(receipt["sealed_features_or_labels_opened"], False)
        self.assertEqual(
            receipt["corgi_plus_status"],
            "TERMINAL_BLOCKED_CURRENT_UPSTREAM_RELEASE",
        )
        self.assertIs(receipt["regular_corgi_relabelled_as_corgi_plus"], False)


if __name__ == "__main__":
    unittest.main()

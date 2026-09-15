from __future__ import annotations

from pathlib import Path
import unittest

from scripts.audit_chrombpnet_full_head_readiness import audit_readiness


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = (
    ROOT
    / "config/artifacts/models/chrombpnet/full_head_production_readiness_20260825.json"
)


class ChromBPNetFullHeadReadinessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.receipt = audit_readiness(ROOT, CONTRACT)

    def test_full_and_nobias_views_have_distinct_scientific_roles(self) -> None:
        receipt = self.receipt
        self.assertEqual(
            receipt["status"],
            "pass_full_head_reconciled_waiting_for_rectangle_and_license",
        )
        self.assertTrue(receipt["full_and_nobias_share_one_fit"])
        self.assertEqual(receipt["primary_biological_view"], "chrombpnet_nobias")
        self.assertEqual(receipt["assay_qc_view"], "chrombpnet_full")
        self.assertFalse(receipt["full_view_is_independent_architecture"])

    def test_five_seed_rectangle_remains_incomplete_and_evaluator_closed(self) -> None:
        receipt = self.receipt
        self.assertEqual(receipt["chrombpnet_expected_fits"], 125)
        self.assertEqual(receipt["chrombpnet_complete_fits_at_freeze"], 19)
        self.assertEqual(receipt["chrombpnet_missing_fits_at_freeze"], 106)
        self.assertEqual(receipt["chrombpnet_expected_prediction_views"], 250)
        self.assertEqual(receipt["chrombpnet_complete_prediction_views_at_freeze"], 38)
        self.assertEqual(receipt["chrombpnet_missing_prediction_views_at_freeze"], 212)
        self.assertFalse(receipt["full_compatible_rectangle"])
        self.assertFalse(receipt["partial_evaluator_access_authorized"])
        self.assertTrue(receipt["legacy_seed_authority_superseded_for_rectangle"])

    def test_license_contradiction_blocks_open_champion_and_redistribution(self) -> None:
        receipt = self.receipt
        self.assertTrue(receipt["license_authority_contradiction_present"])
        self.assertFalse(receipt["weight_redistribution_authorized"])
        self.assertFalse(receipt["open_champion_eligible_now"])
        self.assertFalse(receipt["target_label_unexposed_equals_champion_eligible"])

    def test_existing_dispatch_unit_is_used_without_direct_gpu_submission(self) -> None:
        receipt = self.receipt
        self.assertTrue(receipt["valid_next_gpu_unit_already_registered"])
        self.assertEqual(receipt["next_dispatcher_bundle_id"], "model-training-801")
        self.assertFalse(receipt["new_gpu_unit_registration_needed"])
        self.assertFalse(receipt["direct_gpu_submission_authorized"])
        self.assertFalse(receipt["gpu_job_submitted"])

    def test_outcome_firewall_remains_closed(self) -> None:
        receipt = self.receipt
        self.assertFalse(receipt["prediction_values_read"])
        self.assertFalse(receipt["metric_values_read"])
        self.assertFalse(receipt["development_outcomes_read"])
        self.assertFalse(receipt["sealed_features_read"])
        self.assertFalse(receipt["sealed_labels_or_outcomes_read"])
        self.assertFalse(receipt["thresholds_changed"])


if __name__ == "__main__":
    unittest.main()

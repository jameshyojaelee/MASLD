from __future__ import annotations

from pathlib import Path
import unittest

from scripts.audit_chrombpnet_full_head_rectangle_readiness import audit_readiness


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = (
    ROOT
    / "config/artifacts/models/chrombpnet/full_head_rectangle_execution_readiness_20260825.json"
)


class ChromBPNetFullHeadRectangleReadinessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.receipt = audit_readiness(ROOT, CONTRACT)

    def test_one_fit_produces_two_distinct_chrombpnet_views(self) -> None:
        receipt = self.receipt
        self.assertTrue(receipt["one_fit_two_views"])
        self.assertEqual(receipt["primary_biological_view"], "chrombpnet_nobias")
        self.assertEqual(receipt["assay_qc_view"], "chrombpnet_full")
        self.assertFalse(receipt["full_view_is_independent_architecture"])

    def test_exact_five_seed_rectangle_remains_incomplete(self) -> None:
        receipt = self.receipt
        self.assertEqual(receipt["expected_fits"], 125)
        self.assertEqual(receipt["complete_fits_at_freeze"], 19)
        self.assertEqual(receipt["missing_fits_at_freeze"], 106)
        self.assertEqual(receipt["expected_prediction_views"], 250)
        self.assertEqual(receipt["complete_prediction_views_at_freeze"], 38)
        self.assertEqual(receipt["missing_prediction_views_at_freeze"], 212)
        self.assertFalse(receipt["full_compatible_rectangle"])
        self.assertFalse(receipt["partial_evaluator_access_authorized"])

    def test_all_missing_fits_are_already_in_lane_specific_dispatch_bundles(self) -> None:
        receipt = self.receipt
        self.assertEqual(receipt["missing_fits_already_centrally_registered"], 106)
        self.assertEqual(receipt["central_task_bundles"], 7)
        self.assertEqual(
            receipt["central_dispatcher_bundle_ids"],
            [f"model-training-{value}" for value in range(808, 815)],
        )
        self.assertFalse(receipt["new_gpu_queue_registration_needed"])
        self.assertFalse(receipt["direct_gpu_submission_authorized"])
        self.assertFalse(receipt["gpu_job_submitted"])

    def test_gpu_and_license_gates_remain_fail_closed(self) -> None:
        receipt = self.receipt
        self.assertEqual(receipt["qos"], "nslab")
        self.assertEqual(receipt["gpu_squeue_ceiling"], 5)
        self.assertEqual(receipt["maximum_running_gpu_jobs"], 4)
        self.assertEqual(receipt["maximum_pending_gpu_jobs"], 1)
        self.assertTrue(receipt["license_authority_contradiction_present"])
        self.assertFalse(receipt["weight_redistribution_authorized"])
        self.assertFalse(receipt["open_champion_eligible_now"])

    def test_outcome_firewall_and_partial_ranking_remain_closed(self) -> None:
        receipt = self.receipt
        self.assertFalse(receipt["partial_results_ranked"])
        self.assertFalse(receipt["prediction_values_read"])
        self.assertFalse(receipt["metric_values_read"])
        self.assertFalse(receipt["development_outcomes_read"])
        self.assertFalse(receipt["sealed_features_read"])
        self.assertFalse(receipt["sealed_labels_or_outcomes_read"])
        self.assertFalse(receipt["thresholds_changed"])


if __name__ == "__main__":
    unittest.main()

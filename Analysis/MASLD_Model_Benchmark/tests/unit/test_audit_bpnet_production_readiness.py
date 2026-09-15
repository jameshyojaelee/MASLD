from __future__ import annotations

from pathlib import Path
import unittest

from scripts.audit_bpnet_production_readiness import audit_readiness


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "config/artifacts/models/bpnet/production_readiness_20260825.json"


class BPNetProductionReadinessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.receipt = audit_readiness(ROOT, CONTRACT)

    def test_bpnet_is_not_a_chrombpnet_prediction_view(self) -> None:
        receipt = self.receipt
        self.assertTrue(receipt["bpnet_is_separate_architecture_from_chrombpnet"])
        self.assertFalse(receipt["bpnet_is_chrombpnet_nobias_view"])
        self.assertEqual(receipt["prediction_views_per_bpnet_fit"], 1)

    def test_exact_prospective_rectangle_is_reconciled_without_partial_ranking(self) -> None:
        receipt = self.receipt
        self.assertEqual(receipt["expected_fits"], 125)
        self.assertEqual(receipt["complete_fits_at_freeze"], 20)
        self.assertEqual(receipt["missing_fits_at_freeze"], 105)
        self.assertFalse(receipt["full_compatible_rectangle"])
        self.assertFalse(receipt["partial_evaluator_access_authorized"])
        self.assertFalse(receipt["partial_results_ranked"])

    def test_all_missing_fits_are_already_in_bundled_central_dispatch(self) -> None:
        receipt = self.receipt
        self.assertEqual(receipt["missing_fits_already_centrally_registered"], 105)
        self.assertEqual(receipt["central_task_bundles"], 7)
        self.assertEqual(
            receipt["central_dispatcher_bundle_ids"],
            [f"model-training-{value}" for value in range(801, 808)],
        )
        self.assertFalse(receipt["new_gpu_queue_registration_needed"])
        self.assertFalse(receipt["direct_gpu_submission_authorized"])
        self.assertFalse(receipt["gpu_job_submitted"])

    def test_gpu_governance_and_release_gate_remain_closed(self) -> None:
        receipt = self.receipt
        self.assertEqual(receipt["qos"], "nslab")
        self.assertEqual(receipt["gpu_squeue_ceiling"], 5)
        self.assertEqual(receipt["maximum_running_gpu_jobs"], 4)
        self.assertEqual(receipt["maximum_pending_gpu_jobs"], 1)
        self.assertTrue(receipt["code_and_project_weight_redistribution_authorized"])
        self.assertFalse(receipt["open_champion_eligible_now"])

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

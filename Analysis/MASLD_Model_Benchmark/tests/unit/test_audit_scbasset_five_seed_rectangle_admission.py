from __future__ import annotations

from pathlib import Path
import unittest

from scripts.audit_scbasset_five_seed_rectangle_admission import audit_admission


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "config/campaigns/scbasset_five_seed_rectangle_20260825.json"


class ScBassetFiveSeedRectangleAdmissionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.receipt = audit_admission(ROOT, CONTRACT)

    def test_exact_five_seed_rectangle_is_reconciled(self) -> None:
        receipt = self.receipt
        self.assertEqual(receipt["expected_fits"], 25)
        self.assertEqual(receipt["complete_fits_at_freeze"], 4)
        self.assertEqual(receipt["missing_fits_admitted"], 21)
        self.assertEqual(receipt["expected_valid_prediction_views"], 25)
        self.assertEqual(receipt["complete_valid_prediction_views_at_freeze"], 4)
        self.assertEqual(receipt["missing_valid_prediction_views_admitted"], 21)

    def test_prediction_view_preserves_sequence_only_claim_boundary(self) -> None:
        receipt = self.receipt
        self.assertEqual(
            receipt["prediction_view"],
            "valid_training_cell_state_mean_sequence_profile",
        )
        self.assertEqual(receipt["donor_context"], "none_sequence_only")
        self.assertFalse(receipt["observed_query_atac_consumed"])
        self.assertFalse(receipt["rna_context_consumed"])
        self.assertFalse(receipt["held_cell_embedding_available"])
        self.assertFalse(receipt["observed_multiome_inductive_claim_allowed"])
        self.assertFalse(receipt["rna_conditioned_claim_allowed"])

    def test_missing_coverage_is_exactly_partitioned_across_two_bundles(self) -> None:
        receipt = self.receipt
        self.assertEqual(receipt["existing_bundle_id"], "model-training-604")
        self.assertEqual(receipt["existing_bundle_missing_fits"], 11)
        self.assertEqual(receipt["extension_bundle_id"], "model-training-605")
        self.assertEqual(receipt["extension_bundle_missing_fits"], 10)
        self.assertEqual(receipt["central_bundle_count"], 2)
        self.assertEqual(receipt["logical_tasks"], 42)

    def test_new_gpu_unit_is_dispatcher_only_and_governed(self) -> None:
        receipt = self.receipt
        self.assertTrue(receipt["central_queue_registration_present"])
        self.assertTrue(receipt["central_queue_unclaimed_during_validation"])
        self.assertEqual(receipt["requested_new_gpu_allocations"], 1)
        self.assertEqual(receipt["requested_new_gpu_hours"], 12)
        self.assertEqual(receipt["gpu_squeue_ceiling"], 5)
        self.assertEqual(receipt["maximum_running_gpu_jobs"], 4)
        self.assertEqual(receipt["maximum_pending_gpu_jobs"], 1)
        self.assertFalse(receipt["direct_gpu_submission_authorized"])
        self.assertFalse(receipt["gpu_job_submitted_by_admission"])

    def test_evaluation_and_outcome_firewall_remain_closed(self) -> None:
        receipt = self.receipt
        self.assertFalse(receipt["partial_ranking_authorized"])
        self.assertFalse(receipt["development_evaluator_open"])
        self.assertFalse(receipt["champion_claim_allowed"])
        self.assertFalse(receipt["admission_audit_scbasset_prediction_values_read"])
        self.assertFalse(receipt["metric_values_read"])
        self.assertFalse(receipt["raw_outcomes_read"])
        self.assertFalse(receipt["sealed_features_read"])
        self.assertFalse(receipt["sealed_labels_or_outcomes_read"])

    def test_development_search_incident_is_disclosed_and_unused(self) -> None:
        receipt = self.receipt
        self.assertTrue(receipt["development_search_incident_recorded"])
        self.assertEqual(
            receipt["development_search_incident_id"],
            "scbasset_broad_rg_prediction_line_exposure_20260825",
        )
        self.assertTrue(receipt["development_search_incident_structure_validated"])
        self.assertTrue(receipt["unrelated_prediction_lines_accidentally_emitted"])
        self.assertFalse(receipt["incident_values_used"])


if __name__ == "__main__":
    unittest.main()

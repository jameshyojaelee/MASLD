from __future__ import annotations

from pathlib import Path
import unittest

from scripts.audit_scbasset_three_seed_rectangle_admission import audit_admission


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "config/campaigns/scbasset_three_seed_rectangle_20260825.json"


class ScBassetThreeSeedRectangleAdmissionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.receipt = audit_admission(ROOT, CONTRACT)

    def test_exact_five_fold_three_seed_rectangle_is_reconciled(self) -> None:
        receipt = self.receipt
        self.assertEqual(
            receipt["status"],
            "pass_exact_11_fit_rectangle_admitted_for_central_dispatch",
        )
        self.assertEqual(receipt["expected_fits"], 15)
        self.assertEqual(receipt["complete_fits_at_freeze"], 4)
        self.assertEqual(receipt["missing_fits_admitted"], 11)
        self.assertEqual(receipt["fixed_seeds"], [20260824, 20260825, 20260826])

    def test_legacy_seed_namespace_is_excluded(self) -> None:
        receipt = self.receipt
        self.assertEqual(receipt["legacy_seeds_excluded"], [11, 29, 47, 71, 101])
        self.assertFalse(set(receipt["fixed_seeds"]) & set(receipt["legacy_seeds_excluded"]))

    def test_one_bundled_central_allocation_owns_all_missing_fits(self) -> None:
        receipt = self.receipt
        self.assertEqual(receipt["central_bundle_id"], "model-training-604")
        self.assertEqual(receipt["central_bundle_count"], 1)
        self.assertEqual(receipt["logical_tasks"], 22)
        self.assertEqual(receipt["requested_gpu_allocations"], 1)
        self.assertEqual(receipt["requested_gpu_hours"], 8)
        self.assertTrue(receipt["central_queue_registration_present"])
        self.assertTrue(receipt["central_queue_unclaimed_during_validation"])
        self.assertFalse(receipt["admission_gate_present_during_validation"])
        self.assertFalse(receipt["direct_gpu_submission_authorized"])
        self.assertFalse(receipt["manual_gpu_sbatch_allowed"])
        self.assertFalse(receipt["gpu_job_submitted_by_admission"])

    def test_supersession_does_not_expand_model_or_promotion_scope(self) -> None:
        receipt = self.receipt
        self.assertFalse(receipt["observed_multiome_eligibility_changed"])
        self.assertFalse(receipt["rna_conditioned_model_status_changed"])
        self.assertFalse(receipt["development_evaluator_open"])
        self.assertFalse(receipt["champion_claim_allowed"])
        self.assertFalse(receipt["partial_ranking_authorized"])

    def test_prediction_metric_outcome_and_sealed_firewalls_remain_closed(self) -> None:
        receipt = self.receipt
        self.assertFalse(receipt["prediction_values_read"])
        self.assertFalse(receipt["metric_values_read"])
        self.assertFalse(receipt["raw_outcomes_read"])
        self.assertFalse(receipt["sealed_features_read"])
        self.assertFalse(receipt["sealed_labels_or_outcomes_read"])


if __name__ == "__main__":
    unittest.main()

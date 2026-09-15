from __future__ import annotations

from pathlib import Path
import unittest

from scripts.audit_sequence_control_production_readiness import audit_readiness


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = (
    ROOT
    / "config/artifacts/models/regulatory_native_controls/production_readiness_20260825.json"
)


class SequenceControlProductionReadinessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.receipt = audit_readiness(ROOT, CONTRACT)

    def test_cnn_and_transformer_remain_separate_matched_controls(self) -> None:
        receipt = self.receipt
        self.assertEqual(
            receipt["status"],
            "pass_controls_reconciled_waiting_for_full_rectangles_license_and_external_evaluation",
        )
        self.assertEqual(
            receipt["model_ids"],
            ["sequence_cnn_control", "sequence_transformer_control"],
        )
        self.assertTrue(receipt["controls_are_separate_architectures"])
        self.assertFalse(receipt["controls_share_checkpoint"])
        self.assertTrue(receipt["parameter_budget_matched"])
        self.assertEqual(receipt["cnn_parameter_count"], 3_696_386)
        self.assertEqual(receipt["transformer_parameter_count"], 3_683_330)

    def test_donor_and_whole_contig_isolation_remain_locked(self) -> None:
        receipt = self.receipt
        self.assertEqual(receipt["biological_donors"], 39)
        self.assertEqual(receipt["donor_fold_counts"], [11, 9, 7, 4, 8])
        self.assertEqual(receipt["genomic_outer_unit"], "whole_chromosome_group")
        self.assertFalse(receipt["whole_contigs_shared_across_genomic_folds"])
        self.assertFalse(receipt["linear_boundary_buffer_required"])

    def test_exact_five_by_five_by_five_rectangles_remain_incomplete(self) -> None:
        receipt = self.receipt
        self.assertEqual(receipt["cnn_expected_fits"], 125)
        self.assertEqual(receipt["cnn_complete_fits_at_freeze"], 12)
        self.assertEqual(receipt["cnn_missing_fits_at_freeze"], 113)
        self.assertEqual(receipt["transformer_expected_fits"], 125)
        self.assertEqual(receipt["transformer_complete_fits_at_freeze"], 11)
        self.assertEqual(receipt["transformer_missing_fits_at_freeze"], 114)
        self.assertFalse(receipt["full_compatible_rectangle"])

    def test_partial_diagnostics_cannot_rank_or_open_evaluator(self) -> None:
        receipt = self.receipt
        self.assertFalse(receipt["partial_evaluator_access_authorized"])
        self.assertFalse(receipt["partial_results_ranked"])
        self.assertFalse(receipt["prediction_values_read"])
        self.assertFalse(receipt["metric_values_read"])
        self.assertFalse(receipt["development_outcomes_read"])
        self.assertFalse(receipt["sealed_labels_or_outcomes_read"])

    def test_project_license_and_external_evaluation_block_open_champion(self) -> None:
        receipt = self.receipt
        self.assertEqual(receipt["upstream_code_license"], "MIT")
        self.assertFalse(receipt["project_repository_license_resolved"])
        self.assertFalse(receipt["source_redistribution_authorized"])
        self.assertFalse(receipt["derivative_weight_redistribution_authorized"])
        self.assertFalse(receipt["open_champion_eligible_now"])

    def test_existing_central_queue_owns_remaining_work_and_tasks_stay_separate(self) -> None:
        receipt = self.receipt
        self.assertEqual(receipt["cnn_missing_fits_already_centrally_registered"], 113)
        self.assertEqual(
            receipt["transformer_missing_fits_already_centrally_registered"], 114
        )
        self.assertEqual(receipt["cnn_next_dispatcher_bundle_id"], "model-training-815")
        self.assertEqual(
            receipt["transformer_next_dispatcher_bundle_id"], "model-training-823"
        )
        self.assertFalse(receipt["new_gpu_queue_registration_needed"])
        self.assertFalse(receipt["direct_gpu_submission_authorized"])
        self.assertFalse(receipt["gse281364_fits_fill_gse296875_rectangle"])
        self.assertFalse(receipt["gpu_job_submitted"])


if __name__ == "__main__":
    unittest.main()

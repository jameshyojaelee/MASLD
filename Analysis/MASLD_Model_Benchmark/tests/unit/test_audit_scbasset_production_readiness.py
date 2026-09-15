from __future__ import annotations

from pathlib import Path
import unittest

from scripts.audit_scbasset_production_readiness import audit_readiness


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "config/artifacts/models/scbasset/production_readiness_20260825.json"


class ScBassetProductionReadinessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.receipt = audit_readiness(ROOT, CONTRACT)

    def test_native_semantics_are_training_atac_supervised_but_sequence_only_at_query(self) -> None:
        receipt = self.receipt
        self.assertEqual(
            receipt["status"], "pass_reconciled_development_only_no_gpu_unit"
        )
        self.assertEqual(
            receipt["native_observed_atac_role"], "outer_training_supervision_only"
        )
        self.assertFalse(receipt["observed_query_atac_consumed"])
        self.assertFalse(receipt["inductive_observed_multiome_model"])
        self.assertFalse(receipt["rna_conditioned_model"])
        self.assertEqual(receipt["donor_context"], "none_sequence_only")
        self.assertFalse(receipt["held_cell_embedding_available"])

    def test_all_five_valid_is_not_all_seeds_or_test_role(self) -> None:
        receipt = self.receipt
        self.assertTrue(receipt["five_outer_fold_valid_campaign_present"])
        self.assertFalse(receipt["test_role_predictions_present"])
        self.assertTrue(receipt["mixed_seed_meta_aggregate"])
        self.assertFalse(receipt["uniform_five_fold_prediction_lock_complete"])
        self.assertEqual(receipt["fixed_seed_rectangle_expected_fits"], 15)
        self.assertEqual(receipt["fixed_seed_rectangle_complete_fits"], 4)
        self.assertEqual(receipt["fixed_seed_rectangle_missing_fits"], 11)
        self.assertFalse(receipt["fixed_seed_rectangle_complete"])
        self.assertFalse(receipt["partial_ranking_authorized"])

    def test_fixture_artifact_does_not_override_global_task_registries(self) -> None:
        receipt = self.receipt
        self.assertTrue(receipt["profile_fixture_artifact_exists"])
        self.assertFalse(receipt["global_profile_fixture_admission_complete"])
        self.assertFalse(receipt["observed_multiome_task_eligible"])
        self.assertFalse(receipt["tutorial_weights_admitted"])
        self.assertFalse(receipt["open_champion_eligible_now"])

    def test_later_development_gate_blocks_new_gpu_registration(self) -> None:
        receipt = self.receipt
        self.assertFalse(receipt["development_gate_passed"])
        self.assertTrue(receipt["historical_gpu_item_already_claimed"])
        self.assertFalse(receipt["blocked_plan_candidate_registered"])
        self.assertFalse(receipt["valid_next_centrally_dispatched_gpu_unit_exists"])
        self.assertFalse(receipt["new_gpu_queue_registration_authorized"])
        self.assertFalse(receipt["direct_gpu_submission_authorized"])
        self.assertFalse(receipt["gpu_queue_registered"])
        self.assertFalse(receipt["gpu_job_submitted"])

    def test_prediction_metric_and_sealed_firewalls_remain_closed(self) -> None:
        receipt = self.receipt
        self.assertFalse(receipt["prediction_values_read"])
        self.assertFalse(receipt["metric_values_read"])
        self.assertFalse(receipt["raw_outcomes_read"])
        self.assertFalse(receipt["sealed_features_read"])
        self.assertFalse(receipt["sealed_labels_or_outcomes_read"])
        self.assertFalse(receipt["partial_results_ranked"])
        self.assertFalse(receipt["thresholds_changed"])


if __name__ == "__main__":
    unittest.main()

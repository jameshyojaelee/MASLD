from __future__ import annotations

import unittest

from scripts.evaluate_corgi_regular_crossfit_profile_head import (
    CorgiHeadEvaluationError,
    relative_deviance_skill,
    validate_fit_receipt,
)


class CorgiCrossfitProfileHeadEvaluationTests(unittest.TestCase):
    def test_identical_baseline_and_head_have_zero_skill(self) -> None:
        self.assertEqual(relative_deviance_skill(0.42, 0.42), 0.0)

    def test_positive_gain_has_positive_skill(self) -> None:
        self.assertAlmostEqual(relative_deviance_skill(0.5, 0.4), 0.2)

    def test_receipt_rejects_held_donor_overlap(self) -> None:
        value = self._receipt()
        value["folds"][2]["held_donor_overlap_with_fit"] = 1
        with self.assertRaises(CorgiHeadEvaluationError):
            validate_fit_receipt(value)

    def test_receipt_rejects_prior_evaluation(self) -> None:
        value = self._receipt()
        value["evaluation_performed"] = True
        with self.assertRaises(CorgiHeadEvaluationError):
            validate_fit_receipt(value)

    def test_receipt_accepts_five_folds(self) -> None:
        self.assertEqual(len(validate_fit_receipt(self._receipt())), 5)

    @staticmethod
    def _receipt() -> dict[str, object]:
        folds = []
        for evaluation in range(5):
            folds.append(
                {
                    "evaluation_fold": evaluation,
                    "evaluation_base_outer_fold": (evaluation - 1) % 5,
                    "fit_valid_folds": [fold for fold in range(5) if fold != evaluation],
                    "held_donor_overlap_with_fit": 0,
                    "held_genomic_fold_overlap_with_fit": 0,
                    "same_alpha_applied_to_all_context_arms": True,
                    "held_fold_atac_signal_values_used_during_fit": False,
                    "held_fold_atac_signal_values_used_during_export": False,
                    "evaluation_metrics_computed": False,
                    "test_or_sealed_features_or_outcomes_read": False,
                }
            )
        return {
            "schema_version": "masld-bench-corgi-crossfit-profile-head-fit-v1",
            "status": "pass_training_only_crossfit_head_fit_and_prediction_export",
            "model_id": "corgi_regular",
            "task_id": "rna_conditioned_atac",
            "dataset_id": "gse296875",
            "lineage_id": "hepatocyte",
            "adaptation_rung": "new_assay_head",
            "development_atac_outcomes_read_for_training_only_fit": True,
            "held_fold_atac_signal_values_used_during_each_head_fit": False,
            "histology_or_disease_labels_read": False,
            "test_or_sealed_features_or_outcomes_read": False,
            "benchmark_metrics_computed": False,
            "evaluation_performed": False,
            "promotion_gate_evaluated": False,
            "champion_or_external_claim_allowed": False,
            "global_census_modified": False,
            "folds": folds,
        }


if __name__ == "__main__":
    unittest.main()

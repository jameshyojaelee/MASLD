from __future__ import annotations

import unittest

import numpy as np

from scripts.evaluate_sequence_accessibility import (
    SequenceEvaluationError,
    bind_prediction_model_id,
    fisher_mean,
    multinomial_deviance_per_insertion,
    smooth_distribution,
    spearman_or_zero,
    validate_evaluation_role,
    within_window_deviance_per_insertion,
)


class SequenceAccessibilityEvaluatorTests(unittest.TestCase):
    def test_only_validation_role_is_authorized_before_finalist_lock(self) -> None:
        self.assertEqual(validate_evaluation_role("valid"), "valid")
        with self.assertRaises(SequenceEvaluationError):
            validate_evaluation_role("test")

    def test_model_id_can_bind_from_verified_root_when_subview_omits_it(self) -> None:
        self.assertEqual(bind_prediction_model_id({}, "chrombpnet"), "chrombpnet")

    def test_model_id_rejects_conflicting_subview_declaration(self) -> None:
        with self.assertRaises(SequenceEvaluationError):
            bind_prediction_model_id({"model_id": "other"}, "chrombpnet")

    def test_distribution_smoothing_is_positive_and_normalized(self) -> None:
        values = smooth_distribution([0.0, 2.0, 0.0])
        self.assertTrue(np.all(values > 0))
        self.assertAlmostEqual(float(values.sum()), 1.0)

    def test_perfect_multinomial_prediction_has_zero_deviance(self) -> None:
        truth = np.asarray([1.0, 2.0, 3.0])
        value = multinomial_deviance_per_insertion(truth, truth)
        self.assertLess(value, 1.0e-10)

    def test_within_window_deviance_ignores_between_window_scale(self) -> None:
        truth = np.asarray([[1.0, 3.0], [4.0, 2.0]])
        predicted = truth * np.asarray([[10.0], [0.2]])
        self.assertLess(within_window_deviance_per_insertion(truth, predicted), 1.0e-10)

    def test_spearman_constant_baseline_is_zero(self) -> None:
        self.assertEqual(spearman_or_zero([1, 2, 3], [1, 1, 1]), 0.0)

    def test_fisher_mean_is_symmetric(self) -> None:
        self.assertAlmostEqual(fisher_mean([-0.5, 0.5]), 0.0)


if __name__ == "__main__":
    unittest.main()

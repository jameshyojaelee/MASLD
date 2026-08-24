from __future__ import annotations

import unittest

from masld_bench.evaluators.rna_atac_regularized_development import (
    RNAATACRegularizedEvaluationError,
    relative_deviance_reduction,
    select_strongest_task_native_baseline,
)


class RNAATACRegularizedEvaluatorTests(unittest.TestCase):
    def test_relative_deviance_reduction_has_the_registered_direction(self) -> None:
        self.assertAlmostEqual(relative_deviance_reduction(100.0, 80.0), 0.2)
        self.assertAlmostEqual(relative_deviance_reduction(100.0, 120.0), -0.2)
        with self.assertRaises(RNAATACRegularizedEvaluationError):
            relative_deviance_reduction(0.0, 1.0)

    def test_strongest_baseline_is_selected_only_from_frozen_roster(self) -> None:
        models = {
            "assay_native_pseudobulk": {"total_multinomial_deviance": 90.0},
            "mean_track": {"total_multinomial_deviance": 80.0},
            "shrunken_pseudobulk": {"total_multinomial_deviance": 70.0},
            "trans_only": {"total_multinomial_deviance": 1.0},
        }
        self.assertEqual(
            select_strongest_task_native_baseline(models),
            "shrunken_pseudobulk",
        )
        del models["mean_track"]
        with self.assertRaises(RNAATACRegularizedEvaluationError):
            select_strongest_task_native_baseline(models)


if __name__ == "__main__":
    unittest.main()

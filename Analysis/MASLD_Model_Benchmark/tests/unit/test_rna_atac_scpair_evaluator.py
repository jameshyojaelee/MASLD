from __future__ import annotations

import unittest

from masld_bench.evaluators.rna_atac_scpair_development import (
    RNAATACScPairEvaluationError,
    paired_block_rows,
    relative_deviance_reduction,
    select_strongest_task_native_baseline,
)


class RNAATACScPairEvaluatorTests(unittest.TestCase):
    def test_relative_deviance_reduction_has_the_registered_direction(self) -> None:
        self.assertAlmostEqual(relative_deviance_reduction(100.0, 80.0), 0.2)
        self.assertAlmostEqual(relative_deviance_reduction(100.0, 120.0), -0.2)
        with self.assertRaises(RNAATACScPairEvaluationError):
            relative_deviance_reduction(0.0, 1.0)

    def test_strongest_baseline_is_selected_only_from_frozen_roster(self) -> None:
        models = {
            "assay_native_pseudobulk": {"total_multinomial_deviance": 90.0},
            "mean_track": {"total_multinomial_deviance": 80.0},
            "shrunken_pseudobulk": {"total_multinomial_deviance": 70.0},
            "scpair": {"total_multinomial_deviance": 1.0},
        }
        self.assertEqual(
            select_strongest_task_native_baseline(models),
            "shrunken_pseudobulk",
        )
        del models["mean_track"]
        with self.assertRaises(RNAATACScPairEvaluationError):
            select_strongest_task_native_baseline(models)

    def test_block_pairing_uses_the_exact_common_universe(self) -> None:
        rows = [
            {
                "model_id": "scpair",
                "donor_hash": "donor",
                "block_hash": "block",
                "stratum": "hepatocyte",
                "deviance": "80",
            },
            {
                "model_id": "shrunken_pseudobulk",
                "donor_hash": "donor",
                "block_hash": "block",
                "stratum": "hepatocyte",
                "deviance": "100",
            },
        ]
        paired = paired_block_rows(
            rows,
            candidate_model="scpair",
            baseline_model="shrunken_pseudobulk",
        )
        self.assertEqual(len(paired), 1)
        self.assertAlmostEqual(float(paired[0]["relative_deviance_reduction"]), 0.2)
        rows.pop()
        with self.assertRaises(RNAATACScPairEvaluationError):
            paired_block_rows(
                rows,
                candidate_model="scpair",
                baseline_model="shrunken_pseudobulk",
            )


if __name__ == "__main__":
    unittest.main()

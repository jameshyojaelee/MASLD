"""Unit tests for the outcome-separated scBasset cross-fold evaluator."""

from __future__ import annotations

import inspect
import unittest

from scripts import evaluate_scbasset_valid_crossfold as evaluator
from scripts.evaluate_scbasset_valid_crossfold import (
    ScBassetCrossfoldEvaluationError,
    summarize_unit_rows,
)


LINEAGES = (
    "cholangiocyte",
    "fibroblast",
    "hepatocyte",
    "macrophage",
    "t_cell",
)


class ScBassetCrossfoldEvaluationTests(unittest.TestCase):
    def test_prediction_lock_precedes_outcome_authority_resolution(self) -> None:
        source = inspect.getsource(evaluator.evaluate)
        self.assertLess(
            source.index("load_prediction_lock"),
            source.index('"donor bigWigs"'),
        )
        self.assertEqual(evaluator.EXPECTED_FOLDS, (1, 2, 3, 4))

    def test_summary_selects_strongest_baseline_and_pools_biological_rows(self) -> None:
        rows = []
        for lineage in LINEAGES:
            rows.extend(
                [
                    {"model_id": "scbasset", "stratum": lineage, "deviance_per_insertion": "0.8"},
                    {"model_id": "training_lineage_mean", "stratum": lineage, "deviance_per_insertion": "1.0"},
                    {"model_id": "training_global_mean", "stratum": lineage, "deviance_per_insertion": "1.2"},
                ]
            )
        summary = summarize_unit_rows(rows)
        self.assertEqual(
            summary["strongest_training_only_baseline"], "training_lineage_mean"
        )
        self.assertAlmostEqual(summary["scbasset_relative_deviance_reduction"], 0.2)
        self.assertEqual(summary["development_gate"]["improved_lineages"], 5)
        self.assertTrue(
            summary["development_gate"]["no_lineage_worse_than_threshold"]
        )

    def test_summary_rejects_incomplete_model_roster(self) -> None:
        with self.assertRaisesRegex(
            ScBassetCrossfoldEvaluationError, "model roster differs"
        ):
            summarize_unit_rows(
                [
                    {
                        "model_id": "scbasset",
                        "stratum": "hepatocyte",
                        "deviance_per_insertion": "1.0",
                    }
                ]
            )


if __name__ == "__main__":
    unittest.main()

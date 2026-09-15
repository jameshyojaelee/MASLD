"""Unit tests for the scBasset five-fold frozen-metric aggregate."""

from __future__ import annotations

import unittest

from scripts.aggregate_scbasset_all5_valid import summarize


LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")


class ScBassetAllFoldAggregateTests(unittest.TestCase):
    def test_summary_uses_equal_unit_rows_and_strongest_baseline(self) -> None:
        rows = []
        for lineage in LINEAGES:
            rows.extend(
                [
                    {"model_id": "scbasset", "stratum": lineage, "deviance_per_insertion": "0.8"},
                    {"model_id": "training_lineage_mean", "stratum": lineage, "deviance_per_insertion": "1.0"},
                    {"model_id": "training_global_mean", "stratum": lineage, "deviance_per_insertion": "1.2"},
                ]
            )
        result = summarize(rows)
        self.assertEqual(result["strongest_training_only_baseline"], "training_lineage_mean")
        self.assertAlmostEqual(result["scbasset_relative_deviance_reduction"], 0.2)
        self.assertEqual(result["development_gate"]["improved_lineages"], 5)

    def test_summary_rejects_incomplete_model_roster(self) -> None:
        with self.assertRaisesRegex(ValueError, "model roster differs"):
            summarize([{"model_id": "scbasset", "stratum": "hepatocyte", "deviance_per_insertion": "1"}])


if __name__ == "__main__":
    unittest.main()

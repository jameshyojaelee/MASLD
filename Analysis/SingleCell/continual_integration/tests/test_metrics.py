from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from masld_cl.benchmark import _bundle_centroids, _improved_lineage_count, _validate_harmony_info
from masld_cl.contracts import ContractError
from masld_cl.metrics import (
    bh_adjust,
    donor_centroids,
    donor_distance_spearman,
    macro_f1,
    normalized_shift_control,
    paired_shift_bootstrap,
    standardized_case_control_separation,
    positive_class_f1,
)


class TestMetrics(unittest.TestCase):
    def test_bundle_centroids_does_not_mutate_cell_eligibility(self):
        cells = pd.DataFrame({
            "donor_id": ["r1", "q1", "r2", "q2"],
            "dataset": ["R", "Q", "R", "Q"],
            "audit_cell_type": ["A", "A", "B", "B"],
            "preparation": ["unsorted"] * 4,
            "analysis_eligible": [True] * 4,
            "strict_reference": [True, False, True, False],
            "primary_query": [False, True, False, True],
            "query_control": [False, True, False, True],
        })
        before = cells["analysis_eligible"].copy()
        bundle = ({}, np.arange(8, dtype=float).reshape(4, 2), cells)
        _bundle_centroids(bundle, "A")
        _bundle_centroids(bundle, "B")
        pd.testing.assert_series_equal(cells["analysis_eligible"], before)

    def test_donor_centroids_not_cells_are_units(self):
        latent = np.asarray([[0, 0], [2, 0], [10, 0], [12, 0]], dtype=float)
        donors, centroids = donor_centroids(latent, ["A", "A", "B", "B"])
        self.assertEqual(donors.tolist(), ["A", "B"])
        self.assertTrue(np.array_equal(centroids, [[1, 0], [11, 0]]))

    def test_shift_control_and_disease_separation(self):
        reference = np.asarray([[0, 0], [0, 2], [2, 0]], dtype=float)
        controls = reference + 0.1
        cases = reference + 5
        within = np.mean([
            np.mean([2.0, 2.0, np.sqrt(8.0)]),
            np.mean([2.0, 2.0, np.sqrt(8.0)]),
        ])
        self.assertAlmostEqual(
            normalized_shift_control(reference, controls),
            np.sqrt(0.02) / within,
        )
        self.assertGreater(standardized_case_control_separation(cases, controls), 1)

    def test_paired_bootstrap_detects_improvement(self):
        rng = np.random.default_rng(17)
        reference = rng.normal(size=(12, 3))
        control_baseline = reference[:8] + 2
        control_candidate = reference[:8] + 0.2
        result = paired_shift_bootstrap(
            reference, control_candidate, reference, control_baseline,
            replicates=500, seed=41,
        )
        self.assertGreater(result["estimate"], 0)
        self.assertGreater(result["ci_low"], 0)

    def test_distance_spearman_and_macro_f1(self):
        values = np.asarray([[0, 0], [1, 0], [0, 2], [3, 4]], dtype=float)
        self.assertAlmostEqual(donor_distance_spearman(values, values * 2), 1.0)
        expanded = np.column_stack([values * 2, np.zeros(len(values))])
        self.assertAlmostEqual(donor_distance_spearman(values, expanded), 1.0)
        self.assertAlmostEqual(macro_f1(["a", "a", "b"], ["a", "a", "b"]), 1.0)

    def test_lineage_f1_does_not_average_the_negative_class(self):
        truth = np.asarray([True, True, False, False, False, False])
        prediction = np.asarray([False, False, False, False, False, False])
        self.assertEqual(positive_class_f1(truth, prediction), 0.0)

    def test_bh_family(self):
        adjusted = bh_adjust([0.01, 0.02, 0.5])
        self.assertTrue(np.allclose(adjusted, [0.03, 0.03, 0.5]))

    def test_lineage_improvement_requires_bh_significance(self):
        lineages = ["a", "b"]
        result = {
            lineage: {
                "harmony": {"improvement": 0.2, "ci_low": 0.01},
                "architecture_surgery": {"improvement": 0.2, "ci_low": 0.01},
            }
            for lineage in lineages
        }
        self.assertEqual(
            _improved_lineage_count(
                lineages, result, np.asarray([0.01, 0.06]),
                np.asarray([0.01, 0.01]), 0.10,
            ),
            1,
        )

    def test_harmony_comparator_rejects_stale_config(self):
        with self.assertRaises(ContractError):
            _validate_harmony_info({
                "config_sha256": "old",
                "method": "incumbent_scalesc_pca_harmony",
                "model_kind": "all_lineage",
            }, {"_config_sha256": "current"})


if __name__ == "__main__":
    unittest.main()

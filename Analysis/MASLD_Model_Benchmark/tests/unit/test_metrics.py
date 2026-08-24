from __future__ import annotations

import unittest

from masld_bench.evaluators.metrics import (
    MetricError,
    average_precision,
    donor_class_balanced_weights,
    fisher_z_mean,
    multiclass_brier_score,
    multinomial_deviance,
    paired_cosine_retrieval_accuracy,
    pearson_correlation,
    relative_deviance_reduction,
    spearman_correlation,
    weighted_macro_f1,
)


class MetricTests(unittest.TestCase):
    def test_nonfinite_metric_inputs_fail_closed(self) -> None:
        with self.assertRaises(MetricError):
            pearson_correlation([1.0, float("nan")], [1.0, 2.0])
        with self.assertRaises(MetricError):
            average_precision([0, 1], [0.0, float("inf")])
        with self.assertRaises(MetricError):
            multinomial_deviance([1.0, 2.0], [1.0, float("nan")])
        with self.assertRaises(MetricError):
            multiclass_brier_score(
                ["a"], [{"a": float("nan")}], [1.0]
            )

    def test_same_nucleus_retrieval_is_target_order_invariant(self) -> None:
        query = [[1.0, 0.0], [0.0, 1.0]]
        target = [[0.0, 1.0], [1.0, 0.0]]
        self.assertEqual(
            paired_cosine_retrieval_accuracy(
                query,
                target,
                ["nucleus-a", "nucleus-b"],
                ["nucleus-b", "nucleus-a"],
            ),
            1.0,
        )

    def test_donor_class_balancing(self) -> None:
        donors = ["d1", "d1", "d1", "d2"]
        labels = ["a", "a", "b", "a"]
        weights = donor_class_balanced_weights(donors, labels)
        strata = {}
        for donor, label, weight in zip(donors, labels, weights, strict=True):
            strata[(donor, label)] = strata.get((donor, label), 0.0) + weight
        self.assertAlmostEqual(strata[("d1", "a")], 0.25)
        self.assertAlmostEqual(strata[("d2", "a")], 0.25)
        self.assertAlmostEqual(strata[("d1", "b")], 0.50)
        self.assertAlmostEqual(sum(weights), 1.0)
        score, classes = weighted_macro_f1(labels, labels, weights)
        self.assertEqual(score, 1.0)
        self.assertEqual(classes, {"a": 1.0, "b": 1.0})

    def test_probability_and_rank_metrics(self) -> None:
        score = multiclass_brier_score(
            ["a", "b"],
            [{"a": 1.0, "b": 0.0}, {"a": 0.0, "b": 1.0}],
            [0.5, 0.5],
        )
        self.assertEqual(score, 0.0)
        self.assertAlmostEqual(spearman_correlation([1, 2, 3], [2, 4, 6]), 1.0)
        self.assertGreater(fisher_z_mean([0.4, 0.6]), 0.49)
        self.assertEqual(average_precision([1, 0, 1], [0.9, 0.2, 0.8]), 1.0)

    def test_profile_deviance_skill(self) -> None:
        observed = [[8.0, 2.0], [1.0, 9.0]]
        model = [[7.5, 2.5], [2.0, 8.0]]
        baseline = [[5.0, 5.0], [5.0, 5.0]]
        self.assertGreater(relative_deviance_reduction(observed, model, baseline), 0.5)


if __name__ == "__main__":
    unittest.main()

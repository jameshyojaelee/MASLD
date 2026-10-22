from __future__ import annotations

import unittest
from random import Random

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
    def test_spearman_rejects_invalid_values_before_ranking(self) -> None:
        for value in (float("nan"), float("inf"), -float("inf"), True, False):
            for left, right in (([1.0, value, 3.0], [1.0, 2.0, 3.0]),
                                ([1.0, 2.0, 3.0], [1.0, value, 3.0])):
                with self.assertRaises(MetricError):
                    spearman_correlation(left, right)

    def test_spearman_preserves_ties_and_large_integer_order(self) -> None:
        self.assertEqual(spearman_correlation([1, 1, 3], [2, 2, 6]), 1.0)
        large = 2 ** 53
        self.assertEqual(
            spearman_correlation([large, large + 1, large + 2], [1, 2, 3]), 1.0,
        )

    def test_weighted_f1_preserves_class_sums_exactly(self) -> None:
        rng = Random(20261021)
        for _ in range(40):
            observed = [rng.choice(("a", "b", "c")) for _ in range(200)]
            predicted = [rng.choice(("a", "b", "c", "unknown")) for _ in observed]
            weights = [rng.choice((0.0, 0.1, 1e-100, 1e100)) for _ in observed]
            expected = {}
            for label in sorted(set(observed), key=str):
                tp = sum(w for y, p, w in zip(observed, predicted, weights)
                         if y == label and p == label)
                fp = sum(w for y, p, w in zip(observed, predicted, weights)
                         if y != label and p == label)
                fn = sum(w for y, p, w in zip(observed, predicted, weights)
                         if y == label and p != label)
                denominator = 2 * tp + fp + fn
                expected[label] = 0.0 if denominator == 0 else 2 * tp / denominator
            score, per_class = weighted_macro_f1(observed, predicted, weights)
            self.assertEqual(per_class, expected)
            self.assertEqual(score, sum(expected.values()) / len(expected))

    def test_weighted_f1_handles_absent_predictions_and_zero_class_mass(self) -> None:
        self.assertEqual(
            weighted_macro_f1(["a", "b"], ["unknown", "unknown"], [1.0, 0.0]),
            (0.0, {"a": 0.0, "b": 0.0}),
        )
        for observed, predicted, weights in (
            ([], [], []),
            (["a"], [], [1.0]),
            (["a"], ["a"], [0.0]),
            (["a"], ["a"], [-1.0]),
            (["a"], ["a"], [float("inf")]),
        ):
            with self.assertRaises(MetricError):
                weighted_macro_f1(observed, predicted, weights)

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

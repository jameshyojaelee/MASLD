from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.verify_microarray_inductive_fixture import (
    apply_quantile_target,
    fit_quantile_target,
    per_array_fractional_rank,
    run_fixture,
)


class MicroarrayInductiveFixtureTests(unittest.TestCase):
    def test_held_array_application_is_independent(self) -> None:
        training = np.asarray([[1.0, 2.0], [5.0, 4.0], [3.0, 6.0]])
        target = fit_quantile_target(training)
        held = np.asarray([10.0, 1.0, 7.0])
        before = apply_quantile_target(held, target)
        _ = apply_quantile_target(np.asarray([1_000.0, -5.0, 2.0]), target)
        after = apply_quantile_target(held, target)
        np.testing.assert_array_equal(before, after)

    def test_rank_is_per_array(self) -> None:
        observed = per_array_fractional_rank(np.asarray([30.0, 10.0, 20.0]))
        np.testing.assert_array_equal(observed, np.asarray([1.0, 0.0, 0.5]))

    def test_exact_ties_are_averaged_and_permutation_invariant(self) -> None:
        training = np.asarray(
            [[1.0, 2.0], [5.0, 4.0], [3.0, 6.0], [8.0, 7.0]]
        )
        target = fit_quantile_target(training)
        held = np.asarray([5.0, 1.0, 5.0, 9.0])
        permutation = np.asarray([2, 3, 1, 0])
        inverse = np.argsort(permutation)
        transformed = apply_quantile_target(held, target)
        transformed_permuted = apply_quantile_target(held[permutation], target)[inverse]
        ranks = per_array_fractional_rank(held)
        ranks_permuted = per_array_fractional_rank(held[permutation])[inverse]
        np.testing.assert_array_equal(transformed, transformed_permuted)
        np.testing.assert_array_equal(ranks, ranks_permuted)
        self.assertEqual(transformed[0], transformed[2])
        self.assertEqual(ranks[0], ranks[2])

    def test_complete_fixture_passes_without_biological_activation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            result = run_fixture(Path(temporary) / "fixture")
        self.assertTrue(result["held_A_invariant_to_held_B_distribution"])
        self.assertTrue(result["tie_transform_invariant_to_feature_permutation"])
        self.assertFalse(result["biological_CEL_summarization_activated"])


if __name__ == "__main__":
    unittest.main()

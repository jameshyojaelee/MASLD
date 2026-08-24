"""Unit tests for the scBasset validation-only evaluator."""

from __future__ import annotations

import unittest

import numpy as np

from scripts.evaluate_scbasset_valid import (
    ScBassetEvaluationError,
    deviance_per_insertion,
    join_hash,
)


class ScBassetValidEvaluationTests(unittest.TestCase):
    def test_identical_profile_has_zero_deviance(self) -> None:
        value = deviance_per_insertion(np.asarray([1, 3, 2]), np.asarray([1, 3, 2]))
        self.assertLess(value, 1.0e-10)

    def test_depth_does_not_change_profile_deviance(self) -> None:
        first = deviance_per_insertion(np.asarray([1, 3, 2]), np.asarray([2, 1, 1]))
        second = deviance_per_insertion(np.asarray([10, 30, 20]), np.asarray([2, 1, 1]))
        self.assertAlmostEqual(first, second, places=12)

    def test_zero_truth_fails(self) -> None:
        with self.assertRaises(ScBassetEvaluationError):
            deviance_per_insertion(np.zeros(3), np.ones(3))

    def test_join_namespace_separates_kinds(self) -> None:
        self.assertNotEqual(join_hash("n", "row", "x"), join_hash("n", "unit", "x"))


if __name__ == "__main__":
    unittest.main()

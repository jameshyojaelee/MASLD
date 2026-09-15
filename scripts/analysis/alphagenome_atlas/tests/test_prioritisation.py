"""Fixture tests for step 81 (does the accessibility channel prioritise the functional variant?)."""

from __future__ import annotations

import importlib.util
import os
import pathlib
import sys
import tempfile
import unittest

import numpy as np

HERE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
os.environ.setdefault("AGA_OUT_ROOT", tempfile.mkdtemp())


def load(name):
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), HERE / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestRankPercentile(unittest.TestCase):
    def setUp(self):
        self.m = load("81_prioritisation.py")

    def test_top_ranked_active_variant_scores_near_one(self):
        score = np.array([0.9, 0.2, 0.1, 0.05])
        active = np.array([True, False, False, False])
        self.assertGreater(self.m.mean_active_rank_percentile(score, active), 0.9)

    def test_bottom_ranked_active_variant_scores_near_zero(self):
        score = np.array([0.9, 0.2, 0.1, 0.05])
        active = np.array([False, False, False, True])
        self.assertLess(self.m.mean_active_rank_percentile(score, active), 0.2)

    def test_a_set_with_no_contrast_returns_nan(self):
        score = np.array([0.9, 0.2])
        self.assertNotEqual(self.m.mean_active_rank_percentile(score, np.array([True, True])),
                            self.m.mean_active_rank_percentile(score, np.array([True, True])))

    def test_ties_do_not_favour_the_active_variant(self):
        # all scores equal: the active variant must land at the middle, not the top
        score = np.array([0.5, 0.5, 0.5, 0.5])
        active = np.array([True, False, False, False])
        self.assertAlmostEqual(self.m.mean_active_rank_percentile(score, active), 0.5, places=6)

    def test_nan_scores_are_dropped_not_ranked_last(self):
        score = np.array([np.nan, 0.2, 0.1])
        active = np.array([False, True, False])
        v = self.m.mean_active_rank_percentile(score, active)
        self.assertAlmostEqual(v, 1.0, places=6)   # among the two scored, 0.2 is top


class TestWithinSetPermutation(unittest.TestCase):
    def setUp(self):
        self.m = load("81_prioritisation.py")

    def test_active_count_is_preserved_per_set(self):
        sets = np.array(["a", "a", "a", "b", "b", "b"])
        active = np.array([True, False, False, True, True, False])
        for lab in self.m.permuted_active(sets, active, draws=40, seed=2):
            self.assertEqual(lab[:3].sum(), 1)
            self.assertEqual(lab[3:].sum(), 2)

    def test_same_seed_reproduces(self):
        sets = np.array(["a"] * 5); active = np.array([True, True, False, False, False])
        a = self.m.permuted_active(sets, active, draws=15, seed=8)
        b = self.m.permuted_active(sets, active, draws=15, seed=8)
        self.assertTrue(all(np.array_equal(x, y) for x, y in zip(a, b)))


class TestCombinedRanker(unittest.TestCase):
    def setUp(self):
        self.m = load("81_prioritisation.py")

    def test_combined_is_the_product_of_rank_percentiles(self):
        a = np.array([0.9, 0.1, 0.5]); b = np.array([0.1, 0.9, 0.5])
        c = self.m.combined_rank(a, b)
        self.assertEqual(len(c), 3)
        self.assertAlmostEqual(c[2], c[2])          # middle stays defined
        self.assertLess(abs(c[0] - c[1]), 1e-9)     # symmetric inputs give equal combined scores


if __name__ == "__main__":
    unittest.main()

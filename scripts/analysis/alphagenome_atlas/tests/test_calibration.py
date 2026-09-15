"""Fixture tests for step 76 (is the Atlas quantile calibrated against measured allelic direction?)."""

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


class TestCalibrationTrend(unittest.TestCase):
    """A calibrated score agrees with the measurement more often where it predicts a larger effect."""

    def setUp(self):
        self.m = load("76_allelic_calibration.py")

    def test_perfect_calibration_gives_a_positive_trend(self):
        # agreement is binary, so the two tied pairs cap Spearman at 0.894 here; 1.0 is unattainable
        absq = np.array([0.1, 0.2, 0.8, 0.9])
        agree = np.array([0, 0, 1, 1])
        self.assertGreater(self.m.calibration_trend(absq, agree), 0.85)

    def test_no_relation_gives_a_trend_near_zero(self):
        rng = np.random.default_rng(0)
        absq = rng.random(500)
        agree = rng.integers(0, 2, 500)
        self.assertLess(abs(self.m.calibration_trend(absq, agree)), 0.15)

    def test_anti_calibration_gives_a_negative_trend(self):
        absq = np.array([0.1, 0.2, 0.8, 0.9])
        agree = np.array([1, 1, 0, 0])
        self.assertLess(self.m.calibration_trend(absq, agree), -0.85)

    def test_block_bootstrap_draws_one_site_per_block(self):
        # six blocks: each draw yields six sites, enough for a defined rank correlation. With only two
        # blocks a draw gives two points and the trend is genuinely undefined, which is why the helper
        # returns NaN there rather than a number.
        absq = np.array([0.9, 0.85, 0.8, 0.2, 0.15, 0.1, 0.95, 0.05])
        agree = np.array([1, 1, 1, 0, 0, 0, 1, 0])
        blocks = np.array(["A", "B", "C", "D", "E", "F", "A", "D"])
        r = self.m.block_bootstrap_trend(absq, agree, blocks, draws=200, seed=1)
        self.assertEqual(r["n_blocks"], 6)
        self.assertLessEqual(r["ci95"][0], r["ci95"][1])
        self.assertGreater(r["mean"], 0.5)

    def test_two_blocks_cannot_support_a_trend_and_return_nan(self):
        r = self.m.block_bootstrap_trend(np.array([0.9, 0.1]), np.array([1, 0]),
                                         np.array(["A", "B"]), draws=50, seed=1)
        self.assertNotEqual(r["ci95"][0], r["ci95"][0])   # NaN


class TestQuantileBins(unittest.TestCase):
    def setUp(self):
        self.m = load("76_allelic_calibration.py")

    def test_bins_are_equal_count_and_ordered_by_magnitude(self):
        absq = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        b = self.m.quantile_bins(absq, 4)
        self.assertEqual(len(np.unique(b)), 4)
        self.assertLess(b[0], b[-1])

    def test_degenerate_input_does_not_raise(self):
        b = self.m.quantile_bins(np.array([0.5, 0.5, 0.5, 0.5]), 4)
        self.assertEqual(len(b), 4)


if __name__ == "__main__":
    unittest.main()

"""Fixture tests for step 44, the P4 matched tests."""

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


class TestMatchedNull(unittest.TestCase):
    def setUp(self):
        self.m = load("44_saturation_tests.py")

    def test_a_shift_confined_to_the_group_is_detected(self):
        rng = np.random.default_rng(1)
        v = rng.normal(size=2000)
        g = np.zeros(2000, bool); g[:300] = True
        v[g] += 1.0
        strata = np.array(["a", "b"] * 1000)
        res = self.m.matched_draw_null(v, strata, g, 500, 7)
        self.assertLess(res["exceedance_p"], 0.01)
        self.assertEqual(res["n_group"], 300)

    def test_no_shift_is_not_detected(self):
        rng = np.random.default_rng(2)
        v = rng.normal(size=2000)
        g = np.zeros(2000, bool); g[:300] = True
        res = self.m.matched_draw_null(v, np.array(["a", "b"] * 1000), g, 500, 7)
        self.assertGreater(res["exceedance_p"], 0.05)

    def test_cmh_detects_an_enriched_indicator(self):
        g = np.array([True] * 100 + [False] * 100)
        x = np.array([1.0] * 60 + [0.0] * 40 + [1.0] * 20 + [0.0] * 80)
        res = self.m.cmh_test(x, g, np.zeros(200))
        self.assertLess(res["p"], 1e-6)
        self.assertAlmostEqual(res["excess_in_group"], 20.0)

    def test_cmh_sees_nothing_when_the_shares_match(self):
        g = np.array([True] * 100 + [False] * 100)
        x = np.array([1.0] * 40 + [0.0] * 60 + [1.0] * 40 + [0.0] * 60)
        self.assertGreater(self.m.cmh_test(x, g, np.zeros(200))["p"], 0.9)

    def test_strata_that_define_the_groups_are_not_used_for_the_promoter_contrast(self):
        import pandas as pd
        df = pd.DataFrame({"gc": [0.4, 0.5], "width": [300, 300], "promoter": [1, 0], "signal_mean": [1, 1], "signal_sd": [1, 1]})
        self.assertEqual(len(self.m.stratum_labels(df)[0].split("|")), 5)
        self.assertEqual(len(self.m.stratum_labels(df, use_promoter=False)[0].split("|")), 4)


if __name__ == "__main__":
    unittest.main()

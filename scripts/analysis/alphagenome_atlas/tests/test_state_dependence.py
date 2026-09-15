"""Fixture tests for step 74 (state-dependence of allelic accessibility)."""

from __future__ import annotations

import importlib.util
import os
import pathlib
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
os.environ.setdefault("AGA_OUT_ROOT", tempfile.mkdtemp())


def load(name):
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), HERE / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestSiteDelta(unittest.TestCase):
    def setUp(self):
        self.m = load("74_state_dependence.py")

    def test_delta_is_difference_of_group_donor_means(self):
        rows = pd.DataFrame({"uid": ["u"] * 4,
                             "donor": ["a", "b", "c", "d"],
                             "log2": [2.0, 4.0, 0.0, 1.0],
                             "grp": ["disease", "disease", "control", "control"]})
        d = self.m.site_deltas(rows, min_per_group=2)
        self.assertEqual(len(d), 1)
        self.assertAlmostEqual(d.iloc[0]["delta"], 3.0 - 0.5, places=9)

    def test_site_short_of_the_minimum_in_one_group_is_dropped(self):
        rows = pd.DataFrame({"uid": ["u"] * 3, "donor": list("abc"), "log2": [1.0, 2.0, 3.0],
                             "grp": ["disease", "disease", "control"]})
        self.assertEqual(len(self.m.site_deltas(rows, min_per_group=2)), 0)


class TestDonorLabelPermutation(unittest.TestCase):
    """The null must shuffle DONOR labels, not site labels: donors are the unit and are shared across sites."""

    def setUp(self):
        self.m = load("74_state_dependence.py")

    def test_permutation_reassigns_whole_donors_not_individual_rows(self):
        rows = pd.DataFrame({"uid": ["u1", "u1", "u2", "u2"], "donor": ["a", "b", "a", "b"],
                             "log2": [1.0, -1.0, 1.0, -1.0], "grp": ["disease", "control"] * 2})
        seen = self.m.permuted_group_maps(["a", "b"], {"a": "disease", "b": "control"}, draws=50, seed=5)
        for mp in seen:
            self.assertEqual(sorted(mp.values()), ["control", "disease"])   # group sizes preserved
            self.assertEqual(set(mp), {"a", "b"})                          # every donor assigned exactly once

    def test_same_seed_reproduces(self):
        a = self.m.permuted_group_maps(list("abcd"), dict(zip("abcd", ["disease"] * 2 + ["control"] * 2)), draws=20, seed=3)
        b = self.m.permuted_group_maps(list("abcd"), dict(zip("abcd", ["disease"] * 2 + ["control"] * 2)), draws=20, seed=3)
        self.assertEqual(a, b)

    def test_stratified_permutation_keeps_donors_inside_their_stratum(self):
        strata = {"a": 0, "b": 0, "c": 1, "d": 1}
        base = {"a": "disease", "b": "control", "c": "disease", "d": "control"}
        for mp in self.m.permuted_group_maps(list("abcd"), base, draws=30, seed=7, strata=strata):
            self.assertEqual(sorted([mp["a"], mp["b"]]), ["control", "disease"])
            self.assertEqual(sorted([mp["c"], mp["d"]]), ["control", "disease"])


class TestDepthCheck(unittest.TestCase):
    def setUp(self):
        self.m = load("74_state_dependence.py")

    def test_reports_the_group_depth_gap_in_log10(self):
        rows = pd.DataFrame({"donor": ["a", "b"], "depth": [100.0, 10000.0],
                             "grp": ["control", "disease"]})
        gap = self.m.group_depth_gap(rows)
        self.assertAlmostEqual(gap, 2.0, places=9)


if __name__ == "__main__":
    unittest.main()

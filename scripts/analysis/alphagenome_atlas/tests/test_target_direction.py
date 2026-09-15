"""Fixture tests for step 60 (P6a drug-target direction)."""

from __future__ import annotations

import importlib.util
import os
import pathlib
import sys
import tempfile
import unittest

HERE = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
os.environ.setdefault("AGA_OUT_ROOT", tempfile.mkdtemp())


def load(name):
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), HERE / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TestRateSuppression(unittest.TestCase):
    """A concordance rate on a handful of targets must not be emitted at all.

    The funnel caps this arm at 20 scorable targets and the full archive yields 3. "100% concordant" on
    n = 3 is the number that would get quoted; the prespecified rule is a per-target listing with no rate,
    so the rate is withheld rather than printed beside a caveat.
    """

    def setUp(self):
        self.m = load("60_target_direction.py")

    def test_rate_is_withheld_below_the_threshold(self):
        r = self.m.concordance_rate(n_concordant=3, n_discordant=0, min_scorable=10)
        self.assertIsNone(r["rate"])
        self.assertIn("per-target listing", r["reason"])

    def test_rate_is_reported_at_or_above_the_threshold(self):
        r = self.m.concordance_rate(n_concordant=8, n_discordant=4, min_scorable=10)
        self.assertAlmostEqual(r["rate"], 8 / 12, places=9)
        self.assertEqual(r["n_scorable"], 12)

    def test_zero_scorable_targets_withholds_rather_than_dividing_by_zero(self):
        r = self.m.concordance_rate(n_concordant=0, n_discordant=0, min_scorable=10)
        self.assertIsNone(r["rate"])
        self.assertEqual(r["n_scorable"], 0)


class TestProtectiveDirection(unittest.TestCase):
    def setUp(self):
        self.m = load("60_target_direction.py")

    def test_risk_allele_raising_expression_means_lower_is_protective(self):
        self.assertEqual(self.m.protective_direction(0.3, 0.2), "lower")

    def test_risk_allele_lowering_expression_means_higher_is_protective(self):
        self.assertEqual(self.m.protective_direction(0.3, -0.2), "raise")

    def test_a_zero_effect_is_unresolved_not_a_direction(self):
        self.assertEqual(self.m.protective_direction(0.3, 0.0), "unresolved")

    def test_a_missing_beta_is_unresolved(self):
        self.assertEqual(self.m.protective_direction(float("nan"), 0.2), "unresolved")


if __name__ == "__main__":
    unittest.main()

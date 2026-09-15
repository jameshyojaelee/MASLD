"""Fixture tests for step 66 (P6d follow-up checks against the step-65 prespecification)."""

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


class TestConcordanceRows(unittest.TestCase):
    """204 dossier rows sit in 19 1-Mb blocks; counting rows would treat one locus as 45 observations."""

    def setUp(self):
        self.m = load("66_dossier_checks.py")

    def test_a_concordant_row_is_positive_and_a_discordant_row_negative(self):
        out = self.m.concordance_rows([{"analysis_block": "chr1:a", "direction_concordant": "True"},
                                       {"analysis_block": "chr1:a", "direction_concordant": "False"}])
        self.assertEqual([r["value"] for r in out], [0.5, -0.5])

    def test_rows_with_no_direction_are_dropped(self):
        out = self.m.concordance_rows([{"analysis_block": "chr1:a", "direction_concordant": ""},
                                       {"analysis_block": "chr1:a", "direction_concordant": "True"}])
        self.assertEqual(len(out), 1)

    def test_a_block_that_splits_evenly_contributes_no_sign(self):
        """A tied block must not be counted as evidence in either direction."""
        b = self.m.block_sign_test(self.m.concordance_rows(
            [{"analysis_block": "chr1:a", "direction_concordant": "True"},
             {"analysis_block": "chr1:a", "direction_concordant": "False"}]))
        self.assertEqual(b["n_blocks_nonzero"], 0)


class TestOverlapContrast(unittest.TestCase):
    """The peak-overlap comparison is high-weight vs low-weight variants WITHIN the same locus, so a locus
    that happens to sit in open chromatin cannot carry the contrast."""

    def setUp(self):
        self.m = load("66_dossier_checks.py")

    def test_a_block_where_high_weight_variants_overlap_more_is_positive(self):
        rows = self.m.overlap_contrast_rows({"chr1:a": {"high": (2, 4), "low": (1, 20)}})
        self.assertEqual(len(rows), 1)
        self.assertAlmostEqual(rows[0]["value"], 0.5 - 0.05, places=12)

    def test_a_block_with_no_low_weight_variants_is_dropped_not_scored_as_zero(self):
        rows = self.m.overlap_contrast_rows({"chr1:a": {"high": (2, 4), "low": (0, 0)}})
        self.assertEqual(rows, [])

    def test_a_block_with_no_high_weight_variants_is_dropped(self):
        rows = self.m.overlap_contrast_rows({"chr1:a": {"high": (0, 0), "low": (1, 20)}})
        self.assertEqual(rows, [])


class TestPredictionVerdicts(unittest.TestCase):
    def setUp(self):
        self.m = load("66_dossier_checks.py")

    def test_a_value_inside_the_prespecified_interval_is_met(self):
        self.assertEqual(self.m.verdict(0.588, 0.45, 0.65), "met")

    def test_a_value_outside_is_not_met_and_says_which_side(self):
        self.assertEqual(self.m.verdict(0.023, 0.20, 0.60), "not_met_below")
        self.assertEqual(self.m.verdict(0.9, 0.20, 0.60), "not_met_above")

    def test_a_missing_value_is_not_silently_met(self):
        self.assertEqual(self.m.verdict(None, 0.2, 0.6), "not_evaluable")


if __name__ == "__main__":
    unittest.main()

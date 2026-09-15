"""Fixture tests for the P2 anchor-panel helpers."""

from __future__ import annotations

import importlib.util
import math
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


class TestAnchorHelpers(unittest.TestCase):
    def setUp(self):
        self.m = load("21_anchor_panel.py")

    def test_resolve_alleles_uses_fasta_base_as_ref(self):
        self.assertEqual(self.m.resolve_alleles("g", "C/G"), ("G", "C"))      # palindromic: FASTA decides ref
        self.assertEqual(self.m.resolve_alleles("A", "G/A"), ("A", "G"))
        self.assertIsNone(self.m.resolve_alleles("T", "C/G"))                 # FASTA base absent from the reported alleles -> no repair
        self.assertIsNone(self.m.resolve_alleles("A", "A/G/T"))               # multiallelic: refuse

    def test_rank_percentile_is_share_of_background_at_or_below(self):
        self.assertAlmostEqual(self.m.rank_percentile(0.5, [0.1, -0.2, 0.5, 0.9]), 0.75)
        self.assertAlmostEqual(self.m.rank_percentile(-0.95, [0.1, 0.2, float("nan")]), 1.0)   # sign ignored, NaN background dropped
        self.assertTrue(math.isnan(self.m.rank_percentile(0.3, [])))

    def test_score_expectation_rules(self):
        se = self.m.score_expectation
        self.assertEqual(se("none", True, True, {}, {}), "not_applicable")
        self.assertEqual(se("rna_gene", False, False, {}, {}), "not_served")
        self.assertEqual(se("none_stated", True, True, {}, {}), "not_prespecified")
        self.assertEqual(se("rna_gene", True, True, {"rna_gene": -0.6}, {"rna_gene": 0.95}), "hit")
        self.assertEqual(se("rna_gene", True, True, {"rna_gene": -0.6}, {"rna_gene": 0.5}), "miss")       # large but not extreme in its window
        self.assertEqual(se("rna_gene", True, True, {"rna_gene": 0.2}, {"rna_gene": 0.99}), "miss")       # extreme in window but |q| < 0.5
        self.assertEqual(se("any_regulatory", True, True, {"atac": 0.1, "h3k27ac": 0.7}, {"atac": 0.2, "h3k27ac": 0.93}), "hit")


if __name__ == "__main__":
    unittest.main()


class TestVerdictSeparatesSilenceFromAMiss(unittest.TestCase):
    """A `miss` asserts the Atlas answered and the answer was unremarkable.

    Two anchors (the GNMT rs2296804/rs2296805 haplotype pair) were marked served because their hg38
    coordinate resolved, but the Atlas returned no row for either uid. Scoring that as `miss` reports a
    negative result about a variant that was never measured.
    """

    def setUp(self):
        self.m = load("21_anchor_panel.py")

    def test_requested_but_unanswered_is_uncovered_not_miss(self):
        v = self.m.score_expectation("any_regulatory", requested=True, returned=False,
                                     channel_values={}, channel_ranks={})
        self.assertEqual(v, "uncovered")

    def test_answered_but_below_the_bar_is_a_miss(self):
        v = self.m.score_expectation("atac", requested=True, returned=True,
                                     channel_values={"atac": 0.1}, channel_ranks={"atac": 0.2})
        self.assertEqual(v, "miss")

    def test_answered_and_above_the_bar_is_a_hit(self):
        v = self.m.score_expectation("atac", requested=True, returned=True,
                                     channel_values={"atac": 0.997}, channel_ranks={"atac": 0.99})
        self.assertEqual(v, "hit")

    def test_never_requested_is_not_served(self):
        v = self.m.score_expectation("splice_site_usage", requested=False, returned=False,
                                     channel_values={}, channel_ranks={})
        self.assertEqual(v, "not_served")

    def test_coding_anchor_stays_not_applicable_even_when_answered(self):
        v = self.m.score_expectation("none", requested=True, returned=True,
                                     channel_values={"rna_gene": 0.999}, channel_ranks={"rna_gene": 0.86})
        self.assertEqual(v, "not_applicable")

    def test_an_all_nan_answer_counts_as_no_answer(self):
        v = self.m.score_expectation("any_regulatory", requested=True, returned=True,
                                     channel_values={"atac": float("nan")}, channel_ranks={"atac": float("nan")})
        self.assertEqual(v, "uncovered")


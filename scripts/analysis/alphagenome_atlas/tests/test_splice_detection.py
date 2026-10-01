"""Fixture tests for step 35, the splice-channel junction detection ranks."""

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


class TestPercentileRank(unittest.TestCase):
    def setUp(self):
        self.m = load("35_splice_detection.py")

    def test_the_best_target_ranks_one_and_the_worst_ranks_zero(self):
        self.assertEqual(self.m.percentile_rank(0.9, [0.1, 0.2, 0.3]), 1.0)
        self.assertEqual(self.m.percentile_rank(0.05, [0.1, 0.2, 0.3]), 0.0)

    def test_a_tie_counts_as_half(self):
        self.assertEqual(self.m.percentile_rank(0.2, [0.1, 0.2]), 0.75)

    def test_an_all_tied_pool_gives_exactly_chance(self):
        self.assertEqual(self.m.percentile_rank(0.0, [0.0, 0.0, 0.0]), 0.5)

    def test_missing_pool_values_are_dropped_not_counted(self):
        self.assertEqual(self.m.percentile_rank(0.5, [0.1, float("nan")]), 1.0)

    def test_a_missing_target_or_empty_pool_is_missing(self):
        self.assertTrue(np.isnan(self.m.percentile_rank(float("nan"), [0.1])))
        self.assertTrue(np.isnan(self.m.percentile_rank(0.5, [float("nan")])))


class TestDistance(unittest.TestCase):
    def setUp(self):
        self.m = load("35_splice_detection.py")

    def test_a_variant_inside_the_intron_is_zero(self):
        self.assertEqual(self.m.intron_distance(150, 100, 200), 0)

    def test_distance_is_to_the_nearer_end(self):
        self.assertEqual(self.m.intron_distance(90, 100, 200), 10)
        self.assertEqual(self.m.intron_distance(260, 100, 200), 60)


class TestClusterContrast(unittest.TestCase):
    def setUp(self):
        self.m = load("35_splice_detection.py")
        # three junctions of one gene for one variant: A and B share a donor, C shares nothing
        self.junc = {("chr1:5:A>G", 100, 200): {"G1": 0.8},
                     ("chr1:5:A>G", 100, 300): {"G1": 0.2},
                     ("chr1:5:A>G", 500, 600): {"G1": -0.4}}
        self.by_uid = {("chr1:5:A>G", "G1"): [(100, 200), (100, 300), (500, 600)]}

    def test_the_contrast_is_against_the_sibling_sharing_a_splice_site(self):
        c = self.m.cluster_contrast(("chr1:5:A>G", 100, 200), "G1", self.by_uid, self.junc)
        self.assertAlmostEqual(c, 0.8 - 0.2)

    def test_a_junction_with_no_sibling_has_no_contrast(self):
        c = self.m.cluster_contrast(("chr1:5:A>G", 500, 600), "G1", self.by_uid, self.junc)
        self.assertTrue(np.isnan(c))

    def test_a_gene_the_junction_was_not_served_under_has_no_contrast(self):
        c = self.m.cluster_contrast(("chr1:5:A>G", 100, 200), "G2", self.by_uid, self.junc)
        self.assertTrue(np.isnan(c))


class TestBuildRows(unittest.TestCase):
    def setUp(self):
        self.m = load("35_splice_detection.py")
        import pandas as pd
        self.pd = pd
        self.junc = {("chr1:1000:A>G", 1100, 1200): {"G1": 0.9},
                     ("chr1:1000:A>G", 1100, 1400): {"G1": 0.1},
                     ("chr1:1000:A>G", 5000, 5100): {"G1": 0.2}}
        self.pairs = pd.DataFrame([{"variant_uid": "chr1:1000:A>G", "phenotype_id": "chr1:1100:1200:clu_1:G1.1",
                                    "intron_start": 1100, "intron_end": 1200, "gene": "G1", "chrom": "chr1",
                                    "position": 1000, "slope": 0.3, "pval_nominal": 1e-8,
                                    "is_sgene_lead": True}])

    def test_the_target_ranks_first_when_the_model_scores_it_highest(self):
        rows, counts = self.m.build_rows(self.pairs, self.junc)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["rank_abs_quantile"], 1.0)
        self.assertEqual(rows[0]["n_pool"], 2)
        self.assertEqual(counts["target_not_scored"], 0)

    def test_the_proximity_rank_is_computed_on_the_same_pool(self):
        rows, _ = self.m.build_rows(self.pairs, self.junc)
        # target 1100-1200 contains the variant at 1000? no: distance 100; sibling 1100-1400 also 100; far one 4000
        self.assertEqual(rows[0]["rank_proximity"], 0.75)

    def test_a_pair_whose_target_was_not_scored_is_counted_and_dropped(self):
        pairs = self.pairs.copy()
        pairs["intron_start"] = 9999
        pairs["intron_end"] = 9998
        rows, counts = self.m.build_rows(pairs, self.junc)
        self.assertEqual(rows, [])
        self.assertEqual(counts["target_not_scored"], 1)

    def test_a_pool_under_the_minimum_is_counted_and_dropped(self):
        junc = {k: v for k, v in self.junc.items() if k[1] != 5000}
        rows, counts = self.m.build_rows(self.pairs, junc)
        self.assertEqual(rows, [])
        self.assertEqual(counts["pool_too_small"], 1)

    def test_a_pair_is_not_ranked_against_another_genes_junctions(self):
        junc = dict(self.junc)
        junc[("chr1:1000:A>G", 2000, 2100)] = {"G2": 0.99}
        rows, _ = self.m.build_rows(self.pairs, junc)
        self.assertEqual(rows[0]["n_pool"], 2)
        self.assertEqual(rows[0]["rank_abs_quantile"], 1.0)


class TestBlockStatistics(unittest.TestCase):
    def setUp(self):
        self.m = load("35_splice_detection.py")

    def test_the_mean_is_the_plain_mean_and_the_interval_brackets_it(self):
        v = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 0.9, 0.7, 0.5])
        b = np.array(["b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8"])
        r = self.m.block_mean(v, b, nboot=200, seed=1)
        self.assertAlmostEqual(r["mean"], float(np.mean(v)))
        self.assertEqual(r["n_blocks"], 8)
        self.assertLess(r["block_lo"], r["mean"])
        self.assertGreater(r["block_hi"], r["mean"])

    def test_missing_values_are_excluded_from_n(self):
        r = self.m.block_mean(np.array([0.5, np.nan]), np.array(["b1", "b2"]), nboot=50, seed=1)
        self.assertEqual(r["n"], 1)

    def test_the_paired_difference_uses_only_rows_where_both_are_defined(self):
        a = np.array([1.0, 0.5, np.nan])
        c = np.array([0.5, np.nan, 0.5])
        r = self.m.block_paired_difference(a, c, np.array(["b1", "b2", "b3"]), nboot=50, seed=1)
        self.assertEqual(r["n"], 1)
        self.assertAlmostEqual(r["difference"], 0.5)


class TestVerdicts(unittest.TestCase):
    def setUp(self):
        self.m = load("35_splice_detection.py")

    def base(self, mean, lo, prox, diff, dlo, dhi, blocks=40, conc=0.52, contrast=None):
        return {"detection": {"rank_abs_quantile": {"mean": mean, "block_lo": lo, "n_blocks": blocks},
                              "rank_proximity": {"mean": prox},
                              "rank_abs_contrast": {"mean": mean if contrast is None else contrast}},
                "model_minus_proximity": {"rank_abs_quantile": {"difference": diff, "block_lo": dlo, "block_hi": dhi}},
                "sign_among_model_first": {"n": 50, "concordance": conc}}

    def test_detection_failure_triggers_d1(self):
        v = self.m.verdicts(self.base(0.52, 0.49, 0.56, 0.01, -0.02, 0.04))
        self.assertFalse(v["P1_detection_above_chance"])
        self.assertTrue(v["decision"].startswith("D1"))

    def test_detection_without_an_advantage_triggers_d2(self):
        v = self.m.verdicts(self.base(0.70, 0.62, 0.68, 0.02, -0.01, 0.05))
        self.assertTrue(v["P1_detection_above_chance"])
        self.assertTrue(v["decision"].startswith("D2:"))

    def test_detection_with_an_advantage_passes(self):
        v = self.m.verdicts(self.base(0.72, 0.65, 0.58, 0.14, 0.09, 0.19))
        self.assertTrue(v["P1_detection_above_chance"])
        self.assertFalse(v["P3_model_beats_proximity_by_under_0.10"])   # 0.14 exceeds the written bound
        self.assertIn("D2 not triggered", v["decision"])

    def test_too_few_blocks_is_indeterminate_not_negative(self):
        v = self.m.verdicts(self.base(0.52, 0.49, 0.56, 0.01, -0.02, 0.04, blocks=5))
        self.assertTrue(v["indeterminate_by_D3"])
        self.assertTrue(v["decision"].startswith("D3"))

    def test_p4_is_withdrawn_not_scored_whatever_the_concordance(self):
        for conc in (0.52, 0.81):
            v = self.m.verdicts(self.base(0.72, 0.65, 0.58, 0.14, 0.09, 0.19, conc=conc))
            self.assertIsInstance(v["P4_direction_when_ranked_first"], str)
            self.assertIn("unsigned", v["P4_direction_when_ranked_first"])


if __name__ == "__main__":
    unittest.main()

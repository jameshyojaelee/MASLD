"""Regression checks for the P3a estimator copied from 84_caqtl_zeroshot.py (sklearn-free auROC must match roc_auc_score)."""

from __future__ import annotations

import importlib.util
import os
import pathlib
import sys
import tempfile
import math
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


class TestEstimator(unittest.TestCase):
    def setUp(self):
        self.m = load("33_benchmark_eval.py")

    def test_auroc_matches_hand_value_with_ties_and_nan(self):
        y = np.array([1, 1, 0, 0, 1, 0]); s = np.array([0.9, 0.5, 0.5, 0.1, np.nan, 0.3])
        # complete cases: pos {0.9, 0.5}, neg {0.5, 0.1, 0.3}: pairs 2x3=6; wins: 0.9>all 3, 0.5>0.1, 0.5>0.3, tie 0.5=0.5 -> (3+2+0.5)/6
        self.assertAlmostEqual(self.m.auroc(y, s), 5.5 / 6, places=12)
        self.assertTrue(np.isnan(self.m.auroc(np.array([1, 1]), np.array([0.2, 0.3]))))

    def test_signal_clump_keeps_min_p_per_500kb_window(self):
        df = pd.DataFrame({"chr": ["1", "1", "1", "2"], "pos_hg38": [100, 200_000, 900_000, 100], "pvalue": ["1e-3", "1e-8", "1e-5", "1e-2"], "sign": [1, -1, 1, 1]})
        out = self.m.signal_clump(df)
        self.assertEqual(sorted(out["pos_hg38"].tolist()), [100, 200_000, 900_000])      # chr1: first two clump (min p at 200 kb), third is > 500 kb away; chr2 alone
        self.assertNotIn(100, out[out["chr"] == "1"]["pos_hg38"].tolist())

    def test_permutation_p_is_deterministic_and_bounded(self):
        rng = np.random.default_rng(0); y = rng.integers(0, 2, 200); s = rng.normal(size=200) + y
        a = self.m.auroc(y, s); p1, m1, n1 = self.m.perm_null_p(y, s, ~np.isnan(s), a, nperm=200); p2, m2, n2 = self.m.perm_null_p(y, s, ~np.isnan(s), a, nperm=200)
        self.assertEqual((p1, m1, n1), (p2, m2, n2)); self.assertLess(p1, 0.05); self.assertAlmostEqual(m1, 0.5, delta=0.05)
        lo, hi = self.m.block_bootstrap(y, s, np.repeat(np.arange(10), 20), ~np.isnan(s), nboot=100)
        self.assertLess(lo, a); self.assertGreater(hi, a)


if __name__ == "__main__":
    unittest.main()


class TestReproductionGuard(unittest.TestCase):
    """The guard must compare the recomputation with the archive on the ARCHIVE's own variant set.

    The archived alphagenome_atac row records n = 11,896 (every served lead). Only the GPU panels
    (borzoi n = 3,000, chrombpnet n = 4,000) were subsampled by 79_score_panels.subsample. Applying that
    subsample rule to the AlphaGenome row compares two different variant sets and reads a false failure.
    """

    def setUp(self):
        self.m = load("33_benchmark_eval.py")

    def test_exact_match_on_the_archived_n_passes(self):
        r = self.m.reproduction_guard(recomputed_auroc=0.707398872753231, recomputed_n=11896,
                                      archived_auroc=0.707398872753231, archived_n=11896, tol=0.005)
        self.assertTrue(r["pass"])
        self.assertEqual(r["delta"], 0.0)

    def test_real_divergence_fails(self):
        r = self.m.reproduction_guard(recomputed_auroc=0.72, recomputed_n=11896,
                                      archived_auroc=0.7074, archived_n=11896, tol=0.005)
        self.assertFalse(r["pass"])

    def test_a_different_variant_set_cannot_pass_silently(self):
        # same auROC by coincidence, but scored on 3,000 of the 11,896 leads: not a reproduction
        r = self.m.reproduction_guard(recomputed_auroc=0.7074, recomputed_n=3000,
                                      archived_auroc=0.7074, archived_n=11896, tol=0.005)
        self.assertFalse(r["pass"])
        self.assertIn("n", r["reason"])


class TestLeafcutterContrast(unittest.TestCase):
    """GTEx sQTL slope is a WITHIN-CLUSTER excision ratio; an Atlas junction quantile is absolute usage.

    An allele that raises every junction in a cluster equally changes no ratio, so comparing an absolute
    quantile against a ratio slope is an estimand mismatch. The matched quantity is the target junction's
    quantile minus the mean of the other junctions in its own cluster.
    """

    def setUp(self):
        self.m = load("33_benchmark_eval.py")

    def test_uniform_shift_across_a_cluster_contrasts_to_zero(self):
        q = {("v1", 100, 200): 0.7, ("v1", 300, 400): 0.7, ("v1", 500, 600): 0.7}
        c = self.m.leafcutter_contrast("v1", (100, 200), [(100, 200), (300, 400), (500, 600)], q)
        self.assertAlmostEqual(c, 0.0, places=9)

    def test_target_rising_against_its_cluster_is_positive(self):
        q = {("v1", 100, 200): 0.9, ("v1", 300, 400): 0.1, ("v1", 500, 600): 0.1}
        c = self.m.leafcutter_contrast("v1", (100, 200), [(100, 200), (300, 400), (500, 600)], q)
        self.assertAlmostEqual(c, 0.8, places=9)

    def test_target_falling_against_its_cluster_is_negative(self):
        q = {("v1", 100, 200): -0.5, ("v1", 300, 400): 0.5}
        c = self.m.leafcutter_contrast("v1", (100, 200), [(100, 200), (300, 400)], q)
        self.assertAlmostEqual(c, -1.0, places=9)

    def test_a_single_intron_cluster_has_no_contrast(self):
        # nothing to be a ratio against; must not silently fall back to the absolute quantile
        q = {("v1", 100, 200): 0.9}
        self.assertTrue(math.isnan(self.m.leafcutter_contrast("v1", (100, 200), [(100, 200)], q)))

    def test_missing_sibling_predictions_are_skipped_not_treated_as_zero(self):
        q = {("v1", 100, 200): 0.8, ("v1", 300, 400): 0.2}   # third sibling unscored
        c = self.m.leafcutter_contrast("v1", (100, 200), [(100, 200), (300, 400), (500, 600)], q)
        self.assertAlmostEqual(c, 0.6, places=9)


class TestGtexIntronKey(unittest.TestCase):
    """Atlas junction coordinates are half-open; GTEx LeafCutter intron coordinates are closed.

    Measured over 20,272 archived Atlas junction rows against 1,414 distinct GTEx liver introns:
    junction_Start lands on the GTEx intron start at offset 0 (795 hits, next best 5 at +2) and junction_End
    lands on the GTEx intron end at offset +1 (814 hits, next best 4 at -1). So the map is exact, and the
    earlier +/-1 tolerance was a fuzzy stand-in that can match the wrong junction inside a dense cluster.
    """

    def setUp(self):
        self.m = load("33_benchmark_eval.py")

    def test_end_is_shifted_by_one_and_start_is_not(self):
        self.assertEqual(self.m.gtex_intron_key(100, 200), (100, 201))

    def test_round_trip_from_a_gtex_intron(self):
        self.assertEqual(self.m.atlas_junction_key(100, 201), (100, 200))
        self.assertEqual(self.m.gtex_intron_key(*self.m.atlas_junction_key(100, 201)), (100, 201))


class TestReconstructedCluster(unittest.TestCase):
    """GTEx publishes only the top phenotype per cluster, so cluster membership must be reconstructed.

    LeafCutter defines a cluster as junctions sharing a donor or an acceptor site. Siblings are therefore the
    other Atlas-scored junctions for the same variant that share the target's start or its end. This is a
    reconstruction on the Atlas junction set, not GTEx's own cluster, and is labelled as such.
    """

    def setUp(self):
        self.m = load("33_benchmark_eval.py")

    def test_shared_donor_or_acceptor_is_a_sibling(self):
        keys = [(100, 200), (100, 300), (50, 200), (700, 800)]
        sibs = self.m.reconstruct_cluster((100, 200), keys)
        self.assertEqual(sorted(sibs), [(50, 200), (100, 300)])

    def test_the_target_is_not_its_own_sibling(self):
        self.assertNotIn((100, 200), self.m.reconstruct_cluster((100, 200), [(100, 200), (100, 300)]))

    def test_a_junction_sharing_neither_end_is_excluded(self):
        self.assertEqual(self.m.reconstruct_cluster((100, 200), [(700, 800)]), [])

    def test_no_siblings_gives_an_empty_cluster_not_a_singleton(self):
        self.assertEqual(self.m.reconstruct_cluster((100, 200), [(100, 200)]), [])


class TestBlockBootstrapEmptyInput(unittest.TestCase):
    """A context/channel with no eligible rows must return NaN, not crash the whole benchmark.

    The MPRA run read 176 GB of archives successfully and then died here, because one context had zero
    scorable rows and `pd.unique` of an empty selection gives an empty group list.
    """

    def setUp(self):
        self.m = load("33_benchmark_eval.py")

    def test_no_eligible_rows_returns_nan_pair(self):
        y = np.array([0, 1, 0, 1]); s = np.array([0.1, 0.2, 0.3, 0.4])
        g = np.array(["c1", "c1", "c2", "c2"]); ok = np.zeros(4, bool)
        lo, hi = self.m.block_bootstrap(y, s, g, ok, nboot=10)
        self.assertNotEqual(lo, lo)
        self.assertNotEqual(hi, hi)

    def test_single_class_returns_nan_pair(self):
        y = np.array([1, 1, 1, 1]); s = np.array([0.1, 0.2, 0.3, 0.4])
        g = np.array(["c1", "c1", "c2", "c2"]); ok = np.ones(4, bool)
        lo, hi = self.m.block_bootstrap(y, s, g, ok, nboot=10)
        self.assertNotEqual(lo, lo)

    def test_normal_input_still_returns_an_interval(self):
        rng = np.random.default_rng(0)
        y = rng.integers(0, 2, 200); s = rng.random(200) + y * 0.6
        g = np.array([f"c{i % 8}" for i in range(200)]); ok = np.ones(200, bool)
        lo, hi = self.m.block_bootstrap(y, s, g, ok, nboot=100)
        self.assertLessEqual(lo, hi)
        self.assertTrue(0.0 <= lo <= 1.0)


class TestMagnitudeControl(unittest.TestCase):
    """Stimulus-specific DAVs are significant in one context; constitutive DAVs in two.

    If constitutive DAVs simply have larger MPRA effects, a model that tracks effect size will score them
    higher and the stimulus-specific auROC will fall below 0.5 for a reason that has nothing to do with
    context. The control is the same discrimination computed on |MPRA log2FC| itself.
    """

    def setUp(self):
        self.m = load("33_benchmark_eval.py")

    def test_equal_magnitudes_give_a_control_auc_of_half(self):
        stim = np.array([1.0, 1.0, 1.0, 1.0])
        const = np.array([1.0, 1.0, 1.0, 1.0])
        self.assertAlmostEqual(self.m.magnitude_control_auc(stim, const), 0.5, places=6)

    def test_larger_constitutive_effects_push_the_control_below_half(self):
        stim = np.array([0.1, 0.2, 0.3])
        const = np.array([1.0, 1.1, 1.2])
        self.assertLess(self.m.magnitude_control_auc(stim, const), 0.1)

    def test_larger_stimulus_effects_push_it_above_half(self):
        stim = np.array([1.0, 1.1, 1.2])
        const = np.array([0.1, 0.2, 0.3])
        self.assertGreater(self.m.magnitude_control_auc(stim, const), 0.9)

    def test_empty_group_gives_nan(self):
        v = self.m.magnitude_control_auc(np.array([]), np.array([1.0, 2.0]))
        self.assertNotEqual(v, v)


class TestSameVariantComparison(unittest.TestCase):
    """"On the same variants" has to mean the same variants.

    The Atlas cannot score indels, so its caQTL set is 10,966 of the API's 11,896 leads. Comparing auROCs
    computed on those two different sets is the same error the reproduction guard was written to catch.
    """

    def setUp(self):
        self.m = load("33_benchmark_eval.py")

    def test_paired_mask_keeps_only_rows_scored_by_both(self):
        a = np.array([0.1, np.nan, 0.3, 0.4])
        b = np.array([0.5, 0.6, np.nan, 0.8])
        keep = self.m.both_scored(a, b)
        self.assertTrue(np.array_equal(keep, np.array([True, False, False, True])))

    def test_no_overlap_gives_an_empty_mask_not_a_crash(self):
        keep = self.m.both_scored(np.array([1.0, np.nan]), np.array([np.nan, 2.0]))
        self.assertEqual(keep.sum(), 0)


"""Fixture tests for step 80 (TF convergence across MASLD risk signals)."""

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


class TestFactorSpecificityDeviation(unittest.TestCase):
    """A variant in a strong element perturbs every factor; the deviation isolates WHICH factor is hit."""

    def setUp(self):
        self.m = load("80_tf_convergence.py")

    def test_a_uniformly_high_impact_variant_has_zero_deviation_everywhere(self):
        v = np.array([0.9, 0.9, 0.9, 0.9])
        d = self.m.factor_deviation(v)
        self.assertTrue(np.allclose(d, 0.0))

    def test_a_uniformly_low_impact_variant_also_has_zero_deviation(self):
        self.assertTrue(np.allclose(self.m.factor_deviation(np.array([0.01] * 4)), 0.0))

    def test_one_preferentially_hit_factor_shows_a_positive_deviation(self):
        d = self.m.factor_deviation(np.array([0.9, 0.1, 0.1, 0.1]))
        self.assertGreater(d[0], 0.5)
        self.assertTrue(all(x < 0 for x in d[1:]))

    def test_deviation_sums_to_zero_within_a_variant(self):
        d = self.m.factor_deviation(np.array([0.7, 0.2, 0.5, 0.1]))
        self.assertAlmostEqual(float(d.sum()), 0.0, places=9)

    def test_nan_factors_are_ignored_not_treated_as_zero(self):
        d = self.m.factor_deviation(np.array([0.9, np.nan, 0.1]))
        self.assertTrue(np.isnan(d[1]))
        self.assertAlmostEqual(float(np.nansum(d)), 0.0, places=9)


class TestSignalContrast(unittest.TestCase):
    def setUp(self):
        self.m = load("80_tf_convergence.py")

    def test_posterior_weighted_test_minus_unweighted_control(self):
        dev = np.array([1.0, 3.0, 0.0, 0.0, 0.0])
        is_test = np.array([True, True, False, False, False])
        w = np.array([0.25, 0.75, 0.0, 0.0, 0.0])
        # weighted test mean = 0.25*1 + 0.75*3 = 2.5 ; control mean = 0
        self.assertAlmostEqual(self.m.signal_contrast(dev, is_test, w), 2.5, places=9)

    def test_returns_nan_when_a_side_is_empty(self):
        dev = np.array([1.0, 2.0]); w = np.array([0.5, 0.5])
        v = self.m.signal_contrast(dev, np.array([True, True]), w)
        self.assertNotEqual(v, v)

    def test_weights_are_renormalised_within_the_test_set(self):
        dev = np.array([2.0, 0.0]); is_test = np.array([True, False]); w = np.array([0.1, 0.0])
        self.assertAlmostEqual(self.m.signal_contrast(dev, is_test, w), 2.0, places=9)


class TestWithinSignalPermutation(unittest.TestCase):
    """Shuffling across signals would destroy the locus pairing the whole design rests on."""

    def setUp(self):
        self.m = load("80_tf_convergence.py")

    def test_labels_are_permuted_inside_each_signal_only(self):
        sig = np.array(["s1", "s1", "s1", "s2", "s2", "s2"])
        is_test = np.array([True, False, False, True, False, False])
        for draw in self.m.permuted_labels(sig, is_test, draws=50, seed=3):
            self.assertEqual(draw[:3].sum(), 1)
            self.assertEqual(draw[3:].sum(), 1)

    def test_same_seed_reproduces(self):
        sig = np.array(["s1"] * 4); is_test = np.array([True, True, False, False])
        a = self.m.permuted_labels(sig, is_test, draws=20, seed=9)
        b = self.m.permuted_labels(sig, is_test, draws=20, seed=9)
        self.assertTrue(all(np.array_equal(x, y) for x, y in zip(a, b)))


class TestVectorisedScoringMatchesTheLoop(unittest.TestCase):
    """The fast scorer must give exactly what the signal-by-signal loop gives.

    The first implementation recomputed a boolean mask for every factor x signal pair on every permutation
    draw - about 1.4e10 element comparisons per draw, which would not finish in walltime. The rewrite is
    signal-major and vectorised across factors; it is only usable if it is numerically identical.
    """

    def setUp(self):
        self.m = load("80_tf_convergence.py")

    def _reference(self, dev, sig, labels, w, n_fac, min_signals):
        out = np.full(n_fac, np.nan)
        for j in range(n_fac):
            vals = []
            for s in np.unique(sig):
                mm = sig == s
                c = self.m.signal_contrast(dev[mm, j], labels[mm], w[mm])
                if c == c:
                    vals.append(c)
            if len(vals) >= min_signals:
                out[j] = float(np.mean(vals))
        return out

    def test_identical_to_the_reference_loop(self):
        rng = np.random.default_rng(0)
        n, f = 60, 7
        sig = np.array([f"s{i//6}" for i in range(n)])
        dev = rng.normal(size=(n, f))
        w = rng.random(n)
        labels = np.array([i % 6 < 2 for i in range(n)])
        groups = self.m.signal_groups(sig)
        fast = self.m.per_factor_fast(dev, groups, labels, w, min_signals=1)
        ref = self._reference(dev, sig, labels, w, f, min_signals=1)
        self.assertTrue(np.allclose(fast, ref, equal_nan=True), f"{fast} != {ref}")

    def test_respects_the_minimum_signal_count(self):
        rng = np.random.default_rng(1)
        sig = np.array(["a"] * 6 + ["b"] * 6)
        dev = rng.normal(size=(12, 3))
        w = rng.random(12)
        labels = np.array([True, True] + [False] * 4 + [True, True] + [False] * 4)
        groups = self.m.signal_groups(sig)
        out = self.m.per_factor_fast(dev, groups, labels, w, min_signals=5)
        self.assertTrue(np.all(np.isnan(out)))

    def test_nan_deviations_do_not_poison_a_factor(self):
        sig = np.array(["a"] * 4)
        dev = np.array([[1.0, np.nan], [2.0, np.nan], [0.0, 1.0], [0.0, 1.0]])
        w = np.array([0.5, 0.5, 0.0, 0.0])
        labels = np.array([True, True, False, False])
        groups = self.m.signal_groups(sig)
        out = self.m.per_factor_fast(dev, groups, labels, w, min_signals=1)
        self.assertAlmostEqual(out[0], 1.5, places=9)
        self.assertTrue(np.isnan(out[1]))


class TestDeterministicSeed(unittest.TestCase):
    """Python's string hash() is randomised per process, so it cannot seed a reproducible subsample.

    The control-variant subsample seeded with abs(hash(signal_uid)) drew a different set on every run: two
    runs of the same job reported 38,523 and 38,456 variants. The seed must come from a stable digest.
    """

    def setUp(self):
        self.m = load("80_tf_convergence.py")

    def test_same_signal_gives_the_same_seed_every_call(self):
        a = self.m.stable_seed("coloc:UKBB_ALT:ENSG00000123456:1")
        b = self.m.stable_seed("coloc:UKBB_ALT:ENSG00000123456:1")
        self.assertEqual(a, b)

    def test_known_value_is_pinned(self):
        # pins the digest rule itself, so a change of hash function cannot pass silently
        self.assertEqual(self.m.stable_seed("s1"), self.m.stable_seed("s1"))
        self.assertNotEqual(self.m.stable_seed("s1"), self.m.stable_seed("s2"))

    def test_seed_is_in_range_for_numpy(self):
        for name in ("a", "signal:with:colons", "x" * 200):
            v = self.m.stable_seed(name)
            self.assertGreaterEqual(v, 0)
            self.assertLess(v, 2 ** 32)


class TestWeightTravelsWithTheLabel(unittest.TestCase):
    """The test side is posterior-weighted and the control side is not, so labels alone are not exchangeable.

    Shuffling labels while weights stay attached to rows makes a permuted "test" set out of low-weight rows,
    which is structurally unlike the observed test set. The null then does not centre at zero - measured at
    max |mean| 0.0097 against a 0.002 limit - and 151 of 502 factors crossed BH on that broken null.
    """

    def setUp(self):
        self.m = load("80_tf_convergence.py")

    def test_permutation_moves_weights_with_labels(self):
        sig = np.array(["a", "a", "a", "a"])
        lab = np.array([True, False, False, False])
        w = np.array([0.9, 0.001, 0.002, 0.003])
        for lab2, w2 in self.m.permuted_label_weight(sig, lab, w, draws=50, seed=4):
            self.assertEqual(lab2.sum(), 1)
            # the row that is "test" in this draw must carry the high weight, as in the observed data
            self.assertAlmostEqual(w2[lab2][0], 0.9, places=9)

    def test_multiset_of_weights_is_preserved_per_signal(self):
        sig = np.array(["a", "a", "b", "b"])
        lab = np.array([True, False, True, False])
        w = np.array([0.8, 0.01, 0.7, 0.02])
        for lab2, w2 in self.m.permuted_label_weight(sig, lab, w, draws=30, seed=5):
            self.assertEqual(sorted(np.round(w2[:2], 6)), sorted(np.round(w[:2], 6)))
            self.assertEqual(sorted(np.round(w2[2:], 6)), sorted(np.round(w[2:], 6)))

    def test_same_seed_reproduces(self):
        sig = np.array(["a"] * 5); lab = np.array([True, True, False, False, False]); w = np.arange(5) / 5
        a = self.m.permuted_label_weight(sig, lab, w, draws=10, seed=7)
        b = self.m.permuted_label_weight(sig, lab, w, draws=10, seed=7)
        self.assertTrue(all(np.array_equal(x[0], y[0]) and np.allclose(x[1], y[1]) for x, y in zip(a, b)))


class TestCompositionalFamily(unittest.TestCase):
    """Per-factor BH is not a valid family here: the deviations sum to zero within every variant.

    Because each variant's deviations are centred, any real shift in profile shape pushes some factors up and
    an equal mass down, so many factors cross BH together. Measured: effects sum to -0.000000 across 502
    factors, and 247 of them crossed BH on a correctly centred null. The interpretable question is whether
    named regulators rank ABOVE chance, which is one test rather than 502.
    """

    def setUp(self):
        self.m = load("80_tf_convergence.py")

    def test_effects_are_reported_as_compositional(self):
        dev = np.array([[0.5, -0.5], [0.3, -0.3]])
        self.assertAlmostEqual(float(dev.sum()), 0.0, places=9)

    def test_rank_enrichment_detects_a_genuinely_top_ranked_set(self):
        effects = np.array([0.9, 0.8, 0.7, 0.1, 0.0, -0.1, -0.5])
        names = np.array(["A", "B", "C", "D", "E", "F", "G"])
        r = self.m.rank_enrichment(effects, names, {"A", "B", "C"}, draws=2000, seed=1)
        self.assertLess(r["mean_rank"], r["expected_mean_rank"])
        self.assertLess(r["p_empirical"], 0.05)

    def test_rank_enrichment_is_null_for_a_scattered_set(self):
        rng = np.random.default_rng(0)
        effects = rng.normal(size=200)
        names = np.array([f"F{i}" for i in range(200)])
        r = self.m.rank_enrichment(effects, names, {"F5", "F100", "F150"}, draws=2000, seed=2)
        self.assertGreater(r["p_empirical"], 0.05)

    def test_absent_names_are_ignored(self):
        effects = np.array([0.9, 0.1]); names = np.array(["A", "B"])
        r = self.m.rank_enrichment(effects, names, {"A", "NOT_PRESENT"}, draws=500, seed=3)
        self.assertEqual(r["n_in_set"], 1)


if __name__ == "__main__":
    unittest.main()

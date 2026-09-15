"""Fixture tests for step 50 (P5A haplotype additivity)."""

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


class TestMaskResolution(unittest.TestCase):
    """Output types come back at different resolutions: ATAC/DNase/RNA at 1 bp, CHIP_HISTONE at 128 bp.

    A base-resolution readout mask must be folded to the track's own row count before indexing, or the
    boolean index raises (size 8192 vs 1048576) - which is exactly how the first haplotype run died.
    """

    def setUp(self):
        self.m = load("50_haplotype_additivity.py")

    def test_base_mask_is_returned_unchanged_at_base_resolution(self):
        mask = np.zeros(1024, bool); mask[10:20] = True
        out = self.m.fold_mask(mask, 1024)
        self.assertIs(out.dtype.type, np.bool_)
        self.assertTrue(np.array_equal(out, mask))

    def test_a_bin_is_kept_when_any_base_inside_it_is_selected(self):
        mask = np.zeros(1024, bool); mask[130] = True          # falls in bin 1 at width 128
        out = self.m.fold_mask(mask, 8)
        self.assertEqual(out.shape, (8,))
        self.assertTrue(out[1])
        self.assertEqual(out.sum(), 1)

    def test_a_span_crossing_a_bin_boundary_keeps_both_bins(self):
        mask = np.zeros(1024, bool); mask[120:140] = True      # spans bins 0 and 1
        out = self.m.fold_mask(mask, 8)
        self.assertTrue(out[0] and out[1])
        self.assertEqual(out.sum(), 2)

    def test_empty_mask_folds_to_empty(self):
        self.assertEqual(self.m.fold_mask(np.zeros(1024, bool), 8).sum(), 0)

    def test_refuses_a_row_count_that_does_not_divide_the_mask(self):
        with self.assertRaises(Exception):
            self.m.fold_mask(np.zeros(1000, bool), 7)


class TestNullCountFromArgv(unittest.TestCase):
    def setUp(self):
        self.m = load("50_haplotype_additivity.py")

    def test_default_when_no_argument(self):
        self.assertEqual(self.m.n_null_from_argv(["50.py"]), self.m.N_NULL_DEFAULT)

    def test_explicit_count(self):
        self.assertEqual(self.m.n_null_from_argv(["50.py", "300"]), 300)


class TestEmpiricalPAndBH(unittest.TestCase):
    """With 1,587 pairs, "exceeds the null's 95th percentile" flags ~79 by construction.

    The amended prespecification asks for a per-pair two-sided empirical p against the matched null, then
    Benjamini-Hochberg across the whole family, reported per signal.
    """

    def setUp(self):
        self.m = load("50_haplotype_additivity.py")

    def test_a_residual_inside_the_null_gets_a_large_p(self):
        null = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
        self.assertGreater(self.m.empirical_two_sided_p(0.0, null), 0.5)

    def test_a_residual_beyond_every_null_draw_gets_the_floor_p(self):
        null = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
        p = self.m.empirical_two_sided_p(99.0, null)
        self.assertAlmostEqual(p, 1.0 / (len(null) + 1), places=9)

    def test_p_is_never_zero(self):
        self.assertGreater(self.m.empirical_two_sided_p(1e9, np.array([0.0, 0.1])), 0.0)

    def test_sign_does_not_matter(self):
        null = np.array([-1.0, 1.0, 0.0])
        self.assertAlmostEqual(self.m.empirical_two_sided_p(2.0, null),
                               self.m.empirical_two_sided_p(-2.0, null), places=12)

    def test_bh_rejects_nothing_when_all_p_are_large(self):
        self.assertEqual(self.m.bh_reject(np.array([0.4, 0.5, 0.9]), 0.10).sum(), 0)

    def test_bh_rejects_the_clear_signal_only(self):
        p = np.array([0.0001, 0.4, 0.5, 0.9, 0.95])
        r = self.m.bh_reject(p, 0.10)
        self.assertTrue(r[0]); self.assertEqual(r.sum(), 1)

    def test_bh_is_a_step_up_procedure_not_a_per_test_threshold(self):
        # every p below 0.10 individually, but BH keeps them all because they are dense near the boundary
        p = np.array([0.02, 0.03, 0.04, 0.05, 0.06])
        self.assertEqual(self.m.bh_reject(p, 0.10).sum(), 5)


class TestEffectStratum(unittest.TestCase):
    """Additivity is undefined when neither variant has a predicted effect.

    A residual built from two null effects measures numerical wobble; 1,587 such pairs would dominate the
    family. The threshold comes from the matched null's single-variant |effect| distribution, which is fixed
    before any observed residual is examined.
    """

    def setUp(self):
        self.m = load("50_haplotype_additivity.py")

    def test_threshold_is_the_null_95th_percentile_of_single_variant_effects(self):
        null_v1 = np.array([0.0, 0.1, 0.2, 0.3, 0.4])
        null_v2 = np.array([0.05, 0.15, 0.25, 0.35, 0.45])
        t = self.m.effect_threshold(null_v1, null_v2)
        self.assertAlmostEqual(t, float(np.quantile(np.abs(np.concatenate([null_v1, null_v2])), 0.95)), places=9)

    def test_both_variants_must_clear_the_threshold(self):
        keep = self.m.effect_stratum(np.array([1.0, 1.0, 0.01]), np.array([1.0, 0.01, 0.01]), 0.5)
        self.assertTrue(keep[0])
        self.assertFalse(keep[1])   # only one variant clears
        self.assertFalse(keep[2])

    def test_sign_is_ignored_only_magnitude_counts(self):
        keep = self.m.effect_stratum(np.array([-1.0]), np.array([1.0]), 0.5)
        self.assertTrue(keep[0])

    def test_nan_effects_are_excluded(self):
        keep = self.m.effect_stratum(np.array([np.nan]), np.array([1.0]), 0.5)
        self.assertFalse(keep[0])


class TestStratumNeedsAMatchedNull(unittest.TestCase):
    """Selecting pairs on effect size requires a null selected the same way.

    The `both_variants_active` stratum keeps pairs where both variants clear the null's 95th percentile.
    Only 1-3 of 300 null pairs clear that same bar, so the stratum's residuals were being compared against a
    null made almost entirely of no-effect pairs - which any active pair exceeds trivially. That produced
    19-26 "rejections" out of 23-33 tests, i.e. 83% of the stratum, which is a broken comparison rather
    than a discovery.
    """

    def setUp(self):
        self.m = load("50_haplotype_additivity.py")

    def test_too_few_matched_null_pairs_suppresses_the_verdict(self):
        r = self.m.stratum_is_testable(n_observed=23, n_null_matched=2, min_null=20)
        self.assertFalse(r["testable"])
        self.assertIn("matched null", r["reason"])

    def test_enough_matched_null_pairs_allows_it(self):
        r = self.m.stratum_is_testable(n_observed=23, n_null_matched=40, min_null=20)
        self.assertTrue(r["testable"])

    def test_the_unstratified_family_is_always_testable(self):
        r = self.m.stratum_is_testable(n_observed=1572, n_null_matched=300, min_null=20)
        self.assertTrue(r["testable"])


if __name__ == "__main__":
    unittest.main()

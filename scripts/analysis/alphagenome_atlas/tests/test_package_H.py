"""Fixture tests for Package H helpers (per-fold per-region skill, block bootstrap, partial R2, region keys)."""

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


class TestRegionSkill(unittest.TestCase):
    def setUp(self):
        self.m = load("18_package_H_region_covariates.py")

    def test_perfect_prediction_gives_skill_one_and_training_mean_gives_zero(self):
        rng = np.random.default_rng(0)
        fold = np.array([0, 0, 1, 1, 2, 2])
        target = rng.normal(size=(6, 3))
        # the producer's baseline is the TRAINING-fold mean of the fold's own residualised matrix, supplied per fold
        tm = {f: target[fold != f].mean(0) for f in (0, 1, 2)}
        skill = self.m.per_fold_region_skill(target, target.copy(), fold, tm)
        self.assertEqual(skill.shape, (3, 3))
        np.testing.assert_allclose(skill, 1.0)
        pred_mean = np.vstack([tm[f][None, :].repeat(2, 0) for f in (0, 1, 2)])
        np.testing.assert_allclose(self.m.per_fold_region_skill(target, pred_mean, fold, tm), 0.0, atol=1e-12)

    def test_zero_denominator_region_is_nan_not_inf(self):
        fold = np.array([0, 0, 1, 1])
        target = np.array([[1.0, 5.0], [1.0, 6.0], [1.0, 7.0], [1.0, 8.0]])   # region 0 is constant
        tm = {f: target[fold != f].mean(0) for f in (0, 1)}
        skill = self.m.per_fold_region_skill(target, target + 0.1, fold, tm)
        self.assertTrue(np.isnan(skill[:, 0]).all())
        self.assertTrue(np.isfinite(skill[:, 1]).all())


class TestBlockBootstrap(unittest.TestCase):
    def test_point_estimate_is_spearman_and_ci_is_deterministic_under_seed(self):
        m = load("18_package_H_region_covariates.py")
        rng = np.random.default_rng(1)
        x = rng.normal(size=400); y = x + rng.normal(size=400)
        blocks = np.repeat(np.arange(40), 10)
        a = m.block_bootstrap_spearman(x, y, blocks, draws=200, seed=20260909)
        b = m.block_bootstrap_spearman(x, y, blocks, draws=200, seed=20260909)
        from scipy.stats import spearmanr
        self.assertAlmostEqual(a["rho"], spearmanr(x, y).correlation, places=12)
        self.assertEqual((a["ci_low"], a["ci_high"]), (b["ci_low"], b["ci_high"]))
        self.assertLess(a["ci_low"], a["rho"]); self.assertGreater(a["ci_high"], a["rho"])
        self.assertEqual(a["n_blocks"], 40)

    def test_resamples_blocks_not_rows(self):
        m = load("18_package_H_region_covariates.py")
        # within-block x is constant and y = x, so any block resample keeps rho = 1 exactly
        x = np.repeat(np.arange(20, dtype=float), 5); y = x.copy()
        blocks = np.repeat(np.arange(20), 5)
        a = m.block_bootstrap_spearman(x, y, blocks, draws=50, seed=3)
        self.assertAlmostEqual(a["ci_low"], 1.0, places=12); self.assertAlmostEqual(a["ci_high"], 1.0, places=12)


class TestPartialR2(unittest.TestCase):
    def test_exact_linear_dependence_isolates_the_responsible_column(self):
        m = load("18_package_H_region_covariates.py")
        rng = np.random.default_rng(2)
        X = rng.normal(size=(500, 3)); y = 2.0 * X[:, 1]
        out = m.covariate_model(X, y, ["a", "b", "c"])
        self.assertAlmostEqual(out["r2"], 1.0, places=10)
        self.assertAlmostEqual(out["partial_r2"]["b"], 1.0, places=10)
        self.assertAlmostEqual(out["partial_r2"]["a"], 0.0, places=10)
        self.assertAlmostEqual(out["partial_r2"]["c"], 0.0, places=10)

    def test_nan_rows_are_dropped_and_counted(self):
        m = load("18_package_H_region_covariates.py")
        X = np.array([[1.0, 2.0], [2.0, 1.0], [np.nan, 0.0], [3.0, 4.0], [4.0, 3.0]]); y = np.array([1.0, 2.0, 3.0, 3.0, 4.0])
        out = m.covariate_model(X, y, ["a", "b"])
        self.assertEqual(out["n"], 4)


class TestRegionKeys(unittest.TestCase):
    def test_parse_and_block(self):
        m = load("18_package_H_region_covariates.py")
        self.assertEqual(m.parse_region_key("chr1:96466-96736"), ("chr1", 96466, 96736))
        self.assertEqual(m.mb_block("chr1", 96466, 96736), "chr1:0")
        self.assertEqual(m.mb_block("chr2", 1_999_990, 2_000_010), "chr2:2")   # midpoint decides
        with self.assertRaises(ValueError):
            m.parse_region_key("1:5-6")

    def test_cis_count_bins_are_ordered_and_exhaustive(self):
        m = load("18_package_H_region_covariates.py")
        self.assertEqual([m.cis_bin(k) for k in (0, 1, 2, 3, 4, 5, 9, 10, 50)], ["0", "1", "2", "3-4", "3-4", "5-9", "5-9", "10+", "10+"])


if __name__ == "__main__":
    unittest.main()

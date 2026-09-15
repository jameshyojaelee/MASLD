"""Unit tests for the four-arm aspect-axis instrument.

Every test is written so that it FAILS if the thing it checks stops happening:
the ceiling is checked against a hand-computed value and against a tamper, the
partial association is checked against scipy on residuals, the join guard is
handed a zero join on purpose, and the permutation null is calibrated on data
with a known answer in BOTH directions.
"""
from __future__ import annotations

import pathlib
import sys
import tempfile
import unittest

import numpy as np
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "scripts"))

import aspect_axis_common as C  # noqa: E402


class TestTieCeiling(unittest.TestCase):
    def test_reproduces_the_two_published_constants(self):
        self.assertAlmostEqual(C.tie_ceiling_from_counts([71, 15, 9, 4]),
                               0.791771, places=6)
        self.assertAlmostEqual(C.tie_ceiling_from_counts([41, 16, 1]),
                               0.791170, places=6)

    def test_hand_computed_value(self):
        # n=4, blocks 2/2: midranks 1.5,1.5,3.5,3.5.
        # sd(midranks)=1, sd(1..4)=sqrt(1.25) -> 0.894427
        self.assertAlmostEqual(C.tie_ceiling([0, 0, 1, 1]), 1 / np.sqrt(1.25),
                               places=9)

    def test_untied_label_has_ceiling_one(self):
        self.assertAlmostEqual(C.tie_ceiling(np.arange(50)), 1.0, places=12)

    def test_a_one_participant_tamper_moves_the_value(self):
        base = C.tie_ceiling_from_counts([71, 15, 9, 4])
        moved = C.tie_ceiling_from_counts([72, 15, 9, 4])
        self.assertGreater(abs(base - moved), 1e-4)

    def test_a_constant_label_has_a_zero_ceiling(self):
        # a constant label can reach no association at all
        self.assertEqual(C.tie_ceiling(np.zeros(30)), 0.0)


class TestPower(unittest.TestCase):
    def test_mde_shrinks_with_n_and_grows_with_alpha_strictness(self):
        self.assertLess(C.mde_spearman(200, 2, 0.05),
                        C.mde_spearman(50, 2, 0.05))
        self.assertGreater(C.mde_spearman(58, 2, 0.05 / 4019),
                           C.mde_spearman(58, 2, 0.05))

    def test_mde_is_undefined_when_df_is_exhausted(self):
        self.assertTrue(np.isnan(C.mde_spearman(5, 4, 0.05)))


class TestAssociation(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(11)
        self.n = 60
        self.z = rng.normal(size=self.n)
        self.x = self.z + rng.normal(size=self.n)
        self.Y = np.column_stack([self.z + rng.normal(size=self.n) for _ in range(7)])

    def test_partial_matches_scipy_on_residuals(self):
        Cd = np.column_stack([np.ones(self.n), rankdata(self.z)])
        Ry = np.apply_along_axis(rankdata, 0, self.Y)
        rho, a, frac, ok = C.partial_association(Ry, rankdata(self.x), Cd)
        ry = C.residualise(Ry, Cd)
        rx = C.residualise(rankdata(self.x).reshape(-1, 1), Cd).ravel()
        for j in range(self.Y.shape[1]):
            expected = np.corrcoef(ry[:, j], rx)[0, 1]
            self.assertAlmostEqual(rho[j], expected, places=10)

    def test_marginal_partial_equals_spearman_when_only_the_intercept_adjusts(self):
        Cd = np.ones((self.n, 1))
        Ry = np.apply_along_axis(rankdata, 0, self.Y)
        rho, _, _, _ = C.partial_association(Ry, rankdata(self.x), Cd)
        for j in range(self.Y.shape[1]):
            self.assertAlmostEqual(rho[j], spearmanr(self.x, self.Y[:, j]).statistic,
                                   places=10)

    def test_residual_variance_fraction_is_one_without_adjustment(self):
        Cd = np.ones((self.n, 1))
        Ry = np.apply_along_axis(rankdata, 0, self.Y)
        _, _, frac, _ = C.partial_association(Ry, rankdata(self.x), Cd)
        self.assertAlmostEqual(frac, 1.0, places=10)

    def test_an_exposure_determined_by_the_design_returns_no_test(self):
        x = rankdata(self.z)
        Cd = np.column_stack([np.ones(self.n), x])
        Ry = np.apply_along_axis(rankdata, 0, self.Y)
        rho, a, frac, ok = C.partial_association(Ry, x, Cd)
        self.assertIsNone(rho)
        self.assertLess(frac, 1e-12)


class TestBH(unittest.TestCase):
    def test_bh_matches_a_direct_implementation(self):
        rng = np.random.default_rng(3)
        p = rng.uniform(size=500)
        q = C.bh(p)
        n = p.size
        order = np.argsort(p)
        manual = np.minimum.accumulate(
            (p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
        out = np.empty(n)
        out[order] = manual
        np.testing.assert_allclose(q, np.clip(out, 0, 1))

    def test_all_ones_yield_no_discovery(self):
        self.assertEqual(int((C.bh(np.ones(100)) < 0.05).sum()), 0)


class TestDirectionCalibration(unittest.TestCase):
    """The null must be calibrated in BOTH directions: a planted signal must be
    found, and pure noise must not be."""

    def _run(self, planted, seed):
        rng = np.random.default_rng(seed)
        n, p = 80, 400
        z = rng.normal(size=n)
        x = z + rng.normal(size=n)
        Y = rng.normal(size=(n, p)) + np.outer(z, rng.normal(size=p))
        if planted:
            Y[:, :60] += np.outer(x, rng.normal(loc=1.2, scale=0.1, size=60))
        Cd = np.column_stack([np.ones(n), rankdata(z)])
        Ry = np.apply_along_axis(rankdata, 0, Y)
        return C.run_direction(Ry, rankdata(x), Cd, exposure="x",
                               adjusted_for=["z"],
                               rng=np.random.default_rng(seed + 1), n_perm=300)

    def test_a_planted_axis_is_recovered(self):
        r = self._run(True, 20260901)
        self.assertEqual(r.verdict, "UNIQUE")
        self.assertGreater(r.count_bh05, r.null_count_p95)
        self.assertLess(r.perm_p, 0.05)

    def test_pure_noise_is_not_called_unique(self):
        r = self._run(False, 771)
        self.assertNotEqual(r.verdict, "UNIQUE")
        self.assertGreaterEqual(r.perm_p, 0.05)

    def test_the_null_is_centred_near_the_nominal_false_discovery_rate(self):
        r = self._run(False, 99)
        self.assertLess(r.null_count_mean, 5.0)

    def test_a_collinear_direction_returns_no_count(self):
        rng = np.random.default_rng(5)
        n, p = 40, 100
        x = rng.normal(size=n)
        Y = rng.normal(size=(n, p))
        Cd = np.column_stack([np.ones(n), rankdata(x)])
        Ry = np.apply_along_axis(rankdata, 0, Y)
        r = C.run_direction(Ry, rankdata(x), Cd, exposure="x",
                            adjusted_for=["x_copy"],
                            rng=np.random.default_rng(6), n_perm=50)
        self.assertEqual(r.verdict, "NOT_APPLICABLE_COLLINEAR")
        self.assertEqual(r.count_bh05, 0)


class TestJoinGuards(unittest.TestCase):
    def test_a_zero_join_raises(self):
        with self.assertRaises(C.JoinError):
            C.require_join(0, 58, "deliberate zero join")

    def test_a_wrong_size_join_raises(self):
        with self.assertRaises(C.JoinError):
            C.require_join(57, 58, "deliberate short join")

    def test_a_correct_join_passes(self):
        self.assertEqual(C.require_join(58, 58, "ok"), 58)

    def test_a_ragged_pandas_index_header_raises(self):
        with tempfile.TemporaryDirectory() as d:
            p = pathlib.Path(d) / "ragged.tsv"
            p.write_text("a\tb\n1\t2\t3\n4\t5\t6\n")
            with self.assertRaises(C.JoinError):
                C.read_tsv_strict(p, label="ragged fixture")

    def test_a_well_formed_table_is_marked_not_ragged(self):
        with tempfile.TemporaryDirectory() as d:
            p = pathlib.Path(d) / "clean.tsv"
            p.write_text("a\tb\n1\t2\n3\t4\n")
            frame = C.read_tsv_strict(p, label="clean fixture")
            self.assertFalse(frame.attrs["ragged_header"])
            self.assertEqual(list(frame.columns), ["a", "b"])


class TestMinimumDiscoveryCount(unittest.TestCase):
    """A verdict must not rest on one feature. This plants a signal in exactly
    THREE features out of 800 so the count clears a null that is a point mass at
    zero, and asserts the criterion refuses it."""

    def _run(self, n_planted, seed=4242):
        rng = np.random.default_rng(seed)
        n, p = 70, 800
        z = rng.normal(size=n)
        x = z + rng.normal(size=n)
        Y = rng.normal(size=(n, p))
        if n_planted:
            Y[:, :n_planted] += np.outer(x, np.full(n_planted, 3.0))
        Cd = np.column_stack([np.ones(n), rankdata(z)])
        Ry = np.apply_along_axis(rankdata, 0, Y)
        return C.run_direction(Ry, rankdata(x), Cd, exposure="x",
                               adjusted_for=["z"],
                               rng=np.random.default_rng(seed + 1), n_perm=300)

    def test_three_features_clear_the_null_but_are_refused(self):
        r = self._run(3)
        self.assertLess(r.perm_p, 0.05)
        self.assertGreater(r.count_bh05, r.null_count_p95)
        self.assertLess(r.count_bh05, C.MIN_DISCOVERY_COUNT)
        self.assertEqual(r.verdict, "UNIQUE_BUT_BELOW_MINIMUM_COUNT")
        self.assertTrue(r.notes)

    def test_the_same_construction_with_forty_features_is_accepted(self):
        r = self._run(40)
        self.assertGreaterEqual(r.count_bh05, C.MIN_DISCOVERY_COUNT)
        self.assertEqual(r.verdict, "UNIQUE")


class TestStouffer(unittest.TestCase):
    def test_two_identical_arms_beat_either_alone(self):
        z, p, shares = C.stouffer([0.01, 0.01], [1.0, 1.0])
        self.assertLess(p, 0.01)
        self.assertAlmostEqual(sum(shares), z, places=10)

    def test_a_permutation_p_of_exactly_one_does_not_destroy_the_combination(self):
        z, p, shares = C.stouffer([1.0, 1.0], [1.0, 1.0], n_perm=5000)
        self.assertTrue(np.isfinite(z))
        self.assertTrue(np.isfinite(p))

    def test_the_clamp_uses_the_permutation_resolution(self):
        z_small, _, _ = C.stouffer([0.0, 0.0], [1.0, 1.0], n_perm=100)
        z_large, _, _ = C.stouffer([0.0, 0.0], [1.0, 1.0], n_perm=100000)
        self.assertLess(z_small, z_large)

    def test_the_share_vector_exposes_a_single_arm_driver(self):
        z, p, shares = C.stouffer([1e-4, 0.5], [1.0, 1.0])
        self.assertGreater(shares[0], 0.9 * z)


if __name__ == "__main__":
    unittest.main()

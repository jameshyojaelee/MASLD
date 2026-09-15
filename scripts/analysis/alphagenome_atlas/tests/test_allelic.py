"""Fixture tests for step 73 (same-variant snATAC allelic concordance)."""

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


class TestParseAD(unittest.TestCase):
    """bcftools writes AD over [REF] + ALT alleles, where ALT always ends in the <*> catch-all."""

    def setUp(self):
        self.m = load("73_atac_allelic_analysis.py")

    def test_target_alt_present(self):
        self.assertEqual(self.m.parse_ad("A", "G,<*>", "165,51,0", "G"), (165, 51, 0))

    def test_target_alt_absent_gives_zero_alt_not_a_miscount(self):
        # pileup saw only the reference; the listed ALT is the <*> catch-all, not our variant
        self.assertEqual(self.m.parse_ad("A", "<*>", "185,0", "G"), (185, 0, 0))

    def test_third_allele_goes_to_other_not_to_alt(self):
        # ref C, our variant is C>G, but the pileup also saw A: A's reads must not be counted as alt
        self.assertEqual(self.m.parse_ad("C", "A,G,<*>", "0,124,1,0", "G"), (0, 1, 124))

    def test_star_catchall_reads_are_other(self):
        self.assertEqual(self.m.parse_ad("A", "G,<*>", "10,5,3", "G"), (10, 5, 3))


class TestIndelRecords(unittest.TestCase):
    """bcftools writes a SECOND record at a position where it also sees an indel; it is not our SNV."""

    def setUp(self):
        self.m = load("73_atac_allelic_analysis.py")

    def test_multibase_ref_record_is_skipped(self):
        self.assertTrue(self.m.is_indel_record("CAAAAAAA", "CAAAAA,CAAAAAAAA"))

    def test_snv_record_at_the_same_position_is_kept(self):
        self.assertFalse(self.m.is_indel_record("C", "<*>"))
        self.assertFalse(self.m.is_indel_record("C", "T,<*>"))

    def test_multibase_alt_is_an_indel_even_with_single_base_ref(self):
        self.assertTrue(self.m.is_indel_record("C", "CTT,<*>"))


class TestHetCall(unittest.TestCase):
    def setUp(self):
        self.m = load("73_atac_allelic_analysis.py")

    def test_balanced_deep_site_is_het(self):
        self.assertTrue(self.m.call_het(20, 18, 0))

    def test_shallow_site_is_not_het(self):
        self.assertFalse(self.m.call_het(5, 4, 0))

    def test_homozygous_reference_is_not_het(self):
        self.assertFalse(self.m.call_het(40, 0, 0))

    def test_one_stray_alt_read_is_not_het(self):
        self.assertFalse(self.m.call_het(40, 1, 0))

    def test_contaminated_site_with_many_other_reads_is_refused(self):
        self.assertFalse(self.m.call_het(20, 18, 30))


class TestSiteStatistic(unittest.TestCase):
    def setUp(self):
        self.m = load("73_atac_allelic_analysis.py")

    def test_donor_is_the_unit_not_the_read(self):
        # one deep donor must not outvote three shallow donors pointing the other way
        deep = [(2, 200)]                       # log2 ~ +6.3
        shallow = [(20, 10), (20, 10), (20, 10)]  # log2 ~ -1.0 each
        stat = self.m.site_statistic(deep + shallow)
        self.assertLess(stat["mean_log2"], 1.0)
        self.assertEqual(stat["n_donors"], 4)

    def test_symmetric_counts_give_zero(self):
        self.assertAlmostEqual(self.m.site_statistic([(10, 10), (30, 30)])["mean_log2"], 0.0, places=9)


class TestStratifiedPermutation(unittest.TestCase):
    def setUp(self):
        self.m = load("73_atac_allelic_analysis.py")

    def test_permutation_preserves_predictor_magnitude_within_stratum(self):
        # the null must reassign WHICH site gets a prediction, keeping |quantile| structure
        pred = np.array([0.9, -0.8, 0.1, -0.2])
        meas = np.array([1.0, -1.0, 1.0, -1.0])
        strata = np.array([0, 0, 1, 1])
        p, draws = self.m.stratified_permutation_concordance(pred, meas, strata, draws=200, seed=7)
        self.assertEqual(len(draws), 200)
        self.assertTrue(all(0.0 <= d <= 1.0 for d in draws))
        # within each stratum only 2 orderings exist, so the null takes exactly two values
        self.assertLessEqual(len(set(np.round(draws, 6))), 3)

    def test_identical_seed_reproduces(self):
        pred = np.array([0.9, -0.8, 0.1, -0.2, 0.5, -0.4])
        meas = np.array([1.0, -1.0, 1.0, -1.0, 1.0, 1.0])
        strata = np.array([0, 0, 1, 1, 2, 2])
        a = self.m.stratified_permutation_concordance(pred, meas, strata, draws=50, seed=11)[1]
        b = self.m.stratified_permutation_concordance(pred, meas, strata, draws=50, seed=11)[1]
        self.assertEqual(list(a), list(b))


class TestMarginalExpectation(unittest.TestCase):
    """A shared skew (e.g. both sides favouring the reference allele) manufactures concordance on its own."""

    def setUp(self):
        self.m = load("73_atac_allelic_analysis.py")

    def test_two_all_negative_vectors_agree_perfectly_but_expectation_is_also_one(self):
        pred = np.array([-1.0, -2.0, -3.0])
        meas = np.array([-1.0, -1.0, -1.0])
        e = self.m.marginal_expected_concordance(pred, meas)
        self.assertAlmostEqual(e, 1.0, places=9)

    def test_balanced_signs_give_one_half(self):
        pred = np.array([1.0, -1.0, 1.0, -1.0])
        meas = np.array([1.0, 1.0, -1.0, -1.0])
        self.assertAlmostEqual(self.m.marginal_expected_concordance(pred, meas), 0.5, places=9)


class TestBlockLevelConcordance(unittest.TestCase):
    """Sites inside one 1-Mb block are in LD; the independent unit is the block, not the site."""

    def setUp(self):
        self.m = load("73_atac_allelic_analysis.py")

    def test_one_large_block_cannot_outvote_many_small_ones(self):
        # block A: 100 concordant sites. blocks B..D: 2 discordant sites each.
        pred = np.array([1.0] * 100 + [1.0] * 6)
        meas = np.array([1.0] * 100 + [-1.0] * 6)
        blocks = np.array(["A"] * 100 + ["B", "B", "C", "C", "D", "D"])
        r = self.m.block_level_concordance(pred, meas, blocks)
        self.assertEqual(r["n_blocks"], 4)
        self.assertAlmostEqual(r["mean_per_block"], 0.25, places=9)      # 1 of 4 blocks concordant
        self.assertEqual(r["n_blocks_above_half"], 1)
        self.assertAlmostEqual(r["site_level"], 100 / 106, places=9)     # what the site-level number would say

    def test_sign_test_p_is_exact_binomial_on_blocks(self):
        pred = np.array([1.0, 1.0, 1.0, 1.0])
        meas = np.array([1.0, 1.0, 1.0, -1.0])
        blocks = np.array(["A", "B", "C", "D"])
        r = self.m.block_level_concordance(pred, meas, blocks)
        self.assertEqual(r["n_blocks_above_half"], 3)
        self.assertAlmostEqual(r["sign_test_p"], 5 / 16, places=9)       # P(X >= 3 | n=4, p=.5) = (4 + 1)/16


class TestCohortResolution(unittest.TestCase):
    def setUp(self):
        self.m = load("73_atac_allelic_analysis.py")

    def test_default_is_the_discovery_cohort(self):
        self.assertEqual(self.m.resolve_cohort(["73.py"])[0], "gse281367")

    def test_replication_cohort_selects_its_own_counts_and_prefix(self):
        name, counts_dir, prefix = self.m.resolve_cohort(["73.py", "gse244832"])
        self.assertEqual((name, counts_dir, prefix), ("gse244832", "atac_counts_gse244832", "gse244832_"))

    def test_unknown_cohort_is_refused(self):
        import lib_atlas as la
        with self.assertRaises(la.ContractError):
            self.m.resolve_cohort(["73.py", "gse999999"])


if __name__ == "__main__":
    unittest.main()

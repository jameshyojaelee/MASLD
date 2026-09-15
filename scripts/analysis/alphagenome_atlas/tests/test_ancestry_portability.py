"""Fixture tests for step 62 (P6b: is the Atlas annotation ancestry-portable?)."""

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


class TestMassDecomposition(unittest.TestCase):
    def setUp(self):
        self.m = load("62_ancestry_portability.py")

    def test_shares_sum_to_one(self):
        d = self.m.mass_decomposition({"total_mass": "1.0", "queried_mass": "0.8",
                                       "excluded_mass_below_floor": "0.15", "unmapped_mass": "0.05"})
        self.assertAlmostEqual(d["queried_share"] + d["floor_share"] + d["excluded_share"], 1.0, places=12)
        self.assertAlmostEqual(d["queried_share"], 0.8, places=12)

    def test_components_that_do_not_reconstruct_the_total_are_refused(self):
        """The identity is a guard, not a comment: a row whose parts miss the total is a defect upstream."""
        with self.assertRaises(Exception):
            self.m.mass_decomposition({"total_mass": "1.0", "queried_mass": "0.8",
                                       "excluded_mass_below_floor": "0.05", "unmapped_mass": "0.0"})

    def test_zero_total_mass_is_refused_rather_than_dividing(self):
        with self.assertRaises(Exception):
            self.m.mass_decomposition({"total_mass": "0.0", "queried_mass": "0.0",
                                       "excluded_mass_below_floor": "0.0", "unmapped_mass": "0.0"})


class TestChannelAvailability(unittest.TestCase):
    def setUp(self):
        self.m = load("62_ancestry_portability.py")

    def test_share_is_weighted_by_posterior_not_by_variant_count(self):
        """One heavy served variant and nine light unserved ones is a high share, not 0.1."""
        weights = {"v0": 0.9}
        weights.update({f"v{i}": 0.1 / 9 for i in range(1, 10)})
        served = {"v0": True}
        served.update({f"v{i}": False for i in range(1, 10)})
        self.assertAlmostEqual(self.m.channel_available_share(weights, served), 0.9, places=12)

    def test_missing_variant_counts_as_unserved(self):
        self.assertAlmostEqual(self.m.channel_available_share({"a": 0.5, "b": 0.5}, {"a": True}), 0.5, places=12)


class TestBlockBootstrap(unittest.TestCase):
    def setUp(self):
        self.m = load("62_ancestry_portability.py")

    def test_one_signal_per_block_per_draw(self):
        """Ten signals in one block must not act as ten independent units."""
        rows = [{"analysis_block": "b1", "value": v} for v in [0.1] * 9 + [0.9]]
        drawn = self.m.block_draw(rows, seed=1)
        self.assertEqual(len(drawn), 1)

    def test_interval_of_a_constant_difference_is_that_difference(self):
        a = [{"analysis_block": f"b{i}", "value": 0.5} for i in range(20)]
        b = [{"analysis_block": f"c{i}", "value": 0.3} for i in range(20)]
        res = self.m.block_bootstrap_median_difference(a, b, draws=200, seed=7)
        self.assertAlmostEqual(res["observed"], 0.2, places=12)
        self.assertAlmostEqual(res["lo"], 0.2, places=12)
        self.assertAlmostEqual(res["hi"], 0.2, places=12)

    def test_interval_widens_when_the_blocks_disagree(self):
        a = [{"analysis_block": f"b{i}", "value": 0.1 * (i % 10)} for i in range(20)]
        b = [{"analysis_block": f"c{i}", "value": 0.5} for i in range(20)]
        res = self.m.block_bootstrap_median_difference(a, b, draws=400, seed=7)
        self.assertGreater(res["hi"] - res["lo"], 0.0)


class TestMatching(unittest.TestCase):
    def setUp(self):
        self.m = load("62_ancestry_portability.py")

    def test_comparators_come_from_the_same_size_decile(self):
        pool = [{"signal_uid": f"e{i}", "effective_set_size": n, "analysis_block": f"b{i}"}
                for i, n in enumerate([1, 2, 3, 100, 200, 300])]
        target = {"signal_uid": "t", "effective_set_size": 2, "analysis_block": "z"}
        edges = self.m.decile_edges([r["effective_set_size"] for r in pool])
        got = self.m.matched_comparators(target, pool, edges)
        self.assertTrue(all(r["effective_set_size"] < 50 for r in got), got)
        self.assertTrue(got)

    def test_a_target_with_no_comparator_in_its_decile_returns_empty(self):
        pool = [{"signal_uid": "e1", "effective_set_size": 1, "analysis_block": "b1"}]
        edges = self.m.decile_edges([1])
        got = self.m.matched_comparators({"signal_uid": "t", "effective_set_size": 9999, "analysis_block": "z"}, pool, edges)
        self.assertEqual(got, [])


class TestConclusiveness(unittest.TestCase):
    def setUp(self):
        self.m = load("62_ancestry_portability.py")

    def test_too_few_blocks_cannot_conclude(self):
        v = self.m.contrast_verdict(n_signals=40, n_blocks=4)
        self.assertFalse(v["conclusive"])
        self.assertIn("block", v["reason"])

    def test_too_few_signals_cannot_conclude(self):
        v = self.m.contrast_verdict(n_signals=6, n_blocks=6)
        self.assertFalse(v["conclusive"])

    def test_adequate_arm_is_conclusive(self):
        self.assertTrue(self.m.contrast_verdict(n_signals=98, n_blocks=56)["conclusive"])


class TestExcludedMassIsNotLiftover(unittest.TestCase):
    """`unmapped_mass` in signal_query_coverage is total minus MAPPED mass, and `mapped` excludes indels.

    The Atlas point-query API refuses indels (UNIMPLEMENTED), so the third term of the decomposition is
    dominated by deferred indels, not by liftover failure. Attributing it to coordinates would be wrong.
    """

    def setUp(self):
        self.m = load("62_ancestry_portability.py")

    def test_mass_is_split_by_the_reason_the_variant_was_excluded(self):
        pairs = [{"weight": "0.6", "mapping_status": "excluded", "exclusion_reason": "indel_deferred"},
                 {"weight": "0.1", "mapping_status": "excluded", "exclusion_reason": "liftover_failed"},
                 {"weight": "0.3", "mapping_status": "mapped", "exclusion_reason": ""}]
        got = self.m.excluded_mass_by_reason(pairs)
        self.assertAlmostEqual(got["indel_deferred"], 0.6, places=12)
        self.assertAlmostEqual(got["liftover_failed"], 0.1, places=12)
        self.assertNotIn("", got)

    def test_a_signal_with_no_exclusions_reports_no_reasons(self):
        self.assertEqual(self.m.excluded_mass_by_reason(
            [{"weight": "1.0", "mapping_status": "mapped", "exclusion_reason": ""}]), {})


class TestEffectiveSetSize(unittest.TestCase):
    """n_cs_variants is 0 for every coloc signal, so it cannot match on fine-mapping resolution."""

    def setUp(self):
        self.m = load("62_ancestry_portability.py")

    def test_a_point_mass_has_effective_size_one(self):
        self.assertAlmostEqual(self.m.effective_set_size({"a": 1.0, "b": 0.0, "c": 0.0}), 1.0, places=9)

    def test_uniform_mass_over_k_variants_has_effective_size_k(self):
        self.assertAlmostEqual(self.m.effective_set_size({f"v{i}": 0.25 for i in range(4)}), 4.0, places=9)

    def test_unnormalised_weights_give_the_same_answer(self):
        self.assertAlmostEqual(self.m.effective_set_size({f"v{i}": 7.5 for i in range(4)}), 4.0, places=9)

    def test_empty_or_zero_mass_is_refused(self):
        with self.assertRaises(Exception):
            self.m.effective_set_size({})

    def test_matching_key_is_not_degenerate_across_coloc_signals(self):
        """The defect this replaces: every coloc signal had n_cs_variants == 0, so every EUR signal
        matched every target and the matched contrast equalled the unmatched one."""
        a = self.m.effective_set_size({"x": 0.9, "y": 0.1})
        b = self.m.effective_set_size({f"v{i}": 1.0 for i in range(50)})
        self.assertGreater(b, a * 5)


class TestCoverageTail(unittest.TestCase):
    def setUp(self):
        self.m = load("62_ancestry_portability.py")

    def test_tail_share_counts_signals_below_the_construction_floor(self):
        vals = [1.0, 0.9995, 0.999, 0.8, 0.0]
        self.assertAlmostEqual(self.m.tail_share(vals, 0.999), 0.4, places=12)

    def test_a_ceiling_arm_has_zero_tail(self):
        self.assertAlmostEqual(self.m.tail_share([1.0, 0.9999, 0.999], 0.999), 0.0, places=12)


class TestCoverageState(unittest.TestCase):
    """A per-signal state downstream consumers can gate on.

    Thirteen signals carry their whole posterior on indels the Atlas point-query API will not serve, yet
    still receive a profile row whose TOP-VARIANT fields (avi_top_quantile and friends) are unweighted
    maxima over 2,000 near-zero-weight SNVs. The posterior-weighted columns correctly read ~1e-23, but a
    reader who takes a top-variant quantile as the signal's mechanism would be reading a number that rests
    on no posterior mass. The state names that.
    """

    def setUp(self):
        self.m = load("62_ancestry_portability.py")

    def test_full_coverage_is_covered(self):
        self.assertEqual(self.m.coverage_state(0.9995, 0.0005), "covered")

    def test_all_mass_on_unscorable_indels_is_not_covered(self):
        st = self.m.coverage_state(3.3e-23, 1.0)
        self.assertEqual(st, "uncovered_indel")
        self.assertNotEqual(st, "covered")

    def test_partial_indel_loss_is_flagged_separately_from_full_loss(self):
        self.assertEqual(self.m.coverage_state(0.6, 0.4), "partial_indel")

    def test_loss_that_is_not_indels_is_named_for_what_it_is(self):
        self.assertEqual(self.m.coverage_state(0.2, 0.0), "partial_other")


if __name__ == "__main__":
    unittest.main()

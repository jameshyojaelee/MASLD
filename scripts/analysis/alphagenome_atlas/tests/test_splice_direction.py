"""Fixture tests for step 36, the signed model-API splice statistics."""

from __future__ import annotations

import importlib.util
import os
import pathlib
import sys
import tempfile
import types
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


class Junction:
    def __init__(self, start, end):
        self.start, self.end = start, end


class TestCoordinates(unittest.TestCase):
    def setUp(self):
        self.m = load("36_splice_direction.py")

    def test_an_object_with_start_and_end_is_read_directly(self):
        self.assertEqual(self.m.junction_coords(Junction(10, 20)), (10, 20))

    def test_the_string_form_is_the_fallback(self):
        self.assertEqual(self.m.junction_coords(":43464-55401:+"), (43464, 55401))

    def test_an_unparseable_junction_is_refused(self):
        with self.assertRaises(self.m.la.ContractError):
            self.m.junction_coords("nonsense")


class TestExcisionRatio(unittest.TestCase):
    def setUp(self):
        self.m = load("36_splice_direction.py")
        self.t = (100, 200)
        self.s = (100, 300)

    def test_lifting_the_whole_cluster_equally_changes_no_ratio(self):
        ref = {self.t: np.array([1.0, 2.0]), self.s: np.array([3.0, 4.0])}
        alt = {k: v * 5.0 for k, v in ref.items()}
        self.assertAlmostEqual(self.m.excision_ratio_log2(self.t, [self.s], ref, alt), 0.0, places=6)

    def test_lifting_only_the_target_is_positive_and_hand_checkable(self):
        ref = {self.t: np.array([1.0]), self.s: np.array([1.0])}
        alt = {self.t: np.array([3.0]), self.s: np.array([1.0])}
        # ref share 1/2, alt share 3/4 -> log2(0.75/0.5) = log2(1.5)
        self.assertAlmostEqual(self.m.excision_ratio_log2(self.t, [self.s], ref, alt),
                               float(np.log2(1.5)), places=5)

    def test_losing_the_target_is_negative(self):
        ref = {self.t: np.array([2.0]), self.s: np.array([1.0])}
        alt = {self.t: np.array([0.5]), self.s: np.array([1.0])}
        self.assertLess(self.m.excision_ratio_log2(self.t, [self.s], ref, alt), 0.0)

    def test_a_junction_with_no_cluster_member_has_no_ratio(self):
        ref = {self.t: np.array([1.0])}
        self.assertTrue(np.isnan(self.m.excision_ratio_log2(self.t, [], ref, ref)))

    def test_the_raw_usage_ratio_ignores_the_cluster(self):
        ref = {self.t: np.array([1.0, 1.0])}
        alt = {self.t: np.array([2.0, 4.0])}
        self.assertAlmostEqual(self.m.usage_log2(self.t, ref, alt), (1.0 + 2.0) / 2, places=5)


class TestCluster(unittest.TestCase):
    def setUp(self):
        self.m = load("36_splice_direction.py")

    def test_a_shared_donor_or_acceptor_makes_a_sibling(self):
        keys = [(100, 200), (100, 300), (50, 200), (500, 600)]
        index = self.m.cluster_index(keys)
        self.assertEqual(sorted(self.m.cluster_of((100, 200), index)), [(50, 200), (100, 300)])

    def test_an_unrelated_junction_is_not_a_sibling(self):
        self.assertEqual(self.m.cluster_of((100, 200), self.m.cluster_index([(500, 600)])), [])

    def test_the_target_itself_is_never_its_own_sibling(self):
        keys = [(100, 200), (100, 300), (400, 200)]
        out = self.m.cluster_of((100, 200), self.m.cluster_index(keys))
        self.assertNotIn((100, 200), out)
        self.assertEqual(len(out), len(set(out)))


class TestSiteUsage(unittest.TestCase):
    def setUp(self):
        self.m = load("36_splice_direction.py")

    def out(self, values):
        return types.SimpleNamespace(splice_site_usage=types.SimpleNamespace(values=values))

    def test_the_pad_averages_the_site_and_its_neighbours(self):
        v = np.zeros((10, 1))
        v[4:7, 0] = [1.0, 2.0, 3.0]
        self.assertAlmostEqual(self.m.site_usage_at(self.out(v), [5], pad=1), 2.0)

    def test_both_ends_are_included(self):
        v = np.zeros((20, 1))
        v[5, 0], v[15, 0] = 3.0, 6.0
        self.assertAlmostEqual(self.m.site_usage_at(self.out(v), [5, 15], pad=0), 4.5)

    def test_a_missing_output_is_missing_not_zero(self):
        self.assertTrue(np.isnan(self.m.site_usage_at(types.SimpleNamespace(), [5])))


class TestOffsetMeasurement(unittest.TestCase):
    def setUp(self):
        self.m = load("36_splice_direction.py")

    def frame(self):
        return pd.DataFrame([{"variant_uid": "v1", "intron_start": 1000, "intron_end": 2000},
                             {"variant_uid": "v1", "intron_start": 3000, "intron_end": 4000}])

    def test_the_expected_convention_wins_uniquely_when_the_data_follow_it(self):
        keys = {"v1": {(1000 - 500, 2000 - 1 - 500), (3000 - 500, 4000 - 1 - 500)}}
        r = self.m.measure_offsets(self.frame(), keys, {"v1": 500})
        self.assertEqual(r["best"], "0,-1")
        self.assertEqual(r["n_tied_at_best"], 1)
        self.assertTrue(r["agrees_with_expected"])

    def test_a_different_convention_is_reported_and_refused(self):
        keys = {"v1": {(1000 - 500, 2000 + 1 - 500), (3000 - 500, 4000 + 1 - 500)}}
        r = self.m.measure_offsets(self.frame(), keys, {"v1": 500})
        self.assertEqual(r["best"], "0,1")
        self.assertFalse(r["agrees_with_expected"])

    def test_no_match_anywhere_leaves_a_tie_and_is_refused(self):
        r = self.m.measure_offsets(self.frame(), {"v1": {(7, 9)}}, {"v1": 500})
        self.assertGreater(r["n_tied_at_best"], 1)
        self.assertFalse(r["agrees_with_expected"])


class TestBlockConcordance(unittest.TestCase):
    def setUp(self):
        self.m = load("36_splice_direction.py")

    def test_perfect_agreement_and_its_marginal_expectation(self):
        pred = np.array([1.0, -1.0, 1.0, -1.0])
        meas = np.array([2.0, -2.0, 3.0, -4.0])
        r = self.m.block_concordance(pred, meas, np.array(["b1", "b2", "b3", "b4"]), nboot=100)
        self.assertEqual(r["concordance"], 1.0)
        self.assertAlmostEqual(r["marginal_expectation"], 0.5)
        self.assertEqual(r["blocks_above_half"], 4)

    def test_an_always_positive_predictor_has_no_concordance_above_its_marginal(self):
        pred = np.ones(10)
        meas = np.array([1.0] * 7 + [-1.0] * 3)
        r = self.m.block_concordance(pred, meas, np.array([f"b{i}" for i in range(10)]), nboot=100)
        self.assertAlmostEqual(r["concordance"], 0.7)
        self.assertAlmostEqual(r["marginal_expectation"], 0.7)
        self.assertAlmostEqual(r["above_marginal"], 0.0)

    def test_exact_zero_predictions_are_excluded(self):
        r = self.m.block_concordance(np.array([0.0, 1.0]), np.array([1.0, 1.0]),
                                     np.array(["b1", "b2"]), nboot=50)
        self.assertEqual(r["n"], 1)

    def test_an_empty_input_returns_no_concordance(self):
        r = self.m.block_concordance(np.array([np.nan]), np.array([1.0]), np.array(["b1"]), nboot=10)
        self.assertEqual(r["n"], 0)
        self.assertIsNone(r["concordance"])


class TestColumnGuard(unittest.TestCase):
    def setUp(self):
        self.m = load("36_splice_direction.py")

    def test_an_absent_column_becomes_all_missing_of_the_right_length(self):
        df = pd.DataFrame({"a": [1.0, 2.0]})
        v = self.m.column(df, "not_there")
        self.assertEqual(len(v), 2)
        self.assertTrue(np.all(np.isnan(v)))


class TestVerdicts(unittest.TestCase):
    def setUp(self):
        self.m = load("36_splice_direction.py")

    def base(self, conc, above, p, raw, det, atlas, site, blocks=200):
        return {"direction": {"excision_ratio_log2": {"concordance": conc, "above_marginal": above,
                                                      "block_binomial_p": p, "n_blocks": blocks},
                              "target_usage_log2": {"concordance": raw},
                              "site_usage_delta": {"concordance": site}},
                "detection": {"rank_abs_excision_ratio": {"mean": det},
                              "rank_abs_atlas_quantile": {"mean": atlas}}}

    def test_direction_recovered_triggers_d2(self):
        v = self.m.verdicts(self.base(0.72, 0.18, 1e-6, 0.66, 0.75, 0.49, 0.70))
        self.assertTrue(v["P1_direction_recovered"])
        self.assertTrue(v["decision"].startswith("D2"))
        self.assertTrue(v["P5_model_api_beats_atlas_detection"])

    def test_a_concordance_at_its_marginal_expectation_is_not_direction(self):
        v = self.m.verdicts(self.base(0.72, 0.01, 1e-6, 0.66, 0.75, 0.49, 0.70))
        self.assertFalse(v["P1_direction_recovered"])
        self.assertTrue(v["decision"].startswith("D1"))

    def test_a_nominally_high_concordance_with_a_weak_block_test_fails(self):
        v = self.m.verdicts(self.base(0.72, 0.18, 0.2, 0.66, 0.75, 0.49, 0.70))
        self.assertFalse(v["P1_direction_recovered"])

    def test_too_few_blocks_is_indeterminate(self):
        v = self.m.verdicts(self.base(0.72, 0.18, 1e-6, 0.66, 0.75, 0.49, 0.70, blocks=5))
        self.assertTrue(v["decision"].startswith("D3"))

    def test_the_d1_text_names_the_hsd17b13_consequence(self):
        v = self.m.verdicts(self.base(0.50, 0.00, 0.9, 0.50, 0.50, 0.49, 0.50))
        self.assertIn("HSD17B13", v["decision"])


if __name__ == "__main__":
    unittest.main()

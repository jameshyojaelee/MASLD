"""Fixture tests for step 43, the per-region saturation features."""

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


class TestPositions(unittest.TestCase):
    def setUp(self):
        self.m = load("43_saturation_features.py")

    def test_rows_map_to_sorted_positions(self):
        inv, pos = self.m.position_index(["chr1:12:A>C", "chr1:11:G>T", "chr1:12:A>G"], 10, 20)
        self.assertEqual(pos.tolist(), [11, 12])
        self.assertEqual(inv.tolist(), [1, 0, 1])

    def test_a_variant_outside_the_region_is_refused(self):
        with self.assertRaises(self.m.la.ContractError):
            self.m.position_index(["chr1:10:A>C"], 10, 20)      # start0 = 10 means the first base is 11

    def test_an_indel_string_is_refused(self):
        with self.assertRaises(self.m.la.ContractError):
            self.m.position_index(["chr1:12:A>AC"], 10, 20)

    def test_four_substitutions_at_one_position_is_refused(self):
        with self.assertRaises(self.m.la.ContractError):
            self.m.position_index(["chr1:12:A>C"] * 4, 10, 20)


class TestPerPositionMax(unittest.TestCase):
    def setUp(self):
        self.m = load("43_saturation_features.py")

    def test_the_maximum_over_substitutions_is_taken_per_track(self):
        q = np.array([[0.1, 0.9], [0.5, 0.2], [0.3, 0.4]])
        inv = np.array([0, 0, 1])
        out = self.m.per_position_max(q, inv, 2)
        self.assertEqual(out.tolist(), [[0.5, 0.9], [0.3, 0.4]])

    def test_a_missing_value_does_not_hide_a_present_one(self):
        q = np.array([[np.nan, 0.2], [0.7, np.nan]])
        out = self.m.per_position_max(q, np.array([0, 0]), 1)
        self.assertEqual(out.tolist(), [[0.7, 0.2]])

    def test_a_position_with_only_missing_values_stays_missing(self):
        out = self.m.per_position_max(np.array([[np.nan]]), np.array([0]), 1)
        self.assertTrue(np.isnan(out[0, 0]))

class TestProfileStats(unittest.TestCase):
    def setUp(self):
        self.m = load("43_saturation_features.py")

    def test_basic_summaries(self):
        pos = np.arange(1, 101)
        prof = np.zeros(100)
        prof[10:20] = 1.0                  # 10 bp of signal inside one 50-bp window
        s = self.m.profile_stats(prof, pos, 0.25)
        self.assertAlmostEqual(s["mean"], 0.1)
        self.assertEqual(s["max"], 1.0)
        self.assertEqual(s["sum"], 10.0)
        self.assertEqual(s["share_best50"], 1.0)
        self.assertAlmostEqual(s["frac_above_high"], 0.1)

    def test_spread_signal_has_a_low_window_share(self):
        pos = np.arange(1, 201)
        s = self.m.profile_stats(np.ones(200), pos, 0.25)
        self.assertAlmostEqual(s["share_best50"], 0.25)

    def test_the_window_respects_genomic_gaps(self):
        """Two positions 60 bp apart cannot share one 50-bp window even if adjacent in the array."""
        s = self.m.profile_stats(np.array([1.0, 1.0]), np.array([1, 61]), 0.25)
        self.assertEqual(s["share_best50"], 0.5)

    def test_an_all_missing_profile_gives_missing_summaries(self):
        s = self.m.profile_stats(np.array([np.nan, np.nan]), np.array([1, 2]), 0.25)
        self.assertTrue(all(v != v for v in s.values()))


class TestTrackGroups(unittest.TestCase):
    def setUp(self):
        self.m = load("43_saturation_features.py")
        self.var = {"ontology_curie": ["UBERON:0002107", "EFO:0001187", "CL:0000182", "UBERON:0001114"],
                    "biosample_name": ["liver", "HepG2", "hepatocyte", "right lobe of liver"],
                    "biosample_type": ["tissue", "cell_line", "in_vitro_differentiated_cells", "tissue"],
                    "histone_mark": ["H3K27ac", "H3K27ac", "H3K27ac", "H3K4me1"],
                    "transcription_factor": ["HNF4A", "FOXA2", "EZH2", "CTCF"],
                    "name": ["t0", "t1", "t2", "t3"],
                    "nonzero_mean": ["2.0", "1.0", "1.0", "4.0"]}

    def test_adult_liver_and_hepg2_are_separate_and_hepatocyte_is_in_neither(self):
        g = self.m.group_columns("CHIP_TF", self.var)
        self.assertEqual(g["chip_tf_liver"].tolist(), [0, 3])
        self.assertEqual(g["chip_tf_hepg2"].tolist(), [1])

    def test_h3k27ac_requires_the_mark_and_liver(self):
        self.assertEqual(self.m.group_columns("CHIP_HISTONE", self.var)["liver_h3k27ac"].tolist(), [0])

    def test_an_unknown_scorer_is_refused(self):
        with self.assertRaises(self.m.la.ContractError):
            self.m.group_columns("RNA_SEQ", self.var)

    def test_summarise_ranks_tf_tracks_by_mean_sensitivity(self):
        variants = ["chr1:11:A>C", "chr1:11:A>G", "chr1:12:C>T"]
        q = np.array([[0.1, -0.9, 0.0, 0.2], [-0.4, 0.1, 0.0, 0.3], [0.2, 0.5, 0.9, 0.1]])
        feats, top = self.m.summarise_scorer("CHIP_TF", variants, q * 2.0, q, self.var, 10, 20)
        liver = [t for t in top if t["group"] == "chip_tf_liver"]
        # track 0: positions max |q| 0.4 and 0.2 -> 0.3; track 3: 0.3 and 0.1 -> 0.2
        self.assertEqual([t["transcription_factor"] for t in liver], ["HNF4A", "CTCF"])
        self.assertAlmostEqual(liver[0]["mean_max_abs_quantile"], 0.3)
        # liver profile = mean over tracks 0 and 3 per position: (0.4+0.3)/2, (0.2+0.1)/2
        self.assertAlmostEqual(feats["chip_tf_liver_q_mean"], 0.25)
        self.assertEqual(feats["chip_tf_liver_n_tracks"], 2)
        # relative scale: raw = 2q, divided by track means 2.0 and 4.0 -> |q| and |q|/2
        self.assertAlmostEqual(liver[0]["mean_rel_effect"], 0.3)
        self.assertAlmostEqual(liver[1]["mean_rel_effect"], 0.1)
        self.assertAlmostEqual(feats["chip_tf_liver_rel_mean"], 0.2)

    def test_a_track_without_a_usable_mean_level_is_missing_not_infinite(self):
        rel = self.m.relative_effects(np.array([[1.0, 1.0, 1.0]]), {"nonzero_mean": ["0", "", "2.0"]})
        self.assertTrue(np.isnan(rel[0, 0]) and np.isnan(rel[0, 1]))
        self.assertEqual(rel[0, 2], 0.5)


if __name__ == "__main__":
    unittest.main()

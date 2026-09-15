"""Fixture tests for step 68b, the independent checks of step 68's indel orientations."""

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


class TestStudyConvention(unittest.TestCase):
    def setUp(self):
        self.m = load("68b_check_dbsnp_orientation.py")

    def test_a_study_listing_alt_first_is_labelled(self):
        self.assertEqual(self.m.study_convention(n_allele1_is_ref=5, n_snvs=1000), "allele1_is_alt")

    def test_a_study_listing_ref_first_is_labelled(self):
        self.assertEqual(self.m.study_convention(n_allele1_is_ref=995, n_snvs=1000), "allele1_is_ref")

    def test_effect_allele_ordering_is_mixed(self):
        self.assertEqual(self.m.study_convention(n_allele1_is_ref=300, n_snvs=1000), "mixed")

    def test_a_small_study_is_not_labelled(self):
        self.assertEqual(self.m.study_convention(n_allele1_is_ref=0, n_snvs=199), "too_few")

    def test_the_boundaries_are_inclusive(self):
        self.assertEqual(self.m.study_convention(n_allele1_is_ref=10, n_snvs=1000), "allele1_is_alt")
        self.assertEqual(self.m.study_convention(n_allele1_is_ref=990, n_snvs=1000), "allele1_is_ref")


class TestStudyFamily(unittest.TestCase):
    def test_pan_ukbb_is_not_folded_into_ukbb(self):
        m = load("68b_check_dbsnp_orientation.py")
        self.assertEqual(m.study_family("PanUKBB_CSA_GGT"), "PanUKBB")
        self.assertEqual(m.study_family("UKBB_ALT"), "UKBB")
        self.assertEqual(m.study_family("MVP_NAFLD_EUR"), "MVP")
        self.assertEqual(m.study_family("2021_34128465_PDFF_EUR"), "other")


class TestTruthFromConventions(unittest.TestCase):
    def setUp(self):
        self.m = load("68b_check_dbsnp_orientation.py")

    def test_alt_first_puts_the_reference_second(self):
        self.assertEqual(self.m.truth_orientation("T", "TA", {"allele1_is_alt"}), (("TA", "T"), "truth"))

    def test_ref_first_keeps_the_order(self):
        self.assertEqual(self.m.truth_orientation("T", "TA", {"allele1_is_ref"}), (("T", "TA"), "truth"))

    def test_mixed_and_small_studies_give_no_truth(self):
        self.assertEqual(self.m.truth_orientation("T", "TA", {"mixed", "too_few"}), (None, "no_convention"))

    def test_opposite_conventions_are_excluded(self):
        self.assertEqual(self.m.truth_orientation("T", "TA", {"allele1_is_alt", "allele1_is_ref", "mixed"}),
                         (None, "conflicting_conventions"))


class TestIndependentFreqParser(unittest.TestCase):
    def setUp(self):
        self.m = load("68b_check_dbsnp_orientation.py")

    def test_maximum_across_studies_per_alt(self):
        self.assertEqual(self.m.freq_per_alt("A:0.9,0.1,.|B:0.7,0.05,0.25", 2), [0.1, 0.25])

    def test_missing_field_gives_none(self):
        self.assertEqual(self.m.freq_per_alt(".", 2), [None, None])

    def test_a_short_study_covers_only_its_alleles(self):
        self.assertEqual(self.m.freq_per_alt("A:0.9,0.1", 3), [0.1, None, None])


class TestSiteClassification(unittest.TestCase):
    def setUp(self):
        self.m = load("68b_check_dbsnp_orientation.py")
        self.ours = {("chr1", 10, "T", "TT"): ("T", "TT"), ("chr1", 10, "TT", "T"): ("TT", "T")}

    def test_one_matching_record_resolves(self):
        got = self.m.classify_site(self.ours, {("chr1", 10, "T", "TT"): [None]})
        self.assertEqual(got, ("unique_record", "T", "TT"))

    def test_both_matching_with_one_common_resolves_by_frequency(self):
        got = self.m.classify_site(self.ours, {("chr1", 10, "T", "TT"): [0.003],
                                               ("chr1", 10, "TT", "T"): [0.2, 0.31]})
        self.assertEqual(got, ("frequency", "TT", "T"))

    def test_both_common_is_unresolved(self):
        got = self.m.classify_site(self.ours, {("chr1", 10, "T", "TT"): [0.02], ("chr1", 10, "TT", "T"): [0.2]})
        self.assertEqual(got, ("both_common", None, None))

    def test_a_missing_frequency_is_not_common(self):
        got = self.m.classify_site(self.ours, {("chr1", 10, "T", "TT"): [None], ("chr1", 10, "TT", "T"): [None]})
        self.assertEqual(got, ("neither_common", None, None))

    def test_no_match_is_absent(self):
        self.assertEqual(self.m.classify_site(self.ours, {("chr1", 11, "A", "AT"): [0.4]}), ("absent", None, None))


class TestStep68StateVocabulary(unittest.TestCase):
    def setUp(self):
        self.m = load("68b_check_dbsnp_orientation.py")

    def test_states_map_onto_the_check_vocabulary(self):
        base = {"dbsnp_hg38_ref": "", "dbsnp_hg38_alt": ""}
        self.assertEqual(self.m.step68_state({**base, "dbsnp_hg38_state": "resolved", "frequency_tie_break": "",
                                              "dbsnp_hg38_ref": "T", "dbsnp_hg38_alt": "TT"}),
                         ("unique_record", "T", "TT"))
        self.assertEqual(self.m.step68_state({**base, "dbsnp_hg38_state": "resolved",
                                              "frequency_tie_break": "resolved_by_frequency",
                                              "dbsnp_hg38_ref": "TT", "dbsnp_hg38_alt": "T"}),
                         ("frequency", "TT", "T"))
        self.assertEqual(self.m.step68_state({**base, "dbsnp_hg38_state": "both_in_dbsnp",
                                              "frequency_tie_break": "both_common"}), ("both_common", None, None))
        self.assertEqual(self.m.step68_state({**base, "dbsnp_hg38_state": "absent", "frequency_tie_break": ""}),
                         ("absent", None, None))

    def test_an_unknown_state_is_refused(self):
        with self.assertRaises(self.m.la.ContractError):
            self.m.step68_state({"dbsnp_hg38_state": "resolved", "frequency_tie_break": "whatever",
                                 "dbsnp_hg38_ref": "T", "dbsnp_hg38_alt": "TT"})


if __name__ == "__main__":
    unittest.main()

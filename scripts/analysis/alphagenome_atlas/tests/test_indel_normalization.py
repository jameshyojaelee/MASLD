"""Fixture tests for indel normalisation (step 69).

dbSNP collapses a homopolymer or short tandem repeat into ONE multi-allelic record whose REF spans the whole
repeat, so `C/CTT` at a position whose reference reads `CTTTTTTT` appears as rs3038218 REF=CTTTTTTT with
twenty ALTs. Exact string matching on (REF, ALT) therefore misses it and calls the variant absent from dbSNP.
Both sides must be reduced to a canonical left-aligned, trimmed form before comparison.
"""

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


class TestNormalise(unittest.TestCase):
    def setUp(self):
        self.m = load("69_normalise_indels.py")
        # 1-based: 1='A' 2..9='C' then T-run. Reference reads  A C T T T T T T T G  at 1..10
        self.seq = "AC" + "T" * 7 + "G" + "ACGT" * 8
        self.fetch = lambda chrom, s, e: self.seq[s:e]

    def n(self, pos, ref, alt):
        return self.m.normalise_indel("chr1", pos, ref, alt, self.fetch)

    def test_an_already_minimal_insertion_is_unchanged(self):
        self.assertEqual(self.n(2, "C", "CTT"), (2, "C", "CTT"))

    def test_a_common_suffix_is_trimmed(self):
        """C+7T -> C+9T is the same event as C -> CTT once the shared tail is removed."""
        self.assertEqual(self.n(2, "C" + "T" * 7, "C" + "T" * 9), (2, "C", "CTT"))

    def test_a_deletion_in_the_same_repeat_normalises_to_the_short_form(self):
        self.assertEqual(self.n(2, "C" + "T" * 7, "C" + "T" * 5), (2, "CTT", "C"))

    def test_the_canonical_form_is_left_aligned_not_right_aligned(self):
        """Verified against `bcftools norm -f ref.fa`, which returns (2, C, CT) for both spellings."""
        self.assertEqual(self.n(2, "CT", "CTT"), (2, "C", "CT"))

    def test_a_right_aligned_spelling_converges_on_the_same_left_aligned_form(self):
        self.assertEqual(self.n(3, "T", "TT"), self.n(2, "CT", "CTT"))

    def test_a_whole_repeat_deletion_keeps_its_span(self):
        self.assertEqual(self.n(2, "C" + "T" * 7, "C"), (2, "C" + "T" * 7, "C"))

    def test_an_snv_is_left_alone(self):
        self.assertEqual(self.n(2, "C", "G"), (2, "C", "G"))

    def test_normalisation_is_idempotent(self):
        once = self.n(2, "C" + "T" * 7, "C" + "T" * 9)
        self.assertEqual(self.m.normalise_indel("chr1", *once, self.fetch), once)

    def test_identical_ref_and_alt_are_refused(self):
        with self.assertRaises(Exception):
            self.n(2, "CT", "CT")

    def test_an_empty_allele_is_refused(self):
        with self.assertRaises(Exception):
            self.n(2, "C", "")


class TestMatchAgainstDbsnpRecord(unittest.TestCase):
    """The real failure: our C/CTT and dbSNP's CTTTTTTT-with-twenty-ALTs describe the same events."""

    def setUp(self):
        self.m = load("69_normalise_indels.py")
        self.seq = "AC" + "T" * 7 + "G" + "ACGT" * 8
        self.fetch = lambda chrom, s, e: self.seq[s:e]
        # rs3038218-shaped: REF spans the whole 7-T run, ALTs enumerate other run lengths
        self.rec = (2, "C" + "T" * 7, ["C" + "T" * n for n in (0, 2, 5, 9, 12)])

    def test_our_insertion_matches_a_longer_run_allele(self):
        out = self.m.orientation_from_normalised([self.rec], "chr1", 2, "C", "CTT", self.fetch)
        self.assertEqual(out["state"], "both_in_dbsnp")

    def test_a_record_holding_only_the_longer_run_resolves_the_insertion(self):
        rec = (2, "C" + "T" * 7, ["C" + "T" * 9])
        out = self.m.orientation_from_normalised([rec], "chr1", 2, "C", "CTT", self.fetch)
        self.assertEqual(out["state"], "resolved")
        self.assertEqual((out["ref"], out["alt"]), ("C", "CTT"))

    def test_a_record_holding_only_the_shorter_run_resolves_the_deletion(self):
        rec = (2, "C" + "T" * 7, ["C" + "T" * 5])
        out = self.m.orientation_from_normalised([rec], "chr1", 2, "C", "CTT", self.fetch)
        self.assertEqual(out["state"], "resolved")
        self.assertEqual((out["ref"], out["alt"]), ("CTT", "C"))

    def test_an_unrelated_record_leaves_the_variant_absent(self):
        out = self.m.orientation_from_normalised([(2, "C", ["G"])], "chr1", 2, "C", "CTT", self.fetch)
        self.assertEqual(out["state"], "absent")

    def test_a_plain_non_repeat_record_still_matches_exactly(self):
        out = self.m.orientation_from_normalised([(2, "C", ["CTT"])], "chr1", 2, "C", "CTT", self.fetch)
        self.assertEqual(out["state"], "resolved")
        self.assertEqual((out["ref"], out["alt"]), ("C", "CTT"))


if __name__ == "__main__":
    unittest.main()

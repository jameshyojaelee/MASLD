"""Fixture tests for step 68: resolving indel allele orientation against dbSNP.

The reference cannot orient a prefix-anchored indel (both readings are representable), and the GWAS source
id carries effect/other alleles in no order. dbSNP stores each variant in its canonical REF/ALT form, so a
record whose allele pair matches ours names which allele is the reference.
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


class TestRefSeqAccessions(unittest.TestCase):
    def setUp(self):
        self.m = load("68_dbsnp_orientation.py")

    def test_grch37_autosome_accession_maps_to_a_chr_name(self):
        self.assertEqual(self.m.chrom_of("NC_000001.10"), "chr1")

    def test_grch38_accession_maps_to_the_same_chr_name(self):
        self.assertEqual(self.m.chrom_of("NC_000001.11"), "chr1")

    def test_sex_chromosomes_map(self):
        self.assertEqual(self.m.chrom_of("NC_000023.10"), "chrX")
        self.assertEqual(self.m.chrom_of("NC_000024.10"), "chrY")

    def test_an_unplaced_scaffold_is_not_invented_into_a_chromosome(self):
        self.assertIsNone(self.m.chrom_of("NT_167214.1"))


class TestOrientationFromRecords(unittest.TestCase):
    """dbSNP may hold one orientation, the other, both, or neither at a position."""

    def setUp(self):
        self.m = load("68_dbsnp_orientation.py")

    def test_a_single_matching_record_resolves_the_orientation(self):
        out = self.m.orientation_from_records([("T", ["TA"])], "T", "TA")
        self.assertEqual((out["ref"], out["alt"]), ("T", "TA"))
        self.assertEqual(out["state"], "resolved")

    def test_the_pair_is_matched_regardless_of_the_order_we_hold_it_in(self):
        out = self.m.orientation_from_records([("T", ["TA"])], "TA", "T")
        self.assertEqual((out["ref"], out["alt"]), ("T", "TA"))
        self.assertEqual(out["state"], "resolved")

    def test_a_multiallelic_record_still_matches_on_the_relevant_alt(self):
        out = self.m.orientation_from_records([("T", ["TA", "TAA", "G"])], "T", "TA")
        self.assertEqual((out["ref"], out["alt"]), ("T", "TA"))

    def test_both_orientations_present_is_reported_as_both_not_picked(self):
        out = self.m.orientation_from_records([("T", ["TA"]), ("TA", ["T"])], "T", "TA")
        self.assertEqual(out["state"], "both_in_dbsnp")
        self.assertIsNone(out["ref"])

    def test_no_matching_record_is_absent(self):
        out = self.m.orientation_from_records([("C", ["G"])], "T", "TA")
        self.assertEqual(out["state"], "absent")
        self.assertIsNone(out["ref"])

    def test_an_empty_record_list_is_absent(self):
        self.assertEqual(self.m.orientation_from_records([], "T", "TA")["state"], "absent")

    def test_alleles_are_compared_case_insensitively(self):
        out = self.m.orientation_from_records([("t", ["ta"])], "T", "TA")
        self.assertEqual((out["ref"], out["alt"]), ("T", "TA"))


class TestAgreementWithTheReferenceRule(unittest.TestCase):
    """dbSNP is the arbiter, but a disagreement with the reference-validated orientation is a defect
    somewhere and must be surfaced, never silently overwritten."""

    def setUp(self):
        self.m = load("68_dbsnp_orientation.py")

    def test_agreement_is_reported(self):
        v = self.m.compare_to_current("chr4:87310240:T:TA", "T", "TA")
        self.assertEqual(v, "agrees")

    def test_a_flip_relative_to_the_current_call_is_reported(self):
        v = self.m.compare_to_current("chr4:87310240:TA:T", "T", "TA")
        self.assertEqual(v, "flips")

    def test_a_different_allele_pair_is_not_called_agreement(self):
        """The function is given alleles, never a position, so this is the only mismatch it can detect;
        in the pipeline the site is identical by construction."""
        v = self.m.compare_to_current("chr4:87310240:C:G", "T", "TA")
        self.assertEqual(v, "different_alleles")

    def test_no_current_call_is_reported_as_such(self):
        self.assertEqual(self.m.compare_to_current("", "T", "TA"), "no_current_call")


if __name__ == "__main__":
    unittest.main()


class TestFreqParsing(unittest.TestCase):
    """dbSNP's FREQ is `STUDY:refAF,alt1AF,...|STUDY2:...` with '.' for unreported alleles. In a repeat
    record both of our candidate orientations exist as ALTs, so the allele FREQUENCIES are what say which
    one a GWAS could have tested: a GWAS indel is common by construction."""

    def setUp(self):
        self.m = load("68_dbsnp_orientation.py")

    def test_frequencies_are_indexed_by_alt_position(self):
        f = self.m.parse_freq("1000Genomes:0.6354,.,.,.,.,.,0.3646", n_alts=6)
        self.assertAlmostEqual(f[5], 0.3646, places=6)
        self.assertIsNone(f[0])

    def test_several_studies_are_combined_by_maximum(self):
        f = self.m.parse_freq("1000Genomes:0.9,.,0.1|TOPMED:0.8,0.05,.", n_alts=2)
        self.assertAlmostEqual(f[0], 0.05, places=6)
        self.assertAlmostEqual(f[1], 0.1, places=6)

    def test_a_missing_freq_field_yields_no_frequencies(self):
        self.assertEqual(self.m.parse_freq("", n_alts=3), [None, None, None])

    def test_a_ragged_study_entry_does_not_crash_or_shift_indices(self):
        f = self.m.parse_freq("A:0.5,0.5|B:0.4", n_alts=2)
        self.assertAlmostEqual(f[0], 0.5, places=6)
        self.assertIsNone(f[1])


class TestFrequencyTieBreak(unittest.TestCase):
    def setUp(self):
        self.m = load("68_dbsnp_orientation.py")

    def test_the_common_allele_wins(self):
        out = self.m.tie_break_by_frequency({("C", "CT"): 0.3646, ("CT", "C"): 3.8e-06}, floor=0.01)
        self.assertEqual(out["state"], "resolved_by_frequency")
        self.assertEqual((out["ref"], out["alt"]), ("C", "CT"))

    def test_two_common_alleles_stay_ambiguous(self):
        out = self.m.tie_break_by_frequency({("C", "CT"): 0.30, ("CT", "C"): 0.25}, floor=0.01)
        self.assertEqual(out["state"], "both_common")

    def test_two_rare_alleles_stay_ambiguous(self):
        out = self.m.tie_break_by_frequency({("C", "CT"): 1e-5, ("CT", "C"): 2e-5}, floor=0.01)
        self.assertEqual(out["state"], "neither_common")

    def test_a_missing_frequency_counts_as_not_common_not_as_zero_evidence(self):
        out = self.m.tie_break_by_frequency({("C", "CT"): 0.3, ("CT", "C"): None}, floor=0.01)
        self.assertEqual(out["state"], "resolved_by_frequency")
        self.assertEqual((out["ref"], out["alt"]), ("C", "CT"))

    def test_no_frequencies_at_all_is_not_a_resolution(self):
        out = self.m.tie_break_by_frequency({("C", "CT"): None, ("CT", "C"): None}, floor=0.01)
        self.assertEqual(out["state"], "neither_common")


class TestCrossAssemblyAgreement(unittest.TestCase):
    """A non-answer is not a disagreement. The hg19 arm cannot left-shift (no hg19 FASTA) and has no
    frequency tie-break, so it returns `both_in_dbsnp` wherever hg38 resolves by frequency. Treating that
    as DISAGREE vetoed 66,859 resolutions the primary arm had made."""

    def setUp(self):
        self.m = load("68_dbsnp_orientation.py")

    def R(self, ref="C", alt="CT"):
        return {"state": "resolved", "ref": ref, "alt": alt}

    def U(self, state="both_in_dbsnp"):
        return {"state": state, "ref": None, "alt": None}

    def test_both_resolved_and_equal_is_agreement(self):
        a, c = self.m.reconcile(self.R(), self.R(), True)
        self.assertEqual(a, "agree")
        self.assertEqual((c["ref"], c["alt"]), ("C", "CT"))

    def test_both_resolved_and_different_is_a_refusal(self):
        a, c = self.m.reconcile(self.R("C", "CT"), self.R("CT", "C"), True)
        self.assertEqual(a, "DISAGREE")
        self.assertIsNone(c)

    def test_hg38_resolved_and_hg19_silent_keeps_the_hg38_answer(self):
        a, c = self.m.reconcile(self.U(), self.R(), True)
        self.assertEqual(a, "one_resolved_only")
        self.assertEqual((c["ref"], c["alt"]), ("C", "CT"))

    def test_hg19_resolved_and_hg38_silent_keeps_the_hg19_answer(self):
        a, c = self.m.reconcile(self.R(), self.U("absent"), True)
        self.assertEqual(a, "one_resolved_only")
        self.assertEqual((c["ref"], c["alt"]), ("C", "CT"))

    def test_neither_resolved_yields_nothing(self):
        a, c = self.m.reconcile(self.U(), self.U(), True)
        self.assertEqual(a, "neither_resolved")
        self.assertIsNone(c)

    def test_without_the_hg38_arm_nothing_is_cross_checked(self):
        a, c = self.m.reconcile(self.R(), self.U("absent"), False)
        self.assertEqual(a, "not_checked")
        self.assertEqual((c["ref"], c["alt"]), ("C", "CT"))

    def test_the_unresolved_label_names_the_ARBITER_state_not_the_confirming_arm(self):
        """The label was built from the hg19 state even after hg38 became the arbiter, which made the
        summary counts describe the wrong arm. hg38 is the arbiter, so its state names the outcome."""
        self.assertEqual(self.m.unresolved_label(self.U("both_in_dbsnp"), self.U("absent")),
                         "unresolved_absent")
        self.assertEqual(self.m.unresolved_label(self.U("absent"), self.U("both_in_dbsnp")),
                         "unresolved_both_in_dbsnp")

    def test_the_label_falls_back_to_the_confirming_arm_when_the_arbiter_resolved(self):
        self.assertEqual(self.m.unresolved_label(self.U("both_in_dbsnp"), self.R()),
                         "unresolved_both_in_dbsnp")

    def test_two_resolved_arms_that_differ_are_labelled_as_a_refusal(self):
        """The only way both arms resolve and nothing is chosen is a refusal; the run of 2026-09-14 wrote
        these 3 rows as `unresolved_resolved`, which names no state at all."""
        self.assertEqual(self.m.unresolved_label(self.R("C", "CT"), self.R("CT", "C")),
                         "refused_assemblies_disagree")

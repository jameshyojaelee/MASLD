"""Fixture tests for step 63: score, through the model API, the indels the Atlas point query refuses.

The Atlas serves neither route for indels - point queries return UNIMPLEMENTED and interval saturation
returns only the three substitutions per position (0 multi-base alleles in 104,352 variants over 40
archives). The model API takes an arbitrary sequence, so an indel can be built into the window instead.
`predict_sequence` requires exactly 1 Mb, so the window must stay that length whatever the indel does to it.
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


GENOME = "ACGT" * 64          # 256 bases, deterministic


def fetch(chrom, start0, end):
    """Mimics pysam FastaFile.fetch: 0-based, half-open, over a repeating reference."""
    return "".join(GENOME[i % len(GENOME)] for i in range(start0, end))


class TestIndelWindow(unittest.TestCase):
    def setUp(self):
        self.m = load("63_indel_rescue.py")

    def test_insertion_keeps_the_window_length(self):
        out = self.m.indel_alt_sequence(fetch, "chr1", 0, 64, pos1=10, ref="C", alt="CTTT")
        self.assertEqual(len(out), 64)

    def test_insertion_places_the_inserted_bases_at_the_variant(self):
        out = self.m.indel_alt_sequence(fetch, "chr1", 0, 64, pos1=10, ref="C", alt="CTTT")
        ref = fetch("chr1", 0, 64)
        self.assertEqual(out[:9], ref[:9])
        self.assertEqual(out[9:13], "CTTT")

    def test_insertion_trims_from_the_far_end_not_the_variant(self):
        out = self.m.indel_alt_sequence(fetch, "chr1", 0, 64, pos1=10, ref="C", alt="CTTT")
        ref = fetch("chr1", 0, 64)
        self.assertEqual(out[13:], ref[10:61])

    def test_deletion_keeps_the_window_length_by_pulling_in_downstream_bases(self):
        out = self.m.indel_alt_sequence(fetch, "chr1", 0, 64, pos1=10, ref="CGT", alt="C")
        self.assertEqual(len(out), 64)
        ref_ext = fetch("chr1", 0, 66)
        self.assertEqual(out, ref_ext[:10] + ref_ext[12:66])

    def test_a_reference_mismatch_is_refused_rather_than_silently_scored(self):
        with self.assertRaises(Exception):
            self.m.indel_alt_sequence(fetch, "chr1", 0, 64, pos1=10, ref="AAAA", alt="A")

    def test_an_snv_through_this_path_changes_exactly_one_base(self):
        ref = fetch("chr1", 0, 64)
        out = self.m.indel_alt_sequence(fetch, "chr1", 0, 64, pos1=10, ref=ref[9], alt="G")
        self.assertEqual(len(out), 64)
        self.assertEqual(sum(a != b for a, b in zip(ref, out)), 1)

    def test_the_reference_window_is_returned_unchanged_for_a_null_edit(self):
        ref = fetch("chr1", 0, 64)
        self.assertEqual(self.m.indel_alt_sequence(fetch, "chr1", 0, 64, 10, ref[9], ref[9]), ref)


class TestIndelClassification(unittest.TestCase):
    def setUp(self):
        self.m = load("63_indel_rescue.py")

    def test_insertion_and_deletion_are_named(self):
        self.assertEqual(self.m.indel_class("T", "TA"), "insertion")
        self.assertEqual(self.m.indel_class("TA", "T"), "deletion")
        self.assertEqual(self.m.indel_class("T", "A"), "snv")

    def test_length_change_is_signed(self):
        self.assertEqual(self.m.length_delta("T", "TAA"), 2)
        self.assertEqual(self.m.length_delta("TAA", "T"), -2)


class TestSourceParsing(unittest.TestCase):
    def setUp(self):
        self.m = load("63_indel_rescue.py")

    def test_both_source_id_shapes_parse(self):
        """signal_variant_weights carries `hg19:chr4:88231392:T:TA` and `hg19:4:88231392:T:TA`."""
        for sid in ("hg19:chr4:88231392:T:TA", "hg19:4:88231392:T:TA"):
            got = self.m.parse_source_id(sid)
            self.assertEqual(got["chrom"], "chr4")
            self.assertEqual(got["pos"], 88231392)
            self.assertEqual(got["ref"], "T")
            self.assertEqual(got["alt"], "TA")

    def test_a_malformed_id_is_refused(self):
        with self.assertRaises(Exception):
            self.m.parse_source_id("hg19:chr4:88231392")


class TestResolutionFolding(unittest.TestCase):
    """CHIP_HISTONE is 128-bp binned; a base-resolution mask applied to it silently selects the wrong rows.

    This is the same defect that hit P5A. The mask must be folded to the track's own row count first.
    """

    def setUp(self):
        self.m = load("63_indel_rescue.py")

    def test_a_base_mask_folds_to_the_binned_row_count(self):
        import numpy as np
        mask = np.zeros(1024, bool)
        mask[300:400] = True
        folded = self.m.fold_mask(mask, 8)              # 1024 bases -> 8 rows of 128
        self.assertEqual(folded.size, 8)
        self.assertTrue(folded[2] and folded[3])
        self.assertFalse(folded[0] or folded[7])

    def test_a_mask_already_at_the_track_resolution_is_unchanged(self):
        import numpy as np
        mask = np.array([True, False, True, False])
        self.assertTrue((self.m.fold_mask(mask, 4) == mask).all())


class TestHistoneTrackSelection(unittest.TestCase):
    def setUp(self):
        self.m = load("63_indel_rescue.py")

    def test_only_h3k27ac_tracks_are_kept_from_the_histone_output(self):
        names = ["liver H3K27ac", "liver H3K4me1", "liver H3K27ac rep2", "liver H3K9me3"]
        keep = self.m.h3k27ac_columns(names)
        self.assertEqual(list(keep), [True, False, True, False])

    def test_a_histone_output_with_no_h3k27ac_track_selects_nothing(self):
        self.assertFalse(any(self.m.h3k27ac_columns(["liver H3K4me3"])))


class TestTargetSelection(unittest.TestCase):
    def setUp(self):
        self.m = load("63_indel_rescue.py")

    def test_targets_are_ordered_by_posterior_mass(self):
        mass = {"a": 0.02, "b": 0.9, "c": 0.2}
        self.assertEqual(self.m.select_targets(mass, min_mass=0.01, only=""), ["b", "c", "a"])

    def test_mass_below_the_threshold_is_dropped(self):
        self.assertEqual(self.m.select_targets({"a": 0.001, "b": 0.5}, min_mass=0.01, only=""), ["b"])

    def test_an_explicit_list_overrides_the_threshold_but_not_membership(self):
        mass = {"a": 0.001, "b": 0.5}
        self.assertEqual(self.m.select_targets(mass, min_mass=0.01, only="a"), ["a"])
        self.assertEqual(self.m.select_targets(mass, min_mass=0.01, only="a,zzz"), ["a"])

    def test_an_explicit_member_id_selects_its_group_once(self):
        mass = {"rep": 0.5}
        got = self.m.select_targets(mass, min_mass=0.01, only="member,rep", rep_of={"member": "rep", "rep": "rep"})
        self.assertEqual(got, ["rep"])


class TestOneVariantManyIds(unittest.TestCase):
    """Studies write one biallelic indel with its alleles in either order and with or without 'chr', so
    `T:TCA` and `TCA:T` at one site are ONE variant. Scoring each id as its own target counted 8 variants
    twice among 152 targets and scored 3 of them as both reciprocal variants."""

    def setUp(self):
        self.m = load("63_indel_rescue.py")
        self.xw = {"hg19:10:5:T:TCA": {"hg38_chrom": "chr10", "hg38_position_1based": "7"},
                   "hg19:10:5:TCA:T": {"hg38_chrom": "chr10", "hg38_position_1based": "7"},
                   "hg19:chr10:5:T:TCA": {"hg38_chrom": "chr10", "hg38_position_1based": "7"},
                   "hg19:10:5:T:TC": {"hg38_chrom": "chr10", "hg38_position_1based": "7"},
                   "hg19:10:9:G:GA": {"hg38_chrom": "", "hg38_position_1based": ""},
                   "hg19:10:9:GA:G": {"hg38_chrom": "", "hg38_position_1based": ""}}
        self.mass = {"hg19:10:5:T:TCA": 0.03, "hg19:10:5:TCA:T": 0.04, "hg19:chr10:5:T:TCA": 0.01,
                     "hg19:10:5:T:TC": 0.2, "hg19:10:9:G:GA": 0.1, "hg19:10:9:GA:G": 0.1}
        self.sigs = {k: {f"sig_{k}"} for k in self.mass}

    def test_ids_naming_one_variant_share_a_group_with_summed_mass_and_signals(self):
        g = self.m.variant_groups(self.mass, self.sigs, self.xw)
        rep = "hg19:10:5:TCA:T"                      # the heaviest member
        self.assertAlmostEqual(g[rep]["mass"], 0.08)
        self.assertEqual(g[rep]["members"], ["hg19:10:5:TCA:T", "hg19:10:5:T:TCA", "hg19:chr10:5:T:TCA"])
        self.assertEqual(len(g[rep]["signals"]), 3)

    def test_a_different_allele_pair_at_the_same_site_is_a_different_variant(self):
        g = self.m.variant_groups(self.mass, self.sigs, self.xw)
        self.assertEqual(g["hg19:10:5:T:TC"]["members"], ["hg19:10:5:T:TC"])

    def test_ids_without_an_hg38_coordinate_are_never_merged(self):
        """Without a coordinate there is nothing to establish they are the same variant."""
        g = self.m.variant_groups(self.mass, self.sigs, self.xw)
        self.assertIn("hg19:10:9:G:GA", g)
        self.assertIn("hg19:10:9:GA:G", g)

    def test_every_id_lands_in_exactly_one_group(self):
        g = self.m.variant_groups(self.mass, self.sigs, self.xw)
        members = [m for v in g.values() for m in v["members"]]
        self.assertEqual(sorted(members), sorted(self.mass))

    def test_one_dbsnp_answer_is_taken_for_the_group(self):
        db = {"a": {"dbsnp_ref": "T", "dbsnp_alt": "TCA", "dbsnp_rule": "frequency"},
              "b": {"dbsnp_ref": "T", "dbsnp_alt": "TCA", "dbsnp_rule": "unique_record"}}
        got = self.m.group_dbsnp_orientation(["a", "b", "c"], db)
        self.assertEqual((got["dbsnp_ref"], got["dbsnp_alt"], got["dbsnp_rule"]), ("T", "TCA", "unique_record"))

    def test_members_dbsnp_disagrees_about_are_refused(self):
        db = {"a": {"dbsnp_ref": "T", "dbsnp_alt": "TCA", "dbsnp_rule": "unique_record"},
              "b": {"dbsnp_ref": "TCA", "dbsnp_alt": "T", "dbsnp_rule": "unique_record"}}
        with self.assertRaises(self.m.la.ContractError):
            self.m.group_dbsnp_orientation(["a", "b"], db)

    def test_a_group_dbsnp_does_not_cover_has_no_answer(self):
        self.assertIsNone(self.m.group_dbsnp_orientation(["a"], {}))


class TestPairedComparison(unittest.TestCase):
    """Is the excluded variant class more impactful than the included one AT THE SAME SIGNALS?

    The median ratio over all indels and all controls is not a test: it compares two differently-composed
    sets. Each indel is paired with the SNV controls drawn at its own signals, the difference is taken per
    indel, and 1-Mb blocks are the resampling unit.
    """

    def setUp(self):
        self.m = load("63_indel_rescue.py")

    def test_pairs_are_built_from_the_control_for_key(self):
        indels = [{"source_variant_id": "i1", "chrom": "chr1", "pos_hg38": 1_000, "rna_log2": 0.4}]
        ctls = [{"control_for": "i1", "chrom": "chr1", "pos_hg38": 1_100, "rna_log2": 0.1},
                {"control_for": "i1", "chrom": "chr1", "pos_hg38": 1_200, "rna_log2": 0.3},
                {"control_for": "other", "chrom": "chr1", "pos_hg38": 1_300, "rna_log2": 9.9}]
        got = self.m.paired_rows(indels, ctls, "rna")
        self.assertEqual(len(got), 1)
        self.assertAlmostEqual(got[0]["value"], 0.4 - 0.2, places=12)   # |0.4| - median(|0.1|,|0.3|)

    def test_an_indel_with_no_control_is_excluded_not_scored_against_zero(self):
        indels = [{"source_variant_id": "i1", "chrom": "chr1", "pos_hg38": 1_000, "rna_log2": 0.4}]
        self.assertEqual(self.m.paired_rows(indels, [], "rna"), [])

    def test_the_comparison_is_on_absolute_effect_so_sign_does_not_cancel(self):
        indels = [{"source_variant_id": "i1", "chrom": "chr1", "pos_hg38": 1_000, "rna_log2": -0.5}]
        ctls = [{"control_for": "i1", "chrom": "chr1", "pos_hg38": 1_100, "rna_log2": 0.1}]
        self.assertAlmostEqual(self.m.paired_rows(indels, ctls, "rna")[0]["value"], 0.4, places=12)

    def test_nan_effects_are_dropped_rather_than_propagated(self):
        indels = [{"source_variant_id": "i1", "chrom": "chr1", "pos_hg38": 1_000, "rna_log2": float("nan")}]
        ctls = [{"control_for": "i1", "chrom": "chr1", "pos_hg38": 1_100, "rna_log2": 0.1}]
        self.assertEqual(self.m.paired_rows(indels, ctls, "rna"), [])

    def test_one_mb_blocks_are_the_resampling_unit(self):
        rows = [{"analysis_block": "chr1:0", "value": 0.3}] * 5 + [{"analysis_block": "chr2:7", "value": 0.1}]
        res = self.m.block_sign_test(rows)
        self.assertEqual(res["n_blocks"], 2)
        self.assertEqual(res["n_blocks_positive"], 2)

    def test_a_split_decision_is_reported_as_such(self):
        rows = [{"analysis_block": f"chr1:{i}", "value": (1.0 if i % 2 else -1.0)} for i in range(10)]
        res = self.m.block_sign_test(rows)
        self.assertEqual(res["n_blocks_positive"], 5)
        self.assertGreater(res["p_two_sided"], 0.5)


class TestReferenceOrientationBeforeScoring(unittest.TestCase):
    """Step 63's first run refused 12 of 152 targets because it tested only the SOURCE id's ref field and
    never tried the other order. GWAS source ids carry effect/other alleles, not ref/alt, so the fallback is
    what recovers them; two of the 12 carry the largest posterior mass in the whole target list.

    The source order is kept wherever the reference permits it. For a prefix-anchored indel the reference
    usually permits BOTH orders, and choosing the longer allele instead would call HSD17B13 rs72613567 a
    deletion when it is a known T->TA insertion at the exon-6 donor."""

    def setUp(self):
        self.m = load("63_indel_rescue.py")
        self.seq = "ACGT" * 64          # 1-based 10 is 'C', 11 'G', 12 'T'
        self.fetch = lambda chrom, s, e: self.seq[s:e]

    def test_the_source_order_is_kept_when_the_reference_permits_it(self):
        ref, alt, swap, amb = self.m.oriented_alleles("chr1", 10, "C", "CGT", self.fetch)
        self.assertEqual((ref, alt), ("C", "CGT"))
        self.assertEqual(swap, "none")

    def test_alleles_written_the_wrong_way_round_fall_back_instead_of_being_refused(self):
        ref, alt, swap, amb = self.m.oriented_alleles("chr1", 11, "GG", "G", self.fetch)
        self.assertEqual((ref, alt), ("G", "GG"))
        self.assertEqual(swap, "swapped")
        self.assertFalse(amb)

    def test_a_prefix_ambiguity_is_reported_to_the_caller(self):
        """The caller has to be able to say which rows rest on the source order rather than the reference."""
        ref, alt, swap, amb = self.m.oriented_alleles("chr1", 10, "C", "CGT", self.fetch)
        self.assertTrue(amb)

    def test_a_variant_genuinely_absent_from_the_reference_is_still_refused(self):
        with self.assertRaises(Exception):
            self.m.oriented_alleles("chr1", 10, "A", "AGT", self.fetch)

    def test_an_snv_is_oriented_and_never_ambiguous(self):
        ref, alt, swap, amb = self.m.oriented_alleles("chr1", 10, "T", "C", self.fetch)
        self.assertEqual((ref, alt), ("C", "T"))
        self.assertEqual(swap, "swapped")
        self.assertFalse(amb)

    def test_the_oriented_insertion_round_trips_through_the_window_builder(self):
        ref, alt, _, _ = self.m.oriented_alleles("chr1", 10, "C", "CGT", self.fetch)
        out = self.m.indel_alt_sequence(self.fetch, "chr1", 0, 40, 10, ref, alt)
        self.assertEqual(len(out), 40)
        self.assertEqual(out, (self.seq[:9] + "CGT" + self.seq[10:])[:40])

    def test_a_fallback_deletion_round_trips_too(self):
        ref, alt, _, _ = self.m.oriented_alleles("chr1", 11, "GG", "G", self.fetch)
        out = self.m.indel_alt_sequence(self.fetch, "chr1", 0, 40, 11, ref, alt)
        self.assertEqual(len(out), 40)


class TestDbsnpOrientationTakesPrecedence(unittest.TestCase):
    """dbSNP stores each variant in its canonical REF/ALT form and its rsIDs identify the variant the GWAS
    actually tested, so where dbSNP resolves an indel it outranks both the source's allele order and the
    reference-only rule. Where dbSNP is silent the reference rule still applies, and the row says which."""

    def setUp(self):
        self.m = load("63_indel_rescue.py")
        self.seq = "ACGT" * 64
        self.fetch = lambda chrom, s, e: self.seq[s:e]
        self.dbsnp = {"hg19:1:99:C:CGT": {"dbsnp_ref": "CGT", "dbsnp_alt": "C"}}

    def test_a_dbsnp_resolved_variant_uses_the_dbsnp_orientation(self):
        ref, alt, swap, amb, src = self.m.orientation_for("hg19:1:99:C:CGT", "chr1", 10, "C", "CGT",
                                                          self.fetch, self.dbsnp)
        self.assertEqual((ref, alt), ("CGT", "C"))
        self.assertEqual(src, "dbsnp")

    def test_a_dbsnp_resolved_variant_still_reports_whether_the_REFERENCE_could_orient_it(self):
        """The prefix ambiguity is a property of the alleles and the genome, not of dbSNP. Zeroing the flag
        for dbSNP rows made the `orientation_resolved` stratum identical to `dbsnp_oriented` -- a stratum
        that reports the same numbers as another is not an independent check."""
        *_, amb, src = self.m.orientation_for("hg19:1:99:C:CGT", "chr1", 10, "C", "CGT",
                                              self.fetch, self.dbsnp)
        self.assertEqual(src, "dbsnp")
        self.assertTrue(amb)        # chr1:10 reads CGT, so both readings are representable

    def test_a_dbsnp_resolved_variant_the_reference_can_also_orient_is_not_flagged(self):
        ref, alt, swap, amb, src = self.m.orientation_for("hg19:1:110:GG:G", "chr1", 11, "GG", "G",
                                                          self.fetch, {"hg19:1:110:GG:G":
                                                          {"dbsnp_ref": "G", "dbsnp_alt": "GG"}})
        self.assertEqual(src, "dbsnp")
        self.assertFalse(amb)

    def test_a_variant_dbsnp_does_not_cover_falls_back_to_the_reference_rule(self):
        ref, alt, swap, amb, src = self.m.orientation_for("hg19:1:99:C:CGT", "chr1", 10, "C", "CGT",
                                                          self.fetch, {})
        self.assertEqual((ref, alt), ("C", "CGT"))
        self.assertEqual(src, "reference_or_source_order")
        self.assertTrue(amb)

    def test_a_dbsnp_orientation_not_matching_our_allele_pair_is_refused_not_used(self):
        """A record for a different variant at the same site must never be applied to ours."""
        with self.assertRaises(Exception):
            self.m.orientation_for("hg19:1:99:C:CGT", "chr1", 10, "C", "CGT", self.fetch,
                                   {"hg19:1:99:C:CGT": {"dbsnp_ref": "A", "dbsnp_alt": "AT"}})

    def test_a_dbsnp_reference_allele_absent_from_the_genome_is_refused(self):
        with self.assertRaises(Exception):
            self.m.orientation_for("hg19:1:99:A:AT", "chr1", 10, "A", "AT", self.fetch,
                                   {"hg19:1:99:A:AT": {"dbsnp_ref": "A", "dbsnp_alt": "AT"}})


class TestDbsnpTableLoading(unittest.TestCase):
    """Step 68 writes one row per deferred indel with both assembly arms AND the reconciled answer. Only the
    reconciled answer may be applied: reading an arm's own columns would let a refused or overruled arm
    decide the orientation that gets scored."""

    def setUp(self):
        self.m = load("63_indel_rescue.py")

    def row(self, **kw):
        base = {"source_variant_id": "s1", "dbsnp_hg19_state": "resolved", "dbsnp_hg19_ref": "C",
                "dbsnp_hg19_alt": "CT", "dbsnp_hg38_state": "resolved", "dbsnp_hg38_ref": "C",
                "dbsnp_hg38_alt": "CT", "resolved_variant_uid": "chr1:10:C:CT", "frequency_tie_break": ""}
        base.update(kw)
        return base

    def test_alleles_come_from_the_reconciled_uid_not_from_an_arm(self):
        got = self.m.dbsnp_orientations([self.row(dbsnp_hg19_ref="CT", dbsnp_hg19_alt="C",
                                                  resolved_variant_uid="chr1:10:C:CT")])
        self.assertEqual((got["s1"]["dbsnp_ref"], got["s1"]["dbsnp_alt"]), ("C", "CT"))

    def test_a_row_with_no_reconciled_answer_is_not_applied(self):
        """Two resolved arms that disagree carry an empty uid; neither arm may be used."""
        got = self.m.dbsnp_orientations([self.row(resolved_variant_uid="",
                                                  dbsnp_hg38_ref="CT", dbsnp_hg38_alt="C")])
        self.assertEqual(got, {})

    def test_a_frequency_tie_break_is_labelled_as_such(self):
        got = self.m.dbsnp_orientations([self.row(frequency_tie_break="resolved_by_frequency")])
        self.assertEqual(got["s1"]["dbsnp_rule"], "frequency")

    def test_a_unique_record_resolution_is_labelled_as_such(self):
        for fb in ("", "both_common", "neither_common"):
            got = self.m.dbsnp_orientations([self.row(frequency_tie_break=fb)])
            self.assertEqual(got["s1"]["dbsnp_rule"], "unique_record", fb)

    def test_a_malformed_uid_is_refused(self):
        with self.assertRaises(self.m.la.ContractError):
            self.m.dbsnp_orientations([self.row(resolved_variant_uid="chr1:10:C")])


if __name__ == "__main__":
    unittest.main()

"""Offline unit tests for lib_atlas (protocol §10 fixtures, no network)."""

from __future__ import annotations

import math
import pathlib
import sys
import os
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import lib_atlas as la  # noqa: E402


class FakeFasta:
    """fetch(chrom, start0, end0) over a tiny synthetic reference."""

    def __init__(self, seqs: dict[str, str]):
        self.seqs = seqs

    def fetch(self, chrom, start0, end0):
        return self.seqs[chrom][start0:end0]


REF = FakeFasta({"chr1": "ACGTACGTAC", "chr2": "GGGGCCCCAA"})


class TestVariantUid(unittest.TestCase):
    def test_uid_includes_alleles_so_multiallelic_sites_stay_distinct(self):
        a = la.variant_uid("chr1", 5, "A", "G")
        b = la.variant_uid("chr1", 5, "A", "T")
        self.assertEqual(a, "chr1:5:A:G")
        self.assertNotEqual(a, b)

    def test_uid_rejects_missing_chr_prefix(self):
        with self.assertRaises(ValueError):
            la.variant_uid("1", 5, "A", "G")


class TestCrosswalk(unittest.TestCase):
    def test_correct_reference_allele_maps_forward(self):
        rec = la.crosswalk_variant("chr1", 5, "A", "G", REF.fetch, n_liftover_mappings=1)
        self.assertEqual(rec["mapping_status"], "mapped")
        self.assertEqual((rec["hg38_ref"], rec["hg38_alt"]), ("A", "G"))
        self.assertEqual(rec["allele_swap"], "none")
        self.assertTrue(rec["is_snv"])

    def test_swapped_alleles_are_flagged_not_silently_repaired(self):
        rec = la.crosswalk_variant("chr1", 5, "G", "A", REF.fetch, n_liftover_mappings=1)
        self.assertEqual(rec["mapping_status"], "mapped")
        self.assertEqual((rec["hg38_ref"], rec["hg38_alt"]), ("A", "G"))
        self.assertEqual(rec["allele_swap"], "swapped")

    def test_incorrect_reference_allele_is_excluded(self):
        rec = la.crosswalk_variant("chr1", 5, "G", "C", REF.fetch, n_liftover_mappings=1)
        self.assertEqual(rec["mapping_status"], "excluded")
        self.assertEqual(rec["exclusion_reason"], "hg38_allele_validation_failed")
        self.assertIsNone(rec["variant_uid"])

    def test_strand_complement_is_recorded(self):
        # position 5 is A; alleles T/C are the reverse complement of A/G
        rec = la.crosswalk_variant("chr1", 5, "T", "C", REF.fetch, n_liftover_mappings=1)
        self.assertEqual(rec["mapping_status"], "mapped")
        self.assertEqual(rec["orientation"], "reverse_complement")
        self.assertEqual((rec["hg38_ref"], rec["hg38_alt"]), ("A", "G"))

    def test_palindromic_alleles_cannot_resolve_strand(self):
        # position 5 is A; A/T is palindromic: forward A/T and rc T/A both match
        rec = la.crosswalk_variant("chr1", 5, "A", "T", REF.fetch, n_liftover_mappings=1)
        self.assertTrue(rec["palindromic"])
        self.assertEqual(rec["orientation"], "forward")

    def test_failed_liftover_is_excluded(self):
        rec = la.crosswalk_variant("chr1", 5, "A", "G", REF.fetch, n_liftover_mappings=0)
        self.assertEqual(rec["exclusion_reason"], "liftover_failed")
        rec2 = la.crosswalk_variant("chr1", 5, "A", "G", REF.fetch, n_liftover_mappings=2)
        self.assertEqual(rec2["exclusion_reason"], "multimapped")

    def test_indel_is_deferred(self):
        rec = la.crosswalk_variant("chr1", 5, "A", "AG", REF.fetch, n_liftover_mappings=1)
        self.assertEqual(rec["exclusion_reason"], "indel_deferred")
        self.assertFalse(rec["is_snv"])

    def test_nonstandard_chromosome_is_excluded(self):
        rec = la.crosswalk_variant("chrUn_KI270302v1", 5, "A", "G", REF.fetch, n_liftover_mappings=1)
        self.assertEqual(rec["exclusion_reason"], "nonstandard_chrom")


class TestBlocks(unittest.TestCase):
    def test_single_linkage_one_megabase(self):
        anchors = [("s1", "1", 100), ("s2", "1", 900_000), ("s3", "1", 2_000_001), ("s4", "2", 5)]
        blocks = la.coarse_blocks(anchors)
        self.assertEqual(blocks["s1"], blocks["s2"])
        self.assertNotEqual(blocks["s2"], blocks["s3"])
        self.assertEqual(blocks["s1"], "chr1:component0001")
        self.assertEqual(blocks["s3"], "chr1:component0002")
        self.assertEqual(blocks["s4"], "chr2:component0001")


class TestQueryFloor(unittest.TestCase):
    def test_floor_and_cumulative_rule_report_excluded_mass(self):
        w = {"v1": 0.5, "v2": 0.3, "v3": 0.19, "v4": 0.00995, "v5": 0.00005}
        sel, excluded = la.select_query_variants(w, floor=1e-4, cumulative=0.999, cap=2000)
        self.assertEqual(sel, {"v1", "v2", "v3", "v4"})
        self.assertAlmostEqual(excluded, 0.00005)

    def test_cap_limits_and_reports_mass(self):
        w = {f"v{i}": 1 / 10 for i in range(10)}
        sel, excluded = la.select_query_variants(w, floor=1e-4, cumulative=0.999, cap=3)
        self.assertEqual(len(sel), 3)
        self.assertAlmostEqual(excluded, 0.7)


class TestMassAssertions(unittest.TestCase):
    def test_duplicate_join_inflates_mass_and_raises(self):
        rows = [("sig1", "v1", 0.6), ("sig1", "v2", 0.4), ("sig1", "v2", 0.4)]
        with self.assertRaises(la.ContractError):
            la.assert_signal_mass(rows, tolerance=1e-6)

    def test_valid_mass_passes(self):
        rows = [("sig1", "v1", 0.6), ("sig1", "v2", 0.4)]
        la.assert_signal_mass(rows, tolerance=1e-6)


class TestPosteriorSummary(unittest.TestCase):
    def test_signed_summary(self):
        w = {"v1": 0.5, "v2": 0.3, "v3": 0.2}
        scores = {"v1": 1.0, "v2": -2.0}  # v3 uncovered
        s = la.posterior_summary(w, scores, is_signed=True)
        self.assertAlmostEqual(s["coverage"], 0.8)
        self.assertAlmostEqual(s["signed"], 0.5 * 1.0 + 0.3 * -2.0)
        self.assertAlmostEqual(s["magnitude"], 0.5 * 1.0 + 0.3 * 2.0)
        self.assertAlmostEqual(s["missing_mass"], 0.2)
        self.assertFalse(s["complete_0.95"])
        self.assertTrue(s["conditional"])

    def test_unsigned_scorer_has_no_signed_summary(self):
        w = {"v1": 1.0}
        s = la.posterior_summary(w, {"v1": 3.0}, is_signed=False)
        self.assertTrue(math.isnan(s["signed"]))
        self.assertAlmostEqual(s["magnitude"], 3.0)

    def test_directional_use_of_unsigned_scorer_raises(self):
        with self.assertRaises(la.ContractError):
            la.require_signed("splice_max_abs", is_signed=False)


class TestTrackClass(unittest.TestCase):
    def test_liver_ontology_is_primary_liver(self):
        self.assertEqual(la.track_class("UBERON:0002107", "liver", "tissue"), "primary_liver")
        self.assertEqual(la.track_class("CL:0000182", "hepatocyte", "primary cell"), "hepatocyte")
        self.assertEqual(la.track_class("EFO:0001187", "HepG2", "cell line"), "HepG2")
        self.assertEqual(la.track_class("UBERON:0002048", "lung", "tissue"), "other")


if __name__ == "__main__":
    unittest.main()


class TestGateSet(unittest.TestCase):
    """The 211-variant AlphaGenome gate set: hg38 alleles come from the 65-series score file, never from the hg19 truth columns."""

    def _fetch(self, chrom, start, end):
        return {("chr10", 31449111): "A", ("chr1", 99): "G"}.get((chrom, start), "N")

    def test_swapped_row_uses_hg38_alleles_and_orients_to_effect_allele(self):
        truth = [{"variant_id": "10:31738041:G:A", "chr": "10", "pos_hg38": "31449112", "ref": "G", "alt": "A", "effect_allele": "A",
                  "ensembl": "ENSG00000237036", "eqtl_sign": "-1", "is_signal_lead": "TRUE", "strand_ambiguous": "FALSE"}]
        ag = [{"variant_id": "10:31738041:G:A", "ensembl": "ENSG00000237036.5", "hg38_ref": "A", "hg38_alt": "G", "is_indel": "FALSE"}]
        rows = la.gate_set(truth, ag, self._fetch)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["variant_uid"], "chr10:31449112:A:G")
        self.assertEqual(rows[0]["orientation_sign"], -1)      # effect allele A == hg38_ref -> negate the ref->alt score
        self.assertEqual(rows[0]["label"], 0)

    def test_fasta_mismatch_raises(self):
        truth = [{"variant_id": "1:100:C:T", "chr": "1", "pos_hg38": "100", "ref": "C", "alt": "T", "effect_allele": "T",
                  "ensembl": "ENSG1", "eqtl_sign": "1", "is_signal_lead": "TRUE", "strand_ambiguous": "FALSE"}]
        ag = [{"variant_id": "1:100:C:T", "ensembl": "ENSG1", "hg38_ref": "C", "hg38_alt": "T", "is_indel": "FALSE"}]
        with self.assertRaises(la.ContractError):
            la.gate_set(truth, ag, self._fetch)                  # FASTA says G at chr1:100

    def test_indel_excluded_and_unmatched_or_bad_effect_allele_raise(self):
        base = {"chr": "1", "pos_hg38": "100", "ref": "G", "alt": "T", "effect_allele": "T", "ensembl": "ENSG1", "eqtl_sign": "1",
                "is_signal_lead": "TRUE", "strand_ambiguous": "FALSE"}
        truth = [{"variant_id": "1:100:G:GT", **base}]
        ag = [{"variant_id": "1:100:G:GT", "ensembl": "ENSG1", "hg38_ref": "G", "hg38_alt": "GT", "is_indel": "TRUE"}]
        self.assertEqual(la.gate_set(truth, ag, self._fetch), [])
        with self.assertRaises(la.ContractError):
            la.gate_set([{"variant_id": "1:100:G:T", **base}], [], self._fetch)   # no hg38 row for the lead
        with self.assertRaises(la.ContractError):
            la.gate_set([{"variant_id": "1:100:G:T", **base, "effect_allele": "C"}],
                        [{"variant_id": "1:100:G:T", "ensembl": "ENSG1", "hg38_ref": "G", "hg38_alt": "T", "is_indel": "FALSE"}], self._fetch)


class TestResultRootOverride(unittest.TestCase):
    """Step scripts read Track 0 deposits by absolute path; the root must be switchable and recorded.

    The allelic arm scored against the PREVIEW (direct-only) liver summaries because the full run had not
    finished. Re-running it against the full archive has to be a stated switch, not an edit.
    """

    def test_default_is_the_named_run_root(self):
        import lib_atlas as la
        os.environ.pop("AGA_TRACK0_ROOT", None)
        self.assertTrue(str(la.track0_root()).endswith("run-20260909T153939Z"))

    def test_env_override_wins(self):
        import lib_atlas as la
        os.environ["AGA_TRACK0_ROOT"] = "/tmp/some-other-root"
        try:
            self.assertEqual(str(la.track0_root()), "/tmp/some-other-root")
        finally:
            os.environ.pop("AGA_TRACK0_ROOT", None)



class TestPosteriorKey(unittest.TestCase):
    """A variant the Atlas defers has no hg38 variant_uid, so a posterior keyed on that column alone
    collapses every deferred indel of a signal into one entry and sums their weights."""

    def test_a_mapped_variant_keys_on_its_hg38_uid(self):
        self.assertEqual(la.posterior_key("chr22:43989339:C:CT", "hg19:22:44385219:C:CT"), "chr22:43989339:C:CT")

    def test_an_unmapped_variant_keys_on_its_normalised_source_id(self):
        self.assertEqual(la.posterior_key("", "hg19:22:44385219:C:CT"), "hg19:chr22:44385219:C:CT")

    def test_bare_and_chr_prefixed_source_ids_agree(self):
        self.assertEqual(la.posterior_key("", "hg19:15:42886074:G:A"),
                         la.posterior_key("", "hg19:chr15:42886074:G:A"))

    def test_two_unmapped_variants_get_distinct_keys(self):
        self.assertNotEqual(la.posterior_key("", "hg19:19:44878467:G:T"),
                            la.posterior_key("", "hg19:19:44892687:T:C"))

    def test_a_variant_with_neither_identity_is_refused(self):
        with self.assertRaises(la.ContractError):
            la.posterior_key("", "")

    def test_a_malformed_source_id_is_refused_rather_than_passed_through(self):
        with self.assertRaises(la.ContractError):
            la.posterior_key("", "rs738409")


class TestAssertUniqueIdentities(unittest.TestCase):
    def test_distinct_identities_pass(self):
        la.assert_unique_identities([("s1", "a"), ("s1", "b"), ("s2", "a")])

    def test_a_repeated_identity_within_one_signal_is_refused(self):
        with self.assertRaises(la.ContractError):
            la.assert_unique_identities([("s1", "a"), ("s1", "a")])

    def test_a_blank_identity_is_refused_even_when_it_appears_once(self):
        """A blank key is the silent-merge bug in its first stage; it never reaches a second row."""
        with self.assertRaises(la.ContractError):
            la.assert_unique_identities([("s1", "")])


class TestResolveIndelAlleles(unittest.TestCase):
    """The hg38 coordinate of a deferred indel was always known; only its alleles and uid were never
    written. Resolving them is fail-closed: the reference base must actually be there."""

    def setUp(self):
        self.seq = "ACGT" * 64          # 256 bp; 1-based position 9 is 'A', 10 is 'C', 11 is 'G', 12 is 'T'
        self.fetch = lambda chrom, s, e: self.seq[s:e]

    def test_an_insertion_with_allele1_as_reference_resolves(self):
        r = la.resolve_indel_alleles("chr1", 10, "C", "CTTT", self.fetch)
        self.assertEqual((r["hg38_ref"], r["hg38_alt"]), ("C", "CTTT"))
        self.assertEqual(r["allele_swap"], "none")
        self.assertEqual(r["variant_uid"], "chr1:10:C:CTTT")

    def test_an_insertion_with_the_alleles_the_other_way_round_resolves_and_records_the_swap(self):
        r = la.resolve_indel_alleles("chr1", 10, "CTTT", "C", self.fetch)
        self.assertEqual((r["hg38_ref"], r["hg38_alt"]), ("C", "CTTT"))
        self.assertEqual(r["allele_swap"], "swapped")

    def test_a_deletion_resolves_against_the_full_reference_span(self):
        r = la.resolve_indel_alleles("chr1", 10, "CGT", "C", self.fetch)
        self.assertEqual((r["hg38_ref"], r["hg38_alt"]), ("CGT", "C"))
        self.assertEqual(r["variant_uid"], "chr1:10:CGT:C")

    def test_a_reference_that_is_not_there_is_refused_not_repaired(self):
        r = la.resolve_indel_alleles("chr1", 10, "A", "AT", self.fetch)
        self.assertIsNone(r["variant_uid"])
        self.assertEqual(r["exclusion_reason"], "indel_deferred_ref_mismatch")

    def test_the_reverse_complement_is_never_tried_for_an_indel(self):
        """A multi-base allele's strand flip is not a base-wise complement of the VCF-anchored form;
        matching one would invent a variant that is not at that site."""
        r = la.resolve_indel_alleles("chr1", 10, "G", "GAAA", self.fetch)   # complement of C is G
        self.assertIsNone(r["variant_uid"])

    def test_the_resolved_row_stays_excluded_from_the_atlas_query_set(self):
        r = la.resolve_indel_alleles("chr1", 10, "C", "CTTT", self.fetch)
        self.assertEqual(r["mapping_status"], "excluded")
        self.assertEqual(r["exclusion_reason"], "indel_deferred")


class TestIndelOrientationFollowsTheSourceThenTheReference(unittest.TestCase):
    """The reference alone cannot orient a prefix-anchored indel: if the long allele is present then the
    short one is too, so BOTH readings are representable. The source's own allele order is the only other
    information there is, so it is tried first and the reference is the arbiter only when it contradicts it.
    Choosing the longer allele instead would call HSD17B13 rs72613567 a deletion, and it is a known
    T->TA insertion at the exon-6 donor."""

    def setUp(self):
        self.seq = "ACGT" * 64        # 1-based 10 is 'C', 11 'G', 12 'T'
        self.fetch = lambda chrom, s, e: self.seq[s:e]

    def test_the_source_order_is_kept_when_the_reference_permits_it(self):
        r = la.resolve_indel_alleles("chr1", 10, "C", "CGT", self.fetch)
        self.assertEqual((r["hg38_ref"], r["hg38_alt"]), ("C", "CGT"))
        self.assertEqual(r["allele_swap"], "none")

    def test_the_other_order_is_used_when_the_source_order_is_not_at_the_reference(self):
        """This is the fallback that was missing: it recovers all 12 targets step 63 refused."""
        r = la.resolve_indel_alleles("chr1", 11, "GG", "G", self.fetch)   # chr1:11 reads "GT", not "GG"
        self.assertEqual((r["hg38_ref"], r["hg38_alt"]), ("G", "GG"))
        self.assertEqual(r["allele_swap"], "swapped")

    def test_a_prefix_ambiguity_is_flagged_rather_than_silently_decided(self):
        """Both readings match the reference here; the choice rested on the source order, not on evidence."""
        r = la.resolve_indel_alleles("chr1", 10, "C", "CGT", self.fetch)
        self.assertTrue(r["orientation_ambiguous"])

    def test_a_resolved_orientation_is_not_flagged(self):
        r = la.resolve_indel_alleles("chr1", 11, "GG", "G", self.fetch)
        self.assertFalse(r["orientation_ambiguous"])

    def test_an_snv_is_never_ambiguous(self):
        r = la.resolve_indel_alleles("chr1", 10, "T", "C", self.fetch)
        self.assertEqual((r["hg38_ref"], r["hg38_alt"]), ("C", "T"))
        self.assertFalse(r["orientation_ambiguous"])

    def test_neither_allele_present_is_refused(self):
        r = la.resolve_indel_alleles("chr1", 10, "A", "AGT", self.fetch)
        self.assertIsNone(r["variant_uid"])
        self.assertEqual(r["exclusion_reason"], "indel_deferred_ref_mismatch")


class TestCrosswalkVariantIndelBranch(unittest.TestCase):
    """The crosswalk kept indels out of the Atlas query set, which is right, but also left them with no
    identity at all, which is what let a downstream posterior merge them."""

    def setUp(self):
        self.seq = "ACGT" * 64
        self.fetch = lambda chrom, s, e: self.seq[s:e]

    def test_a_deferred_indel_now_carries_a_uid(self):
        r = la.crosswalk_variant("chr1", 10, "C", "CTTT", self.fetch, 1)
        self.assertEqual(r["variant_uid"], "chr1:10:C:CTTT")
        self.assertEqual(r["exclusion_reason"], "indel_deferred")

    def test_a_deferred_indel_is_still_excluded_from_the_query_set(self):
        r = la.crosswalk_variant("chr1", 10, "C", "CTTT", self.fetch, 1)
        self.assertEqual(r["mapping_status"], "excluded")

    def test_an_indel_whose_reference_is_absent_keeps_no_uid_and_says_so(self):
        r = la.crosswalk_variant("chr1", 10, "A", "AT", self.fetch, 1)
        self.assertIsNone(r["variant_uid"])
        self.assertEqual(r["exclusion_reason"], "indel_deferred_ref_mismatch")

    def test_the_snv_path_is_unchanged(self):
        r = la.crosswalk_variant("chr1", 10, "C", "T", self.fetch, 1)
        self.assertEqual(r["mapping_status"], "mapped")
        self.assertEqual(r["variant_uid"], "chr1:10:C:T")

    def test_a_failed_liftover_still_has_no_uid_because_it_has_no_coordinate(self):
        r = la.crosswalk_variant("chr1", 10, "C", "CTTT", self.fetch, 0)
        self.assertIsNone(r["variant_uid"])
        self.assertEqual(r["exclusion_reason"], "liftover_failed")

    def test_a_multimapped_indel_is_not_resolved(self):
        r = la.crosswalk_variant("chr1", 10, "C", "CTTT", self.fetch, 2)
        self.assertIsNone(r["variant_uid"])
        self.assertEqual(r["exclusion_reason"], "multimapped")

"""Fixture tests for step 67: recovering the hg38 identity of deferred indels in an already-written run."""

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


class TestRepairOneRow(unittest.TestCase):
    def setUp(self):
        self.m = load("67_repair_deferred_indel_uids.py")
        self.seq = "ACGT" * 64
        self.fetch = lambda chrom, s, e: self.seq[s:e]

    def row(self, **kw):
        base = {"source_variant_id": "hg19:1:99:C:CTTT", "hg38_chrom": "chr1", "hg38_position_1based": "10",
                "source_alleles": "C/CTTT", "n_liftover_mappings": "1",
                "mapping_status": "excluded", "exclusion_reason": "indel_deferred", "variant_uid": ""}
        base.update(kw)
        return base

    def test_a_deferred_indel_recovers_its_uid(self):
        out = self.m.repair_row(self.row(), self.fetch)
        self.assertEqual(out["recovered_variant_uid"], "chr1:10:C:CTTT")
        self.assertEqual(out["repair_state"], "recovered")

    def test_the_recovered_row_is_still_not_served_by_the_atlas(self):
        """Recovering an identity does not make the point-query API accept an indel."""
        out = self.m.repair_row(self.row(), self.fetch)
        self.assertEqual(out["mapping_status"], "excluded")
        self.assertEqual(out["atlas_queryable"], "False")

    def test_a_reference_mismatch_is_recorded_not_repaired(self):
        out = self.m.repair_row(self.row(source_alleles="A/AT"), self.fetch)
        self.assertEqual(out["recovered_variant_uid"], "")
        self.assertEqual(out["repair_state"], "ref_mismatch")

    def test_a_row_that_already_has_a_uid_is_left_alone(self):
        out = self.m.repair_row(self.row(variant_uid="chr1:10:C:T", mapping_status="mapped",
                                         exclusion_reason=""), self.fetch)
        self.assertEqual(out["repair_state"], "already_identified")
        self.assertEqual(out["recovered_variant_uid"], "chr1:10:C:T")

    def test_a_row_with_no_coordinate_cannot_be_repaired(self):
        out = self.m.repair_row(self.row(hg38_chrom="", hg38_position_1based="",
                                         exclusion_reason="liftover_failed", n_liftover_mappings="0"), self.fetch)
        self.assertEqual(out["repair_state"], "no_hg38_coordinate")

    def test_a_multimapped_row_is_not_repaired_even_with_a_coordinate(self):
        out = self.m.repair_row(self.row(exclusion_reason="multimapped", n_liftover_mappings="2"), self.fetch)
        self.assertEqual(out["repair_state"], "not_one_to_one")


class TestAlleleParsing(unittest.TestCase):
    def setUp(self):
        self.m = load("67_repair_deferred_indel_uids.py")

    def test_source_alleles_split_on_the_slash(self):
        self.assertEqual(self.m.source_alleles("C/CTTT"), ("C", "CTTT"))

    def test_a_field_that_is_not_a_pair_is_refused(self):
        with self.assertRaises(Exception):
            self.m.source_alleles("C")


class TestIdentityRecovery(unittest.TestCase):
    """The whole point: after the repair, no two deferred indels of a signal share a key."""

    def setUp(self):
        self.m = load("67_repair_deferred_indel_uids.py")

    def test_repaired_keys_separate_two_indels_that_shared_the_empty_uid(self):
        a = self.m.identity_after_repair({"variant_uid": "", "recovered_variant_uid": "chr1:10:C:CT",
                                          "source_variant_id": "hg19:1:99:C:CT"})
        b = self.m.identity_after_repair({"variant_uid": "", "recovered_variant_uid": "chr1:20:G:GA",
                                          "source_variant_id": "hg19:1:109:G:GA"})
        self.assertNotEqual(a, b)

    def test_an_unrecovered_row_still_falls_back_to_its_source_id(self):
        k = self.m.identity_after_repair({"variant_uid": "", "recovered_variant_uid": "",
                                          "source_variant_id": "hg19:1:99:C:CT"})
        self.assertEqual(k, "hg19:chr1:99:C:CT")


class TestCollisionGuard(unittest.TestCase):
    """The crosswalk is keyed on the SOURCE string, and universe A spells a chromosome `10` where universe B
    spells it `chr10`, so one variant reached through both universes has two rows. Those rows must agree on
    the hg38 site; two rows that disagree would be a real defect."""

    def setUp(self):
        self.m = load("67_repair_deferred_indel_uids.py")

    def rows(self, chrom_b="chr10", pos_b="112173248"):
        return [{"identity_after_repair": "hg19:chr10:113933006:ATT:A", "source_variant_id": "hg19:10:113933006:ATT:A",
                 "hg38_chrom": "chr10", "hg38_position_1based": "112173248", "hg38_ref": "ATT", "hg38_alt": "A"},
                {"identity_after_repair": "hg19:chr10:113933006:ATT:A", "source_variant_id": "hg19:chr10:113933006:ATT:A",
                 "hg38_chrom": chrom_b, "hg38_position_1based": pos_b, "hg38_ref": "ATT", "hg38_alt": "A"}]

    def test_two_spellings_of_the_same_variant_are_allowed(self):
        self.assertEqual(self.m.assert_identities_agree(self.rows()), 1)

    def test_two_rows_that_disagree_on_the_site_are_refused(self):
        with self.assertRaises(Exception):
            self.m.assert_identities_agree(self.rows(pos_b="999"))

    def test_two_rows_that_disagree_on_the_chromosome_are_refused(self):
        with self.assertRaises(Exception):
            self.m.assert_identities_agree(self.rows(chrom_b="chr11"))

    def test_distinct_identities_report_no_collisions(self):
        r = self.rows()
        r[1]["identity_after_repair"] = "hg19:chr10:113938330:AAAAC:A"
        self.assertEqual(self.m.assert_identities_agree(r), 0)


class TestSignalIdentities(unittest.TestCase):
    """The guard that matters. Checking that rows sharing a recovered uid agree on the hg38 site is VACUOUS:
    the uid encodes that site. The invariant worth testing is that keying a posterior on the recovered uid
    does not merge two distinct source variants within one signal, which is the bug being repaired."""

    def setUp(self):
        self.m = load("67_repair_deferred_indel_uids.py")
        self.repair = {"hg19:1:99:T:TG": "chr1:50:TG:T", "hg19:1:99:TG:T": "chr1:50:TG:T",
                       "hg19:1:80:C:CA": "chr1:40:C:CA"}

    def test_two_spellings_of_one_variant_in_different_signals_are_fine(self):
        rows = [{"signal_uid": "sA", "variant_uid": "", "source_variant_id": "hg19:1:99:T:TG"},
                {"signal_uid": "sB", "variant_uid": "", "source_variant_id": "hg19:1:99:TG:T"}]
        self.assertEqual(len(self.m.signal_identities(rows, self.repair)), 2)

    def test_two_distinct_variants_in_one_signal_keep_distinct_identities(self):
        rows = [{"signal_uid": "sA", "variant_uid": "", "source_variant_id": "hg19:1:99:T:TG"},
                {"signal_uid": "sA", "variant_uid": "", "source_variant_id": "hg19:1:80:C:CA"}]
        pairs = self.m.signal_identities(rows, self.repair)
        self.assertEqual(len({k for _, k in pairs}), 2)

    def test_a_repair_that_would_merge_two_entries_of_one_signal_is_refused(self):
        rows = [{"signal_uid": "sA", "variant_uid": "", "source_variant_id": "hg19:1:99:T:TG"},
                {"signal_uid": "sA", "variant_uid": "", "source_variant_id": "hg19:1:99:TG:T"}]
        with self.assertRaises(Exception):
            self.m.assert_no_merge(self.m.signal_identities(rows, self.repair))

    def test_a_mapped_row_keeps_its_own_uid_not_a_repaired_one(self):
        rows = [{"signal_uid": "sA", "variant_uid": "chr1:9:C:T", "source_variant_id": "hg19:1:99:T:TG"}]
        self.assertEqual(self.m.signal_identities(rows, self.repair), [("sA", "chr1:9:C:T")])

    def test_an_unrepaired_row_falls_back_to_its_source_id(self):
        rows = [{"signal_uid": "sA", "variant_uid": "", "source_variant_id": "hg19:1:70:G:GG"}]
        self.assertEqual(self.m.signal_identities(rows, self.repair), [("sA", "hg19:chr1:70:G:GG")])


class TestOrientationAmbiguityIsCarried(unittest.TestCase):
    """The recovered uid is only as good as the orientation behind it, and for most deferred indels the
    reference permits both readings. The repair table has to say which rows those are."""

    def setUp(self):
        self.m = load("67_repair_deferred_indel_uids.py")
        self.seq = "ACGT" * 64
        self.fetch = lambda chrom, s, e: self.seq[s:e]

    def row(self, alleles, pos="10"):
        return {"source_variant_id": "hg19:1:99:x", "hg38_chrom": "chr1", "hg38_position_1based": pos,
                "source_alleles": alleles, "n_liftover_mappings": "1", "mapping_status": "excluded",
                "exclusion_reason": "indel_deferred", "variant_uid": ""}

    def test_an_ambiguous_indel_is_flagged(self):
        out = self.m.repair_row(self.row("C/CGT"), self.fetch)
        self.assertEqual(out["orientation_ambiguous"], "True")

    def test_a_reference_resolved_indel_is_not_flagged(self):
        """chr1:11 reads "GT", so only the short allele is there and the reference does the orienting."""
        out = self.m.repair_row(self.row("GG/G", pos="11"), self.fetch)
        self.assertEqual(out["repair_state"], "recovered")
        self.assertEqual(out["orientation_ambiguous"], "False")

    def test_a_reference_mismatch_carries_no_orientation_claim_at_all(self):
        out = self.m.repair_row(self.row("GG/G"), self.fetch)
        self.assertEqual(out["repair_state"], "ref_mismatch")
        self.assertEqual(out["orientation_ambiguous"], "")

    def test_the_flag_is_a_string_so_it_survives_the_tsv_round_trip(self):
        out = self.m.repair_row(self.row("C/CGT"), self.fetch)
        self.assertIsInstance(out["orientation_ambiguous"], str)


if __name__ == "__main__":
    unittest.main()

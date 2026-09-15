"""Fixture tests for step 65 (P6d: the per-variant clinical dossier and the coverage gate it carries)."""

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


class TestSourceIdNormalisation(unittest.TestCase):
    """The weights table writes hg19:15:..., the indel rescue writes hg19:chr22:...; a join on the raw
    string silently drops every indel, which is exactly the class this lane exists to surface."""

    def setUp(self):
        self.m = load("65_clinical_dossier.py")

    def test_bare_chromosome_and_chr_prefixed_forms_agree(self):
        self.assertEqual(self.m.normalise_source_id("hg19:15:42886074:G:A"),
                         self.m.normalise_source_id("hg19:chr15:42886074:G:A"))

    def test_alleles_are_upper_cased_but_not_reordered(self):
        self.assertEqual(self.m.normalise_source_id("hg19:chr22:44385219:c:ct"), "hg19:chr22:44385219:C:CT")

    def test_a_malformed_id_is_refused_rather_than_passed_through(self):
        with self.assertRaises(Exception):
            self.m.normalise_source_id("rs738409")


class TestVariantClass(unittest.TestCase):
    def setUp(self):
        self.m = load("65_clinical_dossier.py")

    def test_single_base_swap_is_an_snv(self):
        self.assertEqual(self.m.variant_class("C", "T"), "snv")

    def test_longer_alt_is_an_insertion(self):
        self.assertEqual(self.m.variant_class("C", "CT"), "insertion")

    def test_longer_ref_is_a_deletion(self):
        self.assertEqual(self.m.variant_class("CGT", "C"), "deletion")

    def test_equal_length_multi_base_is_an_mnv_not_an_snv(self):
        self.assertEqual(self.m.variant_class("CA", "TG"), "mnv")


class TestRankAndCumulative(unittest.TestCase):
    """Rank and cumulative share are properties of the FULL posterior, not of the exported subset."""

    def setUp(self):
        self.m = load("65_clinical_dossier.py")

    def test_rank_is_one_based_and_descending(self):
        r = self.m.rank_and_cumulative({"a": 0.1, "b": 0.6, "c": 0.3})
        self.assertEqual(r["b"]["rank"], 1)
        self.assertEqual(r["c"]["rank"], 2)
        self.assertEqual(r["a"]["rank"], 3)

    def test_cumulative_share_is_over_the_full_vector(self):
        r = self.m.rank_and_cumulative({"a": 0.1, "b": 0.6, "c": 0.3})
        self.assertAlmostEqual(r["b"]["cumulative_share"], 0.6, places=12)
        self.assertAlmostEqual(r["a"]["cumulative_share"], 1.0, places=12)

    def test_a_partial_posterior_is_not_renormalised(self):
        """The weights of a signal need not sum to one; printing a renormalised share would
        overstate how much of the signal a variant carries."""
        r = self.m.rank_and_cumulative({"a": 0.2, "b": 0.3})
        self.assertAlmostEqual(r["b"]["cumulative_share"], 0.3, places=12)
        self.assertAlmostEqual(r["a"]["cumulative_share"], 0.5, places=12)

    def test_an_empty_posterior_is_refused(self):
        with self.assertRaises(Exception):
            self.m.rank_and_cumulative({})


class TestSelection(unittest.TestCase):
    def setUp(self):
        self.m = load("65_clinical_dossier.py")

    def test_variants_above_the_floor_are_selected(self):
        sel = self.m.select_dossier_variants({"a": 0.5, "b": 0.02, "c": 0.001}, floor=0.01)
        self.assertEqual(sel, {"a", "b"})

    def test_the_top_variant_is_always_selected_even_below_the_floor(self):
        """A signal whose whole posterior sits on deferred indels would otherwise vanish from the export."""
        sel = self.m.select_dossier_variants({"a": 0.004, "b": 0.001}, floor=0.01)
        self.assertEqual(sel, {"a"})

    def test_ties_at_the_top_select_one_deterministically(self):
        sel = self.m.select_dossier_variants({"b": 0.004, "a": 0.004}, floor=0.01)
        self.assertEqual(sel, {"a"})


class TestReadability(unittest.TestCase):
    def setUp(self):
        self.m = load("65_clinical_dossier.py")

    def test_a_covered_signal_with_a_served_variant_is_readable(self):
        self.assertEqual(self.m.row_readable("covered", True), ("True", ""))

    def test_an_unserved_variant_is_not_readable_and_says_why(self):
        ok, why = self.m.row_readable("covered", False)
        self.assertEqual(ok, "False")
        self.assertIn("not served", why)

    def test_a_served_variant_in_a_partial_signal_is_not_readable_alone(self):
        """The variant has a prediction, but the signal's summary columns rest on incomplete mass."""
        ok, why = self.m.row_readable("partial_indel", True)
        self.assertEqual(ok, "False")
        self.assertIn("partial_indel", why)

    def test_an_uncovered_signal_is_not_readable(self):
        ok, why = self.m.row_readable("uncovered_indel", False)
        self.assertEqual(ok, "False")


class TestUnservedRowGuard(unittest.TestCase):
    def setUp(self):
        self.m = load("65_clinical_dossier.py")

    def test_an_unserved_row_carrying_a_numeric_quantile_is_refused(self):
        with self.assertRaises(Exception):
            self.m.guard_unserved_row({"atlas_served": "False", "atac_liver_quantile": "0.0",
                                       "dnase_liver_quantile": ""})

    def test_an_unserved_row_with_empty_channels_passes(self):
        self.m.guard_unserved_row({"atlas_served": "False", "atac_liver_quantile": "",
                                   "dnase_liver_quantile": ""})

    def test_a_served_row_with_a_zero_quantile_is_fine(self):
        self.m.guard_unserved_row({"atlas_served": "True", "atac_liver_quantile": "0.0"})

    def test_a_model_api_column_on_an_unserved_row_is_allowed(self):
        """The rescue lane exists precisely to fill these on rows the Atlas refused."""
        self.m.guard_unserved_row({"atlas_served": "False", "atac_liver_quantile": "",
                                   "model_api_splice_log2": "-1.09"})


class TestChannelColumnSeparation(unittest.TestCase):
    def setUp(self):
        self.m = load("65_clinical_dossier.py")

    def test_atlas_and_model_api_column_names_are_disjoint(self):
        self.assertEqual(set(self.m.ATLAS_QUANTILE_COLUMNS) & set(self.m.MODEL_API_COLUMNS), set())

    def test_every_model_api_column_is_named_for_its_service_and_scale(self):
        for c in self.m.MODEL_API_COLUMNS:
            self.assertTrue(c.startswith("model_api_"), c)
            self.assertTrue(c.endswith("_log2"), c)


class TestDirectionConcordance(unittest.TestCase):
    def setUp(self):
        self.m = load("65_clinical_dossier.py")

    def test_a_thin_family_is_withheld_rather_than_printed(self):
        rows = [{"direction_concordant": "True"} for _ in range(5)]
        out = self.m.direction_concordance(rows, min_rows=30)
        self.assertIsNone(out["fraction"])
        self.assertEqual(out["n"], 5)

    def test_a_full_family_reports_the_fraction(self):
        rows = [{"direction_concordant": "True"} for _ in range(20)] + \
               [{"direction_concordant": "False"} for _ in range(20)]
        out = self.m.direction_concordance(rows, min_rows=30)
        self.assertAlmostEqual(out["fraction"], 0.5, places=12)
        self.assertEqual(out["n"], 40)

    def test_rows_with_no_direction_do_not_enter_the_denominator(self):
        rows = [{"direction_concordant": "True"}] * 30 + [{"direction_concordant": ""}] * 100
        out = self.m.direction_concordance(rows, min_rows=30)
        self.assertEqual(out["n"], 30)
        self.assertAlmostEqual(out["fraction"], 1.0, places=12)


class TestCoveredDenominator(unittest.TestCase):
    def setUp(self):
        self.m = load("65_clinical_dossier.py")

    def test_the_share_is_reported_over_the_covered_subset(self):
        out = self.m.measured_peak_share(n_in_peak=3, n_with_peak_data=10, n_rows=100)
        self.assertAlmostEqual(out["share_over_covered"], 0.3, places=12)
        self.assertEqual(out["denominator"], 10)

    def test_thin_peak_coverage_is_flagged_not_extrapolated(self):
        out = self.m.measured_peak_share(n_in_peak=3, n_with_peak_data=10, n_rows=100)
        self.assertTrue(out["coverage_too_thin_to_extrapolate"])

    def test_broad_peak_coverage_is_not_flagged(self):
        out = self.m.measured_peak_share(n_in_peak=30, n_with_peak_data=80, n_rows=100)
        self.assertFalse(out["coverage_too_thin_to_extrapolate"])


class TestRowKey(unittest.TestCase):
    """A variant the Atlas refused has NO hg38 variant_uid: the crosswalk leaves it empty and only the hg19
    source id identifies it. Keying the posterior on variant_uid alone collapses every deferred indel of a
    signal into one entry, silently summing their weights and erasing the class this lane exists to show."""

    def setUp(self):
        self.m = load("65_clinical_dossier.py")

    def test_a_mapped_variant_keys_on_its_hg38_uid(self):
        self.assertEqual(self.m.row_key("chr22:43989339:C:CT", "hg19:22:44385219:C:CT"), "chr22:43989339:C:CT")

    def test_an_excluded_variant_keys_on_its_normalised_source_id(self):
        self.assertEqual(self.m.row_key("", "hg19:22:44385219:C:CT"), "hg19:chr22:44385219:C:CT")

    def test_two_excluded_variants_in_one_signal_get_distinct_keys(self):
        a = self.m.row_key("", "hg19:19:44878467:G:T")
        b = self.m.row_key("", "hg19:19:44892687:T:C")
        self.assertNotEqual(a, b)

    def test_a_variant_with_neither_identity_is_refused(self):
        with self.assertRaises(Exception):
            self.m.row_key("", "")


class TestAllelesForClassing(unittest.TestCase):
    def setUp(self):
        self.m = load("65_clinical_dossier.py")

    def test_alleles_come_from_the_hg38_uid_when_it_exists(self):
        self.assertEqual(self.m.alleles_of("chr22:43989339:C:CT", "hg19:22:44385219:G:GA"), ("C", "CT"))

    def test_alleles_fall_back_to_the_source_id_for_an_excluded_variant(self):
        """Without this an unmapped indel would be classed by splitting an empty string."""
        self.assertEqual(self.m.alleles_of("", "hg19:22:44385219:C:CT"), ("C", "CT"))

    def test_the_deferred_indel_is_classed_as_an_indel_not_an_snv(self):
        ref, alt = self.m.alleles_of("", "hg19:22:44385219:C:CT")
        self.assertEqual(self.m.variant_class(ref, alt), "insertion")


class TestRecoveredAlleles(unittest.TestCase):
    """`alleles_of` falls back to the SOURCE id for a variant with no hg38 uid, and a GWAS source id carries
    effect/other alleles in no particular order. That mislabels a deletion as an insertion. Step 67 recovers
    the reference orientation, and the dossier should prefer it when it exists."""

    def setUp(self):
        self.m = load("65_clinical_dossier.py")

    def test_the_recovered_uid_supplies_the_alleles_when_the_source_order_is_wrong(self):
        ref, alt = self.m.alleles_of("", "hg19:1:99:C:CGT", recovered_uid="chr1:50:CGT:C")
        self.assertEqual((ref, alt), ("CGT", "C"))
        self.assertEqual(self.m.variant_class(ref, alt), "deletion")

    def test_without_a_recovered_uid_it_still_falls_back_to_the_source(self):
        self.assertEqual(self.m.alleles_of("", "hg19:1:99:C:CGT"), ("C", "CGT"))

    def test_the_atlas_uid_still_wins_over_a_recovered_one(self):
        self.assertEqual(self.m.alleles_of("chr1:9:C:T", "hg19:1:99:C:CGT", recovered_uid="chr1:50:CGT:C"),
                         ("C", "T"))


class TestOrientationPrecedenceInTheDossier(unittest.TestCase):
    """Three sources can name a deferred indel's reference allele, and they do not agree: dbSNP (build 157,
    canonical REF/ALT), the reference-only repair (step 67), and the raw source id. dbSNP flips the
    reference-only call on 61% of the cases the reference cannot decide, so it has to win."""

    def setUp(self):
        self.m = load("65_clinical_dossier.py")

    def test_dbsnp_outranks_the_reference_only_repair(self):
        ref, alt = self.m.alleles_of("", "hg19:1:99:C:CGT",
                                     recovered_uid="chr1:50:C:CGT", dbsnp_uid="chr1:50:CGT:C")
        self.assertEqual((ref, alt), ("CGT", "C"))

    def test_the_repair_is_used_when_dbsnp_is_silent(self):
        ref, alt = self.m.alleles_of("", "hg19:1:99:C:CGT", recovered_uid="chr1:50:C:CGT", dbsnp_uid="")
        self.assertEqual((ref, alt), ("C", "CGT"))

    def test_the_source_id_is_the_last_resort(self):
        self.assertEqual(self.m.alleles_of("", "hg19:1:99:C:CGT"), ("C", "CGT"))

    def test_an_atlas_served_uid_still_outranks_everything(self):
        """A served variant has a real Atlas identity; no indel table may override it."""
        self.assertEqual(self.m.alleles_of("chr1:9:C:T", "hg19:1:99:C:CGT",
                                           recovered_uid="chr1:50:C:CGT", dbsnp_uid="chr1:50:CGT:C"),
                         ("C", "T"))

    def test_orientation_source_is_recorded_so_a_reader_can_stratify(self):
        self.assertEqual(self.m.orientation_source("chr1:50:CGT:C", "chr1:50:C:CGT"), "dbsnp")
        self.assertEqual(self.m.orientation_source("", "chr1:50:C:CGT"), "reference_or_source_order")
        self.assertEqual(self.m.orientation_source("", ""), "")


class TestRescueJoin(unittest.TestCase):
    """Step 63 scores one row per distinct variant under its heaviest source id and lists every id of that
    variant. A dossier row can carry any of those ids, and it may only borrow an effect scored on the same
    alleles it reports."""

    def setUp(self):
        self.m = load("65_clinical_dossier.py")
        self.rows = [{"arm": "indel", "state": "scored", "source_variant_id": "hg19:10:5:TCA:T",
                      "member_source_ids": "hg19:10:5:TCA:T;hg19:chr10:5:T:TCA", "ref": "T", "alt": "TCA",
                      "splice_log2": "0.5"},
                     {"arm": "snv_control", "state": "scored", "source_variant_id": "hg19:10:9:A:G",
                      "ref": "A", "alt": "G"}]

    def test_every_member_id_finds_the_scored_row(self):
        idx = self.m.rescue_index(self.rows)
        self.assertIs(idx[self.m.normalise_source_id("hg19:10:5:T:TCA")], self.rows[0])
        self.assertIs(idx[self.m.normalise_source_id("hg19:10:5:TCA:T")], self.rows[0])

    def test_controls_are_not_indexed(self):
        self.assertNotIn(self.m.normalise_source_id("hg19:10:9:A:G"), self.m.rescue_index(self.rows))

    def test_a_row_from_before_grouping_is_indexed_by_its_own_id(self):
        old = [{"arm": "indel", "state": "scored", "source_variant_id": "hg19:1:5:C:CT", "ref": "C", "alt": "CT"}]
        self.assertIn(self.m.normalise_source_id("hg19:1:5:C:CT"), self.m.rescue_index(old))

    def test_an_effect_scored_on_the_same_alleles_is_attached(self):
        idx = self.m.rescue_index(self.rows)
        row, state = self.m.rescue_for_row(idx, "hg19:chr10:5:T:TCA", "T", "TCA")
        self.assertEqual((row["splice_log2"], state), ("0.5", "scored"))

    def test_an_effect_scored_on_the_reciprocal_variant_is_withheld_and_flagged(self):
        idx = self.m.rescue_index(self.rows)
        row, state = self.m.rescue_for_row(idx, "hg19:10:5:TCA:T", "TCA", "T")
        self.assertEqual((row, state), ({}, "orientation_differs"))

    def test_no_rescue_gives_an_empty_state(self):
        self.assertEqual(self.m.rescue_for_row({}, "hg19:10:5:TCA:T", "TCA", "T"), ({}, ""))


if __name__ == "__main__":
    unittest.main()

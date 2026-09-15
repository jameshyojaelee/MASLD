"""Fixture tests for step 64 (the indel-vs-SNV paired comparison and its orientation strata)."""

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


class TestBenjaminiHochberg(unittest.TestCase):
    def setUp(self):
        self.m = load("64_indel_vs_snv_paired.py")

    def test_a_single_test_is_its_own_q(self):
        self.assertAlmostEqual(self.m.benjamini_hochberg([0.03])[0], 0.03, places=12)

    def test_the_largest_p_is_never_shrunk(self):
        q = self.m.benjamini_hochberg([0.01, 0.5])
        self.assertAlmostEqual(q[1], 0.5, places=12)

    def test_q_values_are_monotone_in_p(self):
        q = self.m.benjamini_hochberg([0.001, 0.01, 0.04, 0.2, 0.5])
        self.assertEqual(q, sorted(q))

    def test_the_step_up_enforcement_pulls_an_earlier_q_down(self):
        """Without the running minimum a small-rank test can end up with a larger q than a later one."""
        q = self.m.benjamini_hochberg([0.04, 0.041])
        self.assertLessEqual(q[0], q[1])

    def test_a_nan_p_does_not_silently_become_significant(self):
        q = self.m.benjamini_hochberg([float("nan"), 0.01])
        self.assertTrue(q[0] != q[0])


class TestOrientationStrata(unittest.TestCase):
    """87.7% of this Resource's deferred indels are orientation-ambiguous: the reference permits both
    readings and the choice rests on the source's allele order. The comparison must therefore be reported
    both over all targets and over the minority the reference itself orients."""

    def setUp(self):
        self.m = load("64_indel_vs_snv_paired.py")

    def rows(self):
        return [{"arm": "indel", "state": "scored", "orientation_ambiguous": "False",
                 "orientation_source": "reference_or_source_order", "tag": "a"},
                {"arm": "indel", "state": "scored", "orientation_ambiguous": "True",
                 "orientation_source": "dbsnp", "tag": "b"},
                {"arm": "indel", "state": "scored", "orientation_ambiguous": "",
                 "orientation_source": "", "tag": "c"}]

    def test_the_dbsnp_stratum_keeps_only_dbsnp_oriented_indels(self):
        """dbSNP outranks the reference rule, so a dbSNP-oriented row belongs to that stratum even when the
        reference alone could not have oriented it -- which is the majority of the useful cases."""
        self.assertEqual([r["tag"] for r in self.m.stratum(self.rows(), "dbsnp_oriented")], ["b"])

    def test_the_all_stratum_keeps_every_scored_indel(self):
        self.assertEqual(len(self.m.stratum(self.rows(), "all")), 3)

    def test_the_resolved_stratum_keeps_only_unambiguous_indels(self):
        keep = self.m.stratum(self.rows(), "reference_resolved")
        self.assertEqual([r["tag"] for r in keep], ["a"])

    def test_a_blank_flag_is_not_treated_as_resolved(self):
        """A row from a run that predates the flag must not be counted as evidence either way."""
        self.assertNotIn("c", [r["tag"] for r in self.m.stratum(self.rows(), "reference_resolved")])

    def test_reference_resolved_ids_come_from_the_repair_table(self):
        ids = self.m.reference_resolved_ids([
            {"source_variant_id": "a", "repair_state": "recovered", "orientation_ambiguous": "False"},
            {"source_variant_id": "b", "repair_state": "recovered", "orientation_ambiguous": "True"},
            {"source_variant_id": "c", "repair_state": "ref_mismatch", "orientation_ambiguous": ""}])
        self.assertEqual(ids, {"a"})

    def test_an_unknown_stratum_is_refused(self):
        with self.assertRaises(Exception):
            self.m.stratum(self.rows(), "whatever")

    def test_the_same_variant_scored_twice_is_refused(self):
        """Two rows for one variant are one observation; pairing both would count it twice."""
        rows = [{"arm": "indel", "state": "scored", "chrom": "chr1", "pos_hg38": "7", "ref": "T", "alt": "TCA"},
                {"arm": "indel", "state": "scored", "chrom": "chr1", "pos_hg38": "7", "ref": "TCA", "alt": "T"}]
        with self.assertRaises(self.m.la.ContractError):
            self.m.assert_distinct_variants(rows)

    def test_distinct_variants_pass(self):
        rows = [{"arm": "indel", "state": "scored", "chrom": "chr1", "pos_hg38": "7", "ref": "T", "alt": "TCA"},
                {"arm": "indel", "state": "scored", "chrom": "chr1", "pos_hg38": "7", "ref": "T", "alt": "TC"}]
        self.m.assert_distinct_variants(rows)

    def test_the_unique_record_stratum_drops_frequency_tie_breaks(self):
        """The frequency tie-break is a second, weaker rule; a stratum without it shows whether the result
        depends on that rule."""
        rows = [{"arm": "indel", "state": "scored", "orientation_source": "dbsnp",
                 "dbsnp_rule": "unique_record", "tag": "u"},
                {"arm": "indel", "state": "scored", "orientation_source": "dbsnp",
                 "dbsnp_rule": "frequency", "tag": "f"},
                {"arm": "indel", "state": "scored", "orientation_source": "dbsnp",
                 "dbsnp_rule": "", "tag": "blank"},
                {"arm": "indel", "state": "scored", "orientation_source": "reference_or_source_order",
                 "dbsnp_rule": "", "tag": "r"}]
        self.assertEqual([r["tag"] for r in self.m.stratum(rows, "dbsnp_unique_record")], ["u"])
        self.assertEqual([r["tag"] for r in self.m.stratum(rows, "dbsnp_oriented")], ["u", "f", "blank"])


if __name__ == "__main__":
    unittest.main()

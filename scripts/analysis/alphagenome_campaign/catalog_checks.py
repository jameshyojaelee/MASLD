#!/usr/bin/env python3
"""Scientific identity and refusal checks for the candidate evidence interface."""
import sqlite3
import tempfile
import unittest
from pathlib import Path

from catalog_build import Catalog
from catalog_query import connect, haplotype_query, query, variant_id
from catalog_infer import construct, parse_request


class IdentityChecks(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "fixture.sqlite"
        self.cat = Catalog(self.path)
        self.cat.db.execute("INSERT INTO source VALUES ('fixture','synthetic_test_only','fixture','synthetic','not_scientific_evidence')")
        self.variant = self.cat.entity("variant", "GRCh38:chr1:100:A:G")
        self.g1 = self.cat.entity("gene", "fixture_gene_1")
        self.g2 = self.cat.entity("gene", "fixture_gene_2")
        self.cat.link("fixture", self.variant, self.g1, "candidate_target")
        self.cat.link("fixture", self.variant, self.g2, "candidate_target")
        self.add_effect()
        self.cat.db.commit()

    def add_effect(self, effect=0.0):
        self.cat.evidence("fixture", "record1", "allele_effect", "measured", "reporter", effect,
                          "log2_activity_ratio", "ALT_minus_REF", [(self.variant,"variant")])

    def tearDown(self):
        self.cat.db.close()
        self.tmp.cleanup()

    def test_build_and_alleles_required(self):
        for value in ("rs1", "chr1:100:A:G", "GRCh37:chr1:100:A:G", "hg38:1:100:N:G", "hg38:1:100:A:A"):
            with self.assertRaises(ValueError):
                variant_id(value)
        self.assertEqual(variant_id("hg38:1:100:A:G"), "GRCh38:chr1:100:A:G")

    def test_missing_event_identity_does_not_form_a_shared_entity(self):
        for value in ("", " "):
            with self.assertRaisesRegex(ValueError, "Missing molecular identity"):
                self.cat.entity("event", value)

    def test_unknown_phase_never_becomes_sum(self):
        result = haplotype_query(self.path, ["hg38:1:100:A:G","hg38:1:110:C:T"], "unknown", "")
        self.assertIsNone(result["joint_effect"])
        self.assertIn("unknown_phase", result["unsupported_reason"])

    def test_known_phase_still_requires_evaluated_joint_model(self):
        result = haplotype_query(self.path, ["hg38:1:100:A:G","hg38:1:110:C:T"], "same_haplotype", "synthetic_fixture")
        self.assertIsNone(result["joint_effect"])
        self.assertEqual(result["unsupported_reason"], "joint_effect_not_stored_in_specimen")

    def test_overlapping_haplotype_edits_rejected(self):
        with self.assertRaises(ValueError):
            haplotype_query(self.path, ["hg38:1:100:AC:A","hg38:1:101:C:T"], "same_haplotype", "fixture")

    def test_missing_assay_is_not_zero(self):
        result = query(self.path,"variant","hg38:1:100:A:G", assay="H3K27ac")
        self.assertEqual(result["status"],"unsupported")
        self.assertEqual(result["total_evidence"],0)
        self.assertEqual(result["evidence"],[])

    def test_measured_zero_remains_recorded(self):
        result = query(self.path,"variant","hg38:1:100:A:G")
        self.assertEqual(result["evidence"][0]["effect"],0.0)
        self.assertEqual(result["evidence"][0]["state"],"recorded")

    def test_repeated_evidence_is_one_record(self):
        self.add_effect()
        self.cat.db.commit()
        self.assertEqual(query(self.path,"variant","hg38:1:100:A:G")["total_evidence"],1)
        with self.assertRaisesRegex(ValueError,"Conflicting duplicate"):
            self.add_effect(2)

    def test_conflicting_targets_both_returned(self):
        result = query(self.path,"variant","hg38:1:100:A:G")
        self.assertEqual({r["to_id"] for r in result["relationships"]}, {"fixture_gene_1","fixture_gene_2"})

    def test_connection_is_read_only(self):
        with connect(self.path) as db:
            with self.assertRaises(sqlite3.OperationalError):
                db.execute("DELETE FROM evidence")

    def test_phased_sequence_uses_reference_coordinates(self):
        parsed = [("chr1", 1025, "A", "G"), ("chr1", 1030, "A", "C")]
        fetch = lambda chrom, start, end: "A"*(end-start)
        reference, alternate = construct(fetch, parsed, 1, 2049)
        self.assertEqual(len(reference),len(alternate))
        self.assertEqual(alternate[1023],"G")
        self.assertEqual(alternate[1028],"C")
        self.assertEqual(sum(a != b for a,b in zip(reference,alternate)),2)
        with self.assertRaisesRegex(ValueError,"reference_mismatch"):
            construct(fetch, [("chr1",1025,"T","G")], 1, 2049)

    def test_inference_requires_phase_and_supported_edits(self):
        with self.assertRaisesRegex(ValueError,"unknown_phase"):
            parse_request(["hg38:1:1025:A:G","hg38:1:1030:A:C"],"unknown","",2048)
        with self.assertRaisesRegex(ValueError,"readout_mapping_not_validated"):
            parse_request(["hg38:1:1025:AC:A"],"unknown","",2048)


if __name__ == "__main__":
    unittest.main()

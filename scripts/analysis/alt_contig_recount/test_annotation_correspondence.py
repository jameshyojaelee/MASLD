#!/usr/bin/env python3
"""Coordinate and biological-object invariants for annotation correspondence."""
from pathlib import Path
import tempfile
import unittest
from resolve_annotation_correspondence import (Gene, alignment_blocks, project, union,
    verified_exon_segments, allele_identity_candidates, added_patch_component)
import gzip


class CoordinateTests(unittest.TestCase):
    def alignment(self, strand, target, gap, start=101, end=112):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/"alignment.gff"
            p.write_text(f'chr1\tNCBI\tmatch\t{start}\t{end}\t.\t{strand}\t.\tTarget=ALT {target};Gap={gap}\n')
            return alignment_blocks(p)[1]

    def test_deletion_does_not_project_an_alt_base_into_reference_gap(self):
        blocks = self.alignment("+", "1 10 +", "M5 D2 M5")
        mapped, covered, directions = project([(0,10)], blocks)
        self.assertEqual(mapped, [(100,105),(107,112)])
        self.assertEqual(covered,10)
        self.assertEqual(directions,{1})

    def test_insertion_is_unresolved_for_exonic_bases_in_insert(self):
        blocks = self.alignment("+", "1 14 +", "M5 I2 M7")
        self.assertEqual(project([(0,14)],blocks)[1],12)

    def test_partial_equivalent_rescue_excludes_unaligned_inserted_bases(self):
        blocks = self.alignment("+", "1 14 +", "M5 I2 M7")
        self.assertEqual(verified_exon_segments([(0,14)], blocks), [(0,5),(7,14)])

    def test_reverse_orientation_preserves_half_open_coordinates(self):
        blocks = self.alignment("-", "1 12 +", "M12")
        self.assertEqual(project([(0,3)],blocks), ([(109,112)],3,{-1}))

    def test_overlapping_isoform_exons_are_not_counted_twice(self):
        self.assertEqual(union([(0,10),(5,15),(0,10)]),[(0,15)])

    def test_malformed_gap_cannot_license_a_correspondence(self):
        with self.assertRaises(ValueError):
            self.alignment("+", "1 12 +", "M10")


class IdentityTests(unittest.TestCase):
    def setUp(self):
        self.alt = Gene("ALT", "", "lncRNA", "ALT1", "+", 0, 10)
        self.primary = Gene("P", "", "lncRNA", "chr1", "+", 100, 110)
        self.genes = {g.gene_id: g for g in [self.alt, self.primary]}
        self.groups = {"G": [("ALT", set()), ("P", set())]}

    def candidates(self, bodies=None):
        return allele_identity_candidates(self.alt, [self.primary] if bodies is None else bodies,
                                           self.genes, self.groups, {"ALT": ["G"]})

    def test_source_group_and_independent_locus_both_required(self):
        self.assertEqual(self.candidates(), [self.primary])
        self.assertEqual(self.candidates([]), [])

    def test_conflicting_HGNC_cannot_be_overruled_by_allele_group(self):
        self.alt.hgnc = "HGNC:1"; self.primary.hgnc = "HGNC:2"
        self.assertEqual(self.candidates(), [])

    def test_same_haplotype_copies_remain_distinct(self):
        self.genes["ALT2"] = Gene("ALT2", "", "lncRNA", "ALT1", "+", 20, 30)
        self.groups["G"].append(("ALT2", set()))
        self.assertEqual(self.candidates(), [])

    def test_missing_exact_GTF_member_cannot_look_like_singleton_group(self):
        self.groups["G"].append(("P.other_version", set()))
        self.assertEqual(self.candidates(), [])

    def test_added_patch_requires_positive_gap_and_component_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with gzip.open(root/"alt.scaf.agp.gz", "wt") as out:
                out.write("ALTACC\t1\t100\t1\tW\tNEWACC\t1\t100\t+\n")
            with gzip.open(root/"chr13.agp.gz", "wt") as out:
                out.write("PRIMARY\t1\t1000\t1\tN\t1000\tshort_arm\tno\tna\n")
            seq = {"unit": "PATCHES", "genbank": "ALTACC"}
            placement = {"parent_acc": "PRIMARY", "parent_name": "13", "parent_start": "200", "parent_stop": "299"}
            self.assertTrue(added_patch_component(self.alt, seq, placement, root))
            self.alt.end = 101
            self.assertFalse(added_patch_component(self.alt, seq, placement, root))

    def test_multiple_group_memberships_keep_conflicting_primary_targets(self):
        second = Gene("P2", "", "lncRNA", "chr1", "+", 100, 110)
        self.genes["P2"] = second
        self.groups["G2"] = [("ALT", set()), ("P2", set())]
        selected = allele_identity_candidates(self.alt, [self.primary, second], self.genes,
                                             self.groups, {"ALT": ["G", "G2"]})
        self.assertEqual({p.gene_id for p in selected}, {"P", "P2"})


if __name__ == "__main__": unittest.main()

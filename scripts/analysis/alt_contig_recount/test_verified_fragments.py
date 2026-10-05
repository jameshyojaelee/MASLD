#!/usr/bin/env python3
"""Scientific invariants: complete molecules, gene ambiguity, and pairing."""
import unittest
from count_verified_fragments import molecule_assignment


def record(nh=1,hi=1,gene="GENE1",mate=0,**kwargs):
    result = dict(NH=nh,HI=hi,XT=gene,XS="Assigned",XN=1,mate=mate,
                  paired=bool(mate),mate_unmapped=False,reference_id=0,
                  reference_start=100 if mate != 2 else 140,
                  next_reference_id=0,next_reference_start=140 if mate != 2 else 100)
    result.update(kwargs)
    return result


class MoleculeTests(unittest.TestCase):
    def test_two_representations_count_as_one(self):
        self.assertEqual(molecule_assignment([record(2,1),record(2,2)],False),
                         ("GENE1","accepted_same_gene_multimapper"))

    def test_incomplete_hits_are_never_renormalized(self):
        self.assertEqual(molecule_assignment([record(2,1)],False)[1],"incomplete_NH_alignment_set")

    def test_cross_gene_fragment_is_not_fractionally_distributed(self):
        self.assertIsNone(molecule_assignment([record(2,1),record(2,2,"GENE2")],False)[0])

    def test_pair_is_one_fragment(self):
        self.assertEqual(molecule_assignment([record(mate=1),record(mate=2)],True)[0],"GENE1")

    def test_duplicate_mate_is_not_an_independent_fragment(self):
        self.assertIsNone(molecule_assignment([record(mate=1),record(mate=1)],True)[0])

    def test_missing_gene_assignment_prevents_rescue(self):
        r=record(2,2);r["XS"]="Unassigned_NoFeatures";r["XN"]=-1
        self.assertIsNone(molecule_assignment([record(2,1),r],False)[0])

    def test_supplementary_segment_cannot_supply_an_extra_molecule(self):
        self.assertIsNone(molecule_assignment([record(supplementary=True)],False)[0])

    def test_missing_NH_cannot_license_a_unique_count(self):
        r=record();del r["NH"]
        self.assertIsNone(molecule_assignment([r],False)[0])

    def test_chimeric_primary_segment_cannot_pass_without_supplementary_record(self):
        self.assertEqual(molecule_assignment([record(chimeric=True)],False)[1],"chimeric_alignment")

    def test_cross_contig_pair_does_not_identify_one_gene_object(self):
        self.assertEqual(molecule_assignment([record(mate=1),record(mate=2,reference_id=1)],True)[1],"cross_contig_pair")

    def test_mate_coordinates_must_describe_the_same_hit(self):
        self.assertEqual(molecule_assignment([record(mate=1,next_reference_start=999),record(mate=2)],True)[1],"nonreciprocal_mate_coordinates")

    def test_hit_numbers_cannot_replace_missing_alignments(self):
        self.assertEqual(molecule_assignment([record(2,2),record(2,3)],False)[1],"invalid_HI_alignment_set")


if __name__=="__main__":unittest.main()

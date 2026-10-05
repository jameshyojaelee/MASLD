#!/usr/bin/env python3
"""Scientific invariants: complete molecules, gene ambiguity, and pairing."""
import unittest
from count_verified_fragments import molecule_assignment


def record(nh=1,hi=1,gene="GENE1",mate=0,**kwargs):
    return dict(NH=nh,HI=hi,XT=gene,XS="Assigned",XN=1,mate=mate,**kwargs)


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


if __name__=="__main__":unittest.main()

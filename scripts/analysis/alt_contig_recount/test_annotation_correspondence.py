#!/usr/bin/env python3
"""Coordinate and biological-object invariants for annotation correspondence."""
from pathlib import Path
import tempfile
import unittest
from resolve_annotation_correspondence import alignment_blocks, project, union


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

    def test_reverse_orientation_preserves_half_open_coordinates(self):
        blocks = self.alignment("-", "1 12 +", "M12")
        self.assertEqual(project([(0,3)],blocks), ([(109,112)],3,{-1}))

    def test_overlapping_isoform_exons_are_not_counted_twice(self):
        self.assertEqual(union([(0,10),(5,15),(0,10)]),[(0,15)])

    def test_malformed_gap_cannot_license_a_correspondence(self):
        with self.assertRaises(ValueError):
            self.alignment("+", "1 12 +", "M10")


if __name__ == "__main__": unittest.main()

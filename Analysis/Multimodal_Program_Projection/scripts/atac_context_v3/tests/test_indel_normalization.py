#!/usr/bin/env python3
"""Regression fixtures for reference-preserving indel left normalization."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from genetics_context import normalize_against_reference


class IndelNormalizationTests(unittest.TestCase):
    def normalized(self, sequence, position, ref, alt):
        def fetch(_chrom, start, end):
            self.assertGreaterEqual(start, 0)
            return sequence[start:end]

        result = normalize_against_reference("chr1", position, ref, alt, fetch)
        self.assertIsNotNone(result)
        return result

    def assert_haplotype_preserved(self, sequence, position, ref, alt, expected):
        result = self.normalized(sequence, position, ref, alt)
        self.assertEqual((result.position_1based, result.ref, result.alt), expected)
        original = sequence[:position - 1] + alt + sequence[position - 1 + len(ref):]
        normalized = (sequence[:result.position_1based - 1] + result.alt
                      + sequence[result.position_1based - 1 + len(result.ref):])
        self.assertEqual(normalized, original)

    def test_homopolymer_reaches_preceding_different_anchor(self):
        self.assert_haplotype_preserved("CAAAAG", 3, "AA", "A", (1, "CA", "C"))
        self.assert_haplotype_preserved("CAAAAG", 3, "A", "AA", (1, "C", "CA"))

    def test_tandem_repeat_rotates_alleles_to_leftmost_anchor(self):
        self.assert_haplotype_preserved("GATATATC", 4, "ATA", "A", (1, "GAT", "G"))
        self.assert_haplotype_preserved("GATATATC", 4, "A", "ATA", (1, "G", "GAT"))

    def test_chromosome_start_keeps_nonempty_anchor(self):
        self.assert_haplotype_preserved("AAAAC", 1, "AA", "A", (1, "AA", "A"))
        self.assert_haplotype_preserved("AAAAC", 1, "A", "AA", (1, "A", "AA"))

    def test_mnv_trimming_and_snp_are_unchanged(self):
        self.assert_haplotype_preserved("AACCGGTT", 3, "CC", "TT", (3, "CC", "TT"))
        self.assert_haplotype_preserved("AACCGGTT", 3, "CCG", "CTG", (4, "C", "T"))
        self.assert_haplotype_preserved("AACCGGTT", 3, "C", "A", (3, "C", "A"))

    def test_identifiable_swap_and_reverse_complement(self):
        for alleles in (("A", "AT"), ("AT", "A"), ("T", "AT")):
            result = self.normalized("GACCC", 2, *alleles)
            self.assertEqual((result.position_1based, result.ref, result.alt), (2, "A", "AT"))
        swapped_snp = self.normalized("AACCGGTT", 3, "A", "C")
        complemented = self.normalized("AACCGGTT", 3, "G", "T")
        self.assertEqual((swapped_snp.ref, swapped_snp.alt), ("C", "A"))
        self.assertEqual((complemented.ref, complemented.alt), ("C", "A"))
        self.assertEqual(complemented.orientation, "reverse_complement")

    def test_unmatched_alleles_and_invalid_position_are_rejected(self):
        fetch = lambda _chrom, start, end: "AACCGGTT"[start:end]
        self.assertIsNone(normalize_against_reference("chr1", 3, "A", "T", fetch))
        self.assertIsNone(normalize_against_reference("chr1", 0, "A", "T", fetch))
        self.assertIsNone(normalize_against_reference("chr1", -1, "A", "T", fetch))


if __name__ == "__main__":
    unittest.main(verbosity=2)

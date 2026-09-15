from __future__ import annotations

import unittest

from scripts.build_gse281364_enformer_fixture import (
    CENTER,
    EnformerFixtureError,
    WINDOW,
    allele_window,
    reverse_complement,
)


class _Reference:
    def __init__(self, sequence: str):
        self.sequence = sequence

    def fetch(self, contig: str, start: int, end: int) -> str:
        if contig != "chr1" or start != 10 or end != 10 + WINDOW:
            raise AssertionError("unexpected reference request")
        return self.sequence


class EnformerFixtureTests(unittest.TestCase):
    def test_centered_ref_alt_and_reverse_complement(self) -> None:
        sequence = "A" * CENTER + "G" + "C" * (WINDOW - CENTER - 1)
        row = {
            "variant_pos0": str(10 + CENTER),
            "contig": "chr1",
            "genomic_ref": "G",
            "genomic_alt": "T",
        }
        reference, alternative, start, end = allele_window(_Reference(sequence), row)
        self.assertEqual((start, end), (10, 10 + WINDOW))
        self.assertEqual(reference[CENTER], "G")
        self.assertEqual(alternative[CENTER], "T")
        self.assertEqual(sum(a != b for a, b in zip(reference, alternative)), 1)
        self.assertEqual(reverse_complement("ACGTN"), "NACGT")

    def test_ref_mismatch_fails_closed(self) -> None:
        sequence = "A" * WINDOW
        row = {
            "variant_pos0": str(10 + CENTER),
            "contig": "chr1",
            "genomic_ref": "G",
            "genomic_alt": "T",
        }
        with self.assertRaises(EnformerFixtureError):
            allele_window(_Reference(sequence), row)


if __name__ == "__main__":
    unittest.main()

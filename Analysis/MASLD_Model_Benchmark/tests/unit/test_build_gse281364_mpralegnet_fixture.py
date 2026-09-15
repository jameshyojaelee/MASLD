from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from scripts.alphagenome_sei_build_fixture import IndexedFasta
from scripts.build_gse281364_mpralegnet_fixture import (
    FORWARD_VARIANT_INDEX0,
    REVERSE_VARIANT_INDEX0,
    WINDOW,
    allele_window,
    reverse_complement,
)


class GSE281364MPRALegNetFixtureTests(unittest.TestCase):
    def test_ref_alt_geometry_and_orientation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sequence = "A" * 10_000
            fasta = root / "reference.fa"
            fasta.write_text(">chr1\n" + sequence + "\n")
            fai = root / "reference.fa.fai"
            fai.write_text(
                f"chr1\t{len(sequence)}\t6\t{len(sequence)}\t{len(sequence)+1}\n"
            )
            reference = IndexedFasta(fasta, fai)
            row = {
                "variant_pos0": "5000",
                "contig": "chr1",
                "genomic_ref": "A",
                "genomic_alt": "C",
                "ref_sequence_107bp": "A" * 107,
                "alt_sequence_107bp": "A" * 62 + "C" + "A" * 44,
                "oligo_difference_index0": "62",
            }
            ref, alt, start, end = allele_window(reference, row)
            self.assertEqual((len(ref), len(alt), end - start), (WINDOW, WINDOW, WINDOW))
            self.assertEqual((ref[FORWARD_VARIANT_INDEX0], alt[FORWARD_VARIANT_INDEX0]), ("A", "C"))
            self.assertEqual(sum(left != right for left, right in zip(ref, alt)), 1)
            self.assertEqual(
                reverse_complement(alt)[REVERSE_VARIANT_INDEX0], "G"
            )


if __name__ == "__main__":
    unittest.main()

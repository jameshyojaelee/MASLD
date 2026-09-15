from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from scripts.alphagenome_sei_build_fixture import IndexedFasta
from scripts.build_gse281364_sei_fixture import CENTER, WINDOW, allele_window


class GSE281364SeiFixtureTests(unittest.TestCase):
    def test_ref_alt_geometry(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sequence = "A" * 10_000
            fasta = root / "reference.fa"
            fasta.write_text(">chr1\n" + sequence + "\n")
            fai = root / "reference.fa.fai"
            fai.write_text(f"chr1\t{len(sequence)}\t6\t{len(sequence)}\t{len(sequence)+1}\n")
            reference = IndexedFasta(fasta, fai)
            row = {
                "variant_pos0": "5000",
                "contig": "chr1",
                "genomic_ref": "A",
                "genomic_alt": "C",
            }
            ref, alt, start, end = allele_window(reference, row)
            self.assertEqual((len(ref), len(alt), end - start), (WINDOW, WINDOW, WINDOW))
            self.assertEqual((ref[CENTER], alt[CENTER]), ("A", "C"))
            self.assertEqual(sum(a != b for a, b in zip(ref, alt)), 1)


if __name__ == "__main__":
    unittest.main()

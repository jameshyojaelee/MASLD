from __future__ import annotations

from hashlib import sha256
import tempfile
from pathlib import Path
import unittest

import numpy as np

from scripts.caduceus_predict_gse281364_131k import (
    ALLELES,
    CENTER,
    WINDOW,
    CaduceusPredictionError,
    allele_sequences,
    digest_text,
    load_chunk,
    reverse_complement,
    write_chunk,
)


class _Reference:
    def __init__(self, sequence: str):
        self.sequence = sequence

    def fetch(self, contig: str, start: int, end: int) -> str:
        if (contig, start, end) != ("chr1", 10, 10 + WINDOW):
            raise AssertionError("unexpected reference request")
        return self.sequence


class CaduceusGSE281364PredictionTests(unittest.TestCase):
    def test_allele_and_external_rc_contract(self) -> None:
        reference = "A" * CENTER + "G" + "C" * (WINDOW - CENTER - 1)
        alternative = reference[:CENTER] + "T" + reference[CENTER + 1 :]
        row = {
            "input_start0": "10",
            "input_end0": str(10 + WINDOW),
            "variant_pos0": str(10 + CENTER),
            "contig": "chr1",
            "ref": "G",
            "alt": "T",
            "reference_sequence_sha256": digest_text(reference),
            "alternative_sequence_sha256": digest_text(alternative),
            "reverse_complement_reference_sha256": digest_text(reverse_complement(reference)),
            "reverse_complement_alternative_sha256": digest_text(reverse_complement(alternative)),
        }
        sequences = allele_sequences(_Reference(reference), row)
        self.assertEqual(len(sequences), len(ALLELES))
        self.assertEqual(sequences[2], reverse_complement(sequences[0]))
        self.assertEqual(sequences[3], reverse_complement(sequences[1]))
        self.assertEqual(sequences[2][CENTER - 1], "C")
        self.assertEqual(sequences[3][CENTER - 1], "A")

    def test_ref_or_hash_mismatch_fails_closed(self) -> None:
        reference = "A" * WINDOW
        row = {
            "input_start0": "10",
            "input_end0": str(10 + WINDOW),
            "variant_pos0": str(10 + CENTER),
            "contig": "chr1",
            "ref": "G",
            "alt": "T",
            "reference_sequence_sha256": sha256(reference.encode()).hexdigest(),
        }
        with self.assertRaises(CaduceusPredictionError):
            allele_sequences(_Reference(reference), row)

    def test_chunk_roundtrip_and_tamper_rejection(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            values = np.zeros((2, 4, 256), dtype=np.float32)
            chunk = write_chunk(root, 0, 2, values, 0.0)
            observed, receipt = load_chunk(chunk, 0, 2)
            np.testing.assert_array_equal(observed, values)
            self.assertEqual(receipt["max_external_rc_feature_abs_difference"], 0.0)
            (chunk / "embeddings.npy").write_bytes(b"tampered")
            with self.assertRaises(CaduceusPredictionError):
                load_chunk(chunk, 0, 2)


if __name__ == "__main__":
    unittest.main()

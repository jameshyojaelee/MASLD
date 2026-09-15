from __future__ import annotations

import unittest

import numpy as np
from scipy import sparse

from scripts.build_transcriptformer_unlabeled_fixture import (
    ALLOWED_SOURCE_OBS,
    COUNT_CLIP,
    PAD_TOKEN_ID,
    SEQUENCE_LENGTH,
    _outer_fold,
    _tokenize,
)


class TranscriptFormerUnlabeledFixtureTests(unittest.TestCase):
    def test_native_and_common_gene_orders_are_distinct_and_exact(self) -> None:
        genes = ["ENSG000003", "ENSG000001", "ENSG000002", "ENSG000004"]
        vocabulary = {gene: token for token, gene in enumerate(genes, start=7)}
        matrix = sparse.csr_matrix(np.asarray([[1, 40, 40, 0]], dtype=np.int64))

        native_tokens, native_counts, _ = _tokenize(
            matrix, genes, vocabulary, "native"
        )
        common_tokens, common_counts, _ = _tokenize(
            matrix, genes, vocabulary, "common"
        )

        self.assertEqual(native_tokens.shape, (1, SEQUENCE_LENGTH))
        self.assertEqual(native_tokens[0, :4].tolist(), [7, 8, 9, PAD_TOKEN_ID])
        self.assertEqual(common_tokens[0, :4].tolist(), [8, 9, 7, PAD_TOKEN_ID])
        self.assertEqual(native_counts[0, :3].tolist(), [1.0, COUNT_CLIP, COUNT_CLIP])
        self.assertEqual(common_counts[0, :3].tolist(), [COUNT_CLIP, COUNT_CLIP, 1.0])

    def test_donor_fold_is_deterministic_and_source_outcomes_are_not_allowed(self) -> None:
        self.assertEqual(_outer_fold("donor-1"), _outer_fold("donor-1"))
        self.assertTrue(0 <= _outer_fold("donor-1") < 5)
        self.assertEqual(
            _outer_fold("donor-1", namespace="shared-screen-v1"),
            _outer_fold("donor-1", namespace="shared-screen-v1"),
        )
        with self.assertRaisesRegex(ValueError, "namespace is empty"):
            _outer_fold("donor-1", namespace="")
        self.assertEqual(ALLOWED_SOURCE_OBS, ("_index", "donor_id", "dataset"))

    def test_negative_or_noninteger_counts_are_rejected(self) -> None:
        genes = ["ENSG000001"]
        vocabulary = {genes[0]: 7}
        for value in (-1.0, 1.5):
            with self.subTest(value=value), self.assertRaisesRegex(
                ValueError, "raw nonnegative integer"
            ):
                _tokenize(
                    sparse.csr_matrix(np.asarray([[value]], dtype=np.float32)),
                    genes,
                    vocabulary,
                    "common",
                )


if __name__ == "__main__":
    unittest.main()

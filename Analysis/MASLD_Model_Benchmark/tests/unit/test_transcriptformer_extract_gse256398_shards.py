from __future__ import annotations

import unittest

import numpy as np

from scripts.transcriptformer_extract_gse256398_shards import (
    GSE256398TranscriptFormerEmbeddingError,
    SEQUENCE_LENGTH,
    validate_shard_arrays,
)


class TranscriptFormerGSE256398ShardExtractionTests(unittest.TestCase):
    def test_shard_contract_accepts_valid_arrays(self) -> None:
        validate_shard_arrays(
            tokens=np.ones((2, SEQUENCE_LENGTH), dtype=np.int32),
            counts=np.ones((2, SEQUENCE_LENGTH), dtype=np.float32),
            assay=np.zeros((2, 1), dtype=np.int64),
            folds=np.asarray([0, 0], dtype=np.int8),
            row_ids=np.asarray(["a", "b"]),
            source_mask=np.asarray([True, False], dtype=bool),
        )

    def test_shard_contract_rejects_duplicate_rows(self) -> None:
        with self.assertRaises(GSE256398TranscriptFormerEmbeddingError):
            validate_shard_arrays(
                tokens=np.ones((2, SEQUENCE_LENGTH), dtype=np.int32),
                counts=np.ones((2, SEQUENCE_LENGTH), dtype=np.float32),
                assay=np.zeros((2, 1), dtype=np.int64),
                folds=np.asarray([0, 0], dtype=np.int8),
                row_ids=np.asarray(["a", "a"]),
                source_mask=np.asarray([True, False], dtype=bool),
            )


if __name__ == "__main__":
    unittest.main()

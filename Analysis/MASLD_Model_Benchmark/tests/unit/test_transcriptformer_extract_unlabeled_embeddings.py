from __future__ import annotations

import csv
from pathlib import Path
import tempfile
import unittest

from scripts.transcriptformer_extract_unlabeled_embeddings import (
    TranscriptFormerEmbeddingError,
    _read_rows,
)


class TranscriptFormerEmbeddingContractTests(unittest.TestCase):
    def test_row_contract_accepts_configured_cardinality_and_donor_folds(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "rows.tsv"
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=(
                        "row_position",
                        "row_id",
                        "donor_id",
                        "dataset_id",
                        "outer_fold",
                    ),
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(
                    (
                        {
                            "row_position": 0,
                            "row_id": "cell-0",
                            "donor_id": "donor-a",
                            "dataset_id": "development",
                            "outer_fold": 1,
                        },
                        {
                            "row_position": 1,
                            "row_id": "cell-1",
                            "donor_id": "donor-a",
                            "dataset_id": "development",
                            "outer_fold": 1,
                        },
                    )
                )
            self.assertEqual(len(_read_rows(path, expected_rows=2)), 2)
            with self.assertRaises(TranscriptFormerEmbeddingError):
                _read_rows(path, expected_rows=3)


if __name__ == "__main__":
    unittest.main()

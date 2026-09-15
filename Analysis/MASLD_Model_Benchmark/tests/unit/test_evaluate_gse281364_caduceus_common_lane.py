from __future__ import annotations

import csv
from pathlib import Path
import tempfile
import unittest

from scripts.evaluate_gse281364_caduceus_common_lane import (
    CommonLaneError,
    load_fixture,
)


class CaduceusCommonLaneEvaluationTests(unittest.TestCase):
    def _fixture(self, root: Path, rows: int = 1_033) -> None:
        target = root / "fixture"
        target.mkdir()
        fields = (
            "fixture_id",
            "element_id",
            "outer_locus_sequence_group_id",
            "outer_fold",
            "input_length_bp",
            "ref",
            "alt",
        )
        with (target / "sequence_manifest.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for index in range(rows):
                writer.writerow(
                    {
                        "fixture_id": f"f{index}",
                        "element_id": f"e{index}",
                        "outer_locus_sequence_group_id": f"g{index}",
                        "outer_fold": index % 5,
                        "input_length_bp": 131_072,
                        "ref": "A",
                        "alt": "C",
                    }
                )

    def test_fixture_builds_same_allele_identity_baseline(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._fixture(root)
            rows, features = load_fixture(root)
            self.assertEqual(len(rows), 1_033)
            self.assertEqual(features.shape, (1_033, 16))
            self.assertTrue((features.sum(axis=1) == 1).all())

    def test_fixture_census_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._fixture(root, rows=10)
            with self.assertRaises(CommonLaneError):
                load_fixture(root)


if __name__ == "__main__":
    unittest.main()

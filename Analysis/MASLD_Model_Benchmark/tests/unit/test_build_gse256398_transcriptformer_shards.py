from __future__ import annotations

import unittest

from scripts.build_gse256398_transcriptformer_shards import donor_sort_key, outer_fold


class GSE256398TranscriptFormerShardTests(unittest.TestCase):
    def test_numeric_donor_sorting(self) -> None:
        self.assertEqual(sorted(["S20", "S3", "S1"], key=donor_sort_key), ["S1", "S3", "S20"])

    def test_outer_fold_is_deterministic_and_bounded(self) -> None:
        self.assertEqual(outer_fold("S1"), outer_fold("S1"))
        self.assertIn(outer_fold("S1"), range(5))


if __name__ == "__main__":
    unittest.main()

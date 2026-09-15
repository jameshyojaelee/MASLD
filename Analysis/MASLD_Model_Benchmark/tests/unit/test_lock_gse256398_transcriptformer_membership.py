from __future__ import annotations

import unittest

from scripts.lock_gse256398_transcriptformer_membership import canonical, counts_by_donor


class GSE256398TranscriptFormerMembershipLockTests(unittest.TestCase):
    def test_canonical_membership_uses_numeric_donor_then_barcode_order(self) -> None:
        rows = [
            {"row_id": "gse256398:S10:B", "gsm": "G2", "source_sample_id": "S10", "barcode": "B"},
            {"row_id": "gse256398:S2:C", "gsm": "G1", "source_sample_id": "S2", "barcode": "C"},
            {"row_id": "gse256398:S2:A", "gsm": "G1", "source_sample_id": "S2", "barcode": "A"},
        ]
        self.assertEqual(
            [row["row_id"] for row in canonical(rows)],
            ["gse256398:S2:A", "gse256398:S2:C", "gse256398:S10:B"],
        )

    def test_counts_by_donor_preserves_biological_unit(self) -> None:
        rows = [
            {"source_sample_id": "S1"},
            {"source_sample_id": "S1"},
            {"source_sample_id": "S2"},
        ]
        self.assertEqual(counts_by_donor(rows), {"S1": 2, "S2": 1})


if __name__ == "__main__":
    unittest.main()

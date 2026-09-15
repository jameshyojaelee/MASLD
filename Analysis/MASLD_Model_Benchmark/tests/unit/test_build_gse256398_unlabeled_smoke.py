from __future__ import annotations

import unittest

from scripts.build_gse256398_unlabeled_smoke import priority, select_balanced


class GSE256398UnlabeledSmokeTests(unittest.TestCase):
    def test_selection_is_exactly_forty_per_donor_and_deterministic(self) -> None:
        rows = []
        for donor in range(1, 27):
            for nucleus in range(50):
                row_id = f"gse256398:S{donor}:bc-{nucleus}"
                rows.append(
                    {
                        "row_id": row_id,
                        "gsm": f"GSM{donor}",
                        "source_sample_id": f"S{donor}",
                        "barcode": f"bc-{nucleus}",
                    }
                )
        selected = select_balanced(rows, {row["row_id"] for row in rows})
        self.assertEqual(len(selected), 1040)
        self.assertEqual(len({row["row_id"] for row in selected}), 1040)
        counts = {}
        for row in selected:
            counts[row["source_sample_id"]] = counts.get(row["source_sample_id"], 0) + 1
        self.assertEqual(set(counts.values()), {40})
        s1 = [row for row in selected if row["source_sample_id"] == "S1"]
        self.assertEqual(
            [row["row_id"] for row in s1],
            [row["row_id"] for row in sorted(s1, key=lambda row: priority(row["row_id"]))],
        )


if __name__ == "__main__":
    unittest.main()

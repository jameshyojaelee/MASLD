from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import unittest

from scripts.build_gse281364_outcome_blind_splits import (
    SplitContractError,
    assign_groups,
)


def row(group: int, contig: str, position: int) -> dict[str, str]:
    return {
        "outer_locus_sequence_group_id": f"outer_{group:020x}",
        "contig": contig,
        "variant_pos0": str(position),
    }


class GSE281364OutcomeBlindSplitTests(unittest.TestCase):
    def test_assignment_is_deterministic_group_safe_and_balanced(self) -> None:
        rows: list[dict[str, str]] = []
        for group in range(1_033):
            size = 1 + (group % 7 == 0)
            rows.extend(row(group, f"chr{1 + group % 4}", group * 10_000 + index) for index in range(size))
        first = assign_groups(rows, 5, "fixture-seed")
        second = assign_groups(list(reversed(rows)), 5, "fixture-seed")
        self.assertEqual(first, second)
        self.assertEqual(set(first.values()), set(range(5)))
        counts = defaultdict(int)
        for item in rows:
            counts[first[item["outer_locus_sequence_group_id"]]] += 1
        self.assertLessEqual(max(counts.values()) - min(counts.values()), 3)

    def test_rejects_wrong_number_of_groups(self) -> None:
        with self.assertRaises(SplitContractError):
            assign_groups([row(group, "chr1", group * 10_000) for group in range(10)], 5, "x")


if __name__ == "__main__":
    unittest.main()

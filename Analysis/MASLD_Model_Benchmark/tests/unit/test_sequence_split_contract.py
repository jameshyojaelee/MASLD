from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from scripts import build_sequence_split_contract as builder


class SequenceSplitContractTests(unittest.TestCase):
    def test_balanced_fold_assignment_is_exact_and_deterministic(self) -> None:
        sizes = [(f"chr{index}", 1000 - index) for index in range(1, 23)]
        sizes.extend((("chrX", 500), ("chrY", 200)))
        first = builder.assign_balanced_folds(sizes)
        second = builder.assign_balanced_folds(tuple(reversed(sizes)))
        self.assertEqual(first, second)
        self.assertEqual(set(first), {contig for contig, _ in sizes})
        self.assertEqual(set(first.values()), set(range(5)))

    def test_ccre_window_geometry_and_hash_are_stable(self) -> None:
        folds = {"chr1": 3}
        sizes = {"chr1": 10_000}
        fields = [
            "chr1",
            "4500",
            "4700",
            "EH38E0000001",
            "0",
            ".",
            "4500",
            "4700",
            "255,0,0",
            "PLS",
        ]
        value = builder._candidate_from_fields(
            fields,
            fold_by_contig=folds,
            chrom_size_by_contig=sizes,
            seed=builder.SEED,
        )
        self.assertIsNotNone(value)
        assert value is not None
        self.assertEqual(value[6], 3)
        self.assertEqual((value[7], value[8]), (4100, 5100))
        self.assertEqual((value[9], value[10]), (3543, 5657))
        self.assertEqual(value[0], builder.digest_integer("PLS", "EH38E0000001"))

    def test_overlap_guard_uses_full_input_window(self) -> None:
        occupied: dict[str, list[tuple[int, int]]] = {}
        first = (1, "a", "chr1", 0, 1, "PLS", 0, 0, 1000, 100, 2214)
        overlap = (2, "b", "chr1", 0, 1, "pELS", 0, 0, 1000, 2213, 4327)
        adjacent = (3, "c", "chr1", 0, 1, "dELS", 0, 0, 1000, 2214, 4328)
        self.assertTrue(builder._add_if_nonoverlapping(occupied, first))
        self.assertFalse(builder._add_if_nonoverlapping(occupied, overlap))
        self.assertTrue(builder._add_if_nonoverlapping(occupied, adjacent))

    def test_ambiguous_input_windows_are_excluded_before_selection(self) -> None:
        valid = (1, "valid", "chr1", 0, 1, "PLS", 0, 0, 1000, 100, 2214)
        ambiguous = (2, "ambiguous", "chr1", 0, 1, "PLS", 0, 0, 1000, 3000, 5114)
        heaps = {
            (0, "PLS"): [
                (-1, "valid", valid),
                (-2, "ambiguous", ambiguous),
            ]
        }

        def fetch_sequence(_contig: str, start: int, end: int) -> str:
            sequence = "A" * (end - start)
            return sequence if start == 100 else sequence[:-1] + "N"

        filtered, excluded = builder.filter_candidate_heaps_by_sequence(
            heaps, fetch_sequence=fetch_sequence
        )
        self.assertEqual([item[1] for item in filtered[(0, "PLS")]], ["valid"])
        self.assertEqual(excluded[(0, "PLS")], 1)

    def test_chrom_size_contract_rejects_chrM(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "sizes.tsv"
            rows = [(f"chr{index}", 10_000) for index in range(1, 23)]
            rows.extend((("chrX", 10_000), ("chrY", 10_000), ("chrM", 16_569)))
            path.write_text(
                "".join(f"{contig}\t{size}\n" for contig, size in rows),
                encoding="utf-8",
            )
            with self.assertRaises(builder.SequenceSplitError):
                builder.read_chrom_sizes(path)


if __name__ == "__main__":
    unittest.main()

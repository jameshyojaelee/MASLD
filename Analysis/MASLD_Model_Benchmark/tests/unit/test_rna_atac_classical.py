from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from masld_bench.adapters.rna_atac_classical import (
    RNAATACAdapterError,
    deterministic_peak_indices,
    fold_index,
    join_hash,
)


class RNAATACClassicalContractTests(unittest.TestCase):
    def test_donor_folds_are_deterministic_and_grouped(self) -> None:
        donors = ["donor-1", "donor-1", "donor-2", "donor-3"]
        first = [fold_index(item, seed=20260821, outer_folds=5) for item in donors]
        second = [fold_index(item, seed=20260821, outer_folds=5) for item in donors]
        self.assertEqual(first, second)
        self.assertEqual(first[0], first[1])
        self.assertTrue(all(0 <= item < 5 for item in first))
        with self.assertRaises(RNAATACAdapterError):
            fold_index("donor", seed=1, outer_folds=1)

    def test_smoke_peak_selection_uses_only_stable_peak_ids(self) -> None:
        peak_ids = [f"chr1:{index}-{index + 100}" for index in range(100)]
        first = deterministic_peak_indices(peak_ids, 20)
        second = deterministic_peak_indices(list(peak_ids), 20)
        self.assertEqual(first, second)
        self.assertEqual(len(first), 20)
        self.assertEqual(tuple(sorted(first)), first)
        with self.assertRaises(RNAATACAdapterError):
            deterministic_peak_indices(["duplicate", "duplicate"], 1)

    def test_join_hash_separates_rows_donors_blocks_and_namespaces(self) -> None:
        values = {
            join_hash("rna-atac:v1", "row", "d1\0hep\0p1"),
            join_hash("rna-atac:v1", "unit", "d1"),
            join_hash("rna-atac:v1", "block", "chr1"),
            join_hash("rna-atac:v2", "unit", "d1"),
        }
        self.assertEqual(len(values), 4)
        self.assertTrue(all(len(value) == 64 for value in values))


if __name__ == "__main__":
    unittest.main()

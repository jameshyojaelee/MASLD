from __future__ import annotations

import unittest

from scripts.build_atac_transport_membership import (
    ATACMembershipError,
    fold_index,
    resolve_raw_barcodes,
)


class ATACTransportMembershipTests(unittest.TestCase):
    def test_resolves_anndata_make_unique_suffixes(self) -> None:
        raw = ["AAAC-1", "TTTT-1"]
        self.assertEqual(
            resolve_raw_barcodes(["AAAC-1-2", "TTTT-1-7"], raw), raw
        )

    def test_exact_raw_ids_are_retained(self) -> None:
        raw = ["cell_A", "cell_B"]
        self.assertEqual(resolve_raw_barcodes(raw, raw), raw)

    def test_rejects_nonbijective_cell_universe(self) -> None:
        with self.assertRaises(ATACMembershipError):
            resolve_raw_barcodes(["cell_A"], ["cell_A", "cell_B"])

    def test_fold_is_deterministic_and_cohort_namespaced(self) -> None:
        first = fold_index("gse244832", "D01")
        second = fold_index("gse281367", "D01")
        self.assertEqual(first, fold_index("gse244832", "D01"))
        self.assertIn(first, range(5))
        self.assertIn(second, range(5))


if __name__ == "__main__":
    unittest.main()

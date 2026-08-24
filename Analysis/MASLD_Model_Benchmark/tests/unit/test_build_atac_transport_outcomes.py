from __future__ import annotations

import gzip
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.build_atac_transport_outcomes import (
    ATACOutcomeError,
    WindowIndex,
    process_donor,
    read_windows,
    tn5_positions,
    windows_at,
)


class ATACTransportOutcomeTests(unittest.TestCase):
    def _windows(self) -> WindowIndex:
        rows = (
            ("chr1", 100, 1100, "w1", "test"),
            ("chr1", 600, 1600, "w2", "test"),
        )
        return WindowIndex(rows, {"chr1": ((100, 600), (1100, 1600), (0, 1))})

    def test_overlap_lookup_is_half_open(self) -> None:
        windows = self._windows()
        self.assertEqual(sorted(windows_at(windows, "chr1", 700)), [0, 1])
        self.assertEqual(list(windows_at(windows, "chr1", 1600)), [])

    def test_assay_native_tn5_coordinate_contracts(self) -> None:
        self.assertEqual(tn5_positions("gse281367", 100, 120), (100, 119))
        self.assertEqual(tn5_positions("gse244832", 100, 120), (104, 115))
        with self.assertRaises(ATACOutcomeError):
            tn5_positions("gse244832", 100, 105)

    def test_reads_fixed_windows_and_rejects_wrong_role(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "regions.bed"
            path.write_text(
                "chr1\t100\t1100\tfixed_test|w1|PLS\t0\t.\t0\t-1\t-1\t500\n",
                encoding="utf-8",
            )
            parsed = read_windows(
                path,
                role_by_contig={"chr1": "test"},
                reference_lengths={"chr1": 2000},
                expected_per_role=None,
            )
            self.assertEqual(len(parsed.rows), 1)
            with self.assertRaises(ATACOutcomeError):
                read_windows(
                    path,
                    role_by_contig={"chr1": "valid"},
                    reference_lengths={"chr1": 2000},
                    expected_per_role=None,
                )

    def test_fragment_unit_and_multiplicity_are_separate(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fragments.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("chr1\t96\t111\tcell_a\t7\n")
                handle.write("chr1\t700\t720\tcell_b\t3\n")
                handle.write("chr1\t800\t820\tnot_selected\t4\n")
            result = process_donor(
                dataset_id="gse244832",
                donor_id="D01",
                fragment_path=path,
                expected_size=path.stat().st_size,
                primary_by_barcode={"cell_a": "hepatocyte", "cell_b": "hepatocyte"},
                all_barcodes={"cell_a", "cell_b", "not_selected"},
                windows=self._windows(),
                reference_lengths={"chr1": 2000},
                profile_bins=20,
            )
            hepatocyte = 2
            self.assertEqual(result.genome_unit_fragments[hepatocyte], 2)
            self.assertEqual(result.genome_multiplicity_fragments[hepatocyte], 10)
            self.assertEqual(result.unit_total[hepatocyte, 0], 4)
            self.assertEqual(result.unit_total[hepatocyte, 1], 2)
            self.assertEqual(result.multiplicity_total[hepatocyte, 0], 20)
            self.assertEqual(result.multiplicity_total[hepatocyte, 1], 6)
            np.testing.assert_array_equal(
                result.unit_profile.sum(axis=2, dtype=np.uint64), result.unit_total
            )

    def test_missing_primary_barcode_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fragments.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("chr1\t100\t120\tcell_a\t1\n")
            with self.assertRaises(ATACOutcomeError):
                process_donor(
                    dataset_id="gse281367",
                    donor_id="Z01",
                    fragment_path=path,
                    expected_size=path.stat().st_size,
                    primary_by_barcode={"cell_a": "hepatocyte", "cell_b": "hepatocyte"},
                    all_barcodes={"cell_a", "cell_b"},
                    windows=self._windows(),
                    reference_lengths={"chr1": 2000},
                    profile_bins=20,
                )


if __name__ == "__main__":
    unittest.main()

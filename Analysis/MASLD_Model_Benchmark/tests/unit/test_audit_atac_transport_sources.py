from __future__ import annotations

import csv
import gzip
from pathlib import Path
import tempfile
import unittest

from scripts.audit_atac_transport_sources import audit_peak_names, audit_pseudobulk


class ATACTransportSourceAuditTests(unittest.TestCase):
    def test_peak_name_audit_detects_duplicates_and_bad_intervals(self) -> None:
        audit = audit_peak_names(
            ["chr1:10-20", "chr1:10-20", "chrX:30-40", "chr1:50-40", "bad"]
        )
        self.assertEqual(audit["feature_columns"], 5)
        self.assertEqual(audit["unique_feature_columns"], 4)
        self.assertEqual(audit["duplicate_feature_instances"], 1)
        self.assertEqual(audit["invalid_intervals"], 1)
        self.assertEqual(audit["malformed_feature_ids"], 1)

    def test_pseudobulk_rejects_duplicate_feature_columns(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            counts = root / "counts.tsv.gz"
            coldata = root / "coldata.tsv"
            with gzip.open(counts, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["donor_id", "chr1:10-20", "chr1:10-20"])
                writer.writerow(["D1", "2", "3"])
                writer.writerow(["D2", "0", "4"])
            coldata.write_text(
                "donor_id\tcondition\tlibsize\nD1\tA\t5\nD2\tB\t4\n",
                encoding="utf-8",
            )
            audit = audit_pseudobulk(counts, coldata)
            self.assertEqual(audit["donor_rows"], 2)
            self.assertEqual(audit["duplicate_feature_instances"], 1)
            self.assertFalse(audit["usable_as_benchmark_outcome"])

    def test_pseudobulk_accepts_unique_integer_matrix(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            counts = root / "counts.tsv.gz"
            coldata = root / "coldata.tsv"
            with gzip.open(counts, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["donor_id", "chr1:10-20", "chrX:30-40"])
                writer.writerow(["D1", "2", "3"])
            coldata.write_text(
                "donor_id\tcondition\tlibsize\nD1\tA\t5\n", encoding="utf-8"
            )
            audit = audit_pseudobulk(counts, coldata)
            self.assertEqual(audit["duplicate_feature_instances"], 0)
            self.assertTrue(audit["usable_as_benchmark_outcome"])


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import csv
import gzip
from pathlib import Path
import tempfile
import unittest

from scripts.audit_gse267145_rna_measurement import RNAMeasurementError, audit_matrix


class GSE267145RNAMeasurementTests(unittest.TestCase):
    def test_fractional_nonnegative_matrix_passes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "matrix.tsv.gz"
            with gzip.open(path, "wt", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["ensembl_gene_id", "A", "B"])
                writer.writerow(["ENSG00000000001", "1.5", "2"])
                writer.writerow(["ENSG00000000002", "0", "3"])
            result = audit_matrix(path, expected_rows=2, expected_samples=2)
        self.assertFalse(result["all_integer"])
        self.assertEqual(result["fractional_values"], 1)

    def test_integer_only_matrix_fails_semantic_gate(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "matrix.tsv.gz"
            with gzip.open(path, "wt", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["ensembl_gene_id", "A"])
                writer.writerow(["ENSG00000000001", "1"])
            with self.assertRaises(RNAMeasurementError):
                audit_matrix(path, expected_rows=1, expected_samples=1)


if __name__ == "__main__":
    unittest.main()

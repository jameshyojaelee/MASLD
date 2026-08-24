from __future__ import annotations

import gzip
from pathlib import Path
import tempfile
import unittest

from scripts.audit_gse267145_qc_rights import matrix_qc, robust_z


class GSE267145QCRightsTests(unittest.TestCase):
    def test_rna_header_feature_id_is_not_counted_twice(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "rna.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("ensembl_gene_id\tP1\tP2\n")
                for index in range(43_285):
                    handle.write(f"ENSG{index:011d}\t1\t2\n")
            observed = matrix_qc(path, {"P1", "P2"}, h3=False)
        self.assertEqual(observed["P1"]["feature_rows"], 43_285)
        self.assertEqual(observed["P2"]["library_sum"], 86_570.0)

    def test_robust_z_uses_median_absolute_deviation(self) -> None:
        values = {"a": 1.0, "b": 2.0, "c": 3.0}
        observed = robust_z(values)
        self.assertAlmostEqual(observed["b"], 0.0)
        self.assertLess(observed["a"], 0.0)
        self.assertGreater(observed["c"], 0.0)

    def test_constant_values_have_zero_diagnostic_score(self) -> None:
        self.assertEqual(robust_z({"a": 1.0, "b": 1.0}), {"a": 0.0, "b": 0.0})


if __name__ == "__main__":
    unittest.main()

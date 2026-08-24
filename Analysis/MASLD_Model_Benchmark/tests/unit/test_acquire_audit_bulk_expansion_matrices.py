from __future__ import annotations

import unittest

from scripts.acquire_audit_bulk_expansion_matrices import https_url, numeric_profile


class BulkExpansionMatrixAuditTests(unittest.TestCase):
    def test_numeric_profile_does_not_round_fractional_measurements(self) -> None:
        observed = numeric_profile([0.0, 1.5, 2.0, float("nan")])
        self.assertEqual(observed["numeric_cells"], 4)
        self.assertEqual(observed["finite_numeric_cells"], 3)
        self.assertEqual(observed["nonnegative_integer_cells"], 2)
        self.assertEqual(observed["nonfinite_numeric_cells"], 1)

    def test_ncbi_ftp_is_promoted_to_https(self) -> None:
        self.assertEqual(
            https_url("ftp://ftp.ncbi.nlm.nih.gov/a/matrix.gz"),
            "https://ftp.ncbi.nlm.nih.gov/a/matrix.gz",
        )


if __name__ == "__main__":
    unittest.main()

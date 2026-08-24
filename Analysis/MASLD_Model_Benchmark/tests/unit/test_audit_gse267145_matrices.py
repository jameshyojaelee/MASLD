from __future__ import annotations

import gzip
import unittest

from scripts.audit_gse267145_matrices import MatrixAuditError, audit_matrix


class GSE267145MatrixAuditTests(unittest.TestCase):
    def test_gene_count_matrix(self) -> None:
        payload = gzip.compress(b"gene\tD1\tD2\nG1\t0\t2\nG2\t3\t4\n")
        audit = audit_matrix(payload, expected_candidate_ids={"D1", "D2"})
        self.assertEqual(audit["feature_rows"], 2)
        self.assertEqual(audit["matched_candidate_sample_columns"], 2)
        self.assertTrue(audit["consistent_row_width"])
        self.assertTrue(audit["all_sample_values_nonnegative_integers"])

    def test_peak_count_matrix_with_verbose_titles(self) -> None:
        payload = gzip.compress(
            b"chr\tstart\tend\tD1_nor_H3K27ac\tD2_nash_H3K27ac\n"
            b"chr1\t10\t20\t1\t0\n"
        )
        audit = audit_matrix(payload, expected_candidate_ids={"D1", "D2"})
        self.assertEqual(audit["matched_candidate_ids"], ["D1", "D2"])
        self.assertEqual(audit["duplicate_feature_rows"], 0)

    def test_r_write_table_whitespace_with_implicit_row_names(self) -> None:
        payload = gzip.compress(
            b'"D1_nor_H3K27ac" "D2_nash_H3K27ac"\n'
            b'"chr1:10-20" 1 0\n'
            b'"chr1:30-40" 2 3\n'
        )
        audit = audit_matrix(payload, expected_candidate_ids={"D1", "D2"})
        self.assertEqual(audit["delimiter"], "shell_whitespace")
        self.assertTrue(audit["implicit_row_name_column"])
        self.assertEqual(audit["feature_rows"], 2)
        self.assertEqual(audit["matched_candidate_sample_columns"], 2)

    def test_rejects_negative_or_nonnumeric_counts(self) -> None:
        with self.assertRaises(MatrixAuditError):
            audit_matrix(
                gzip.compress(b"gene\tD1\nG1\t-1\n"),
                expected_candidate_ids={"D1"},
            )
        with self.assertRaises(MatrixAuditError):
            audit_matrix(
                gzip.compress(b"gene\tD1\nG1\tmissing\n"),
                expected_candidate_ids={"D1"},
            )


if __name__ == "__main__":
    unittest.main()

"""Scientific equality checks for bounded growing ridge prediction tables."""
import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location("equivalence", Path(__file__).with_name("check_ridge_prediction_equivalence.py"))
review = importlib.util.module_from_spec(spec)
spec.loader.exec_module(review)


class RidgePredictionEquivalenceTests(unittest.TestCase):
    def rows(self, genes):
        return [dict(gene_id=gene, outer_fold=3, grid=grid, endpoint="rna", method="promoter_only",
                     donor_id=donor, standardized_squared_error="0.123456789", ridge_alpha="100.0")
                for gene in genes for grid in ("original_replay", "wide") for donor in ("A", "B")]

    def table(self, path, rows):
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def pair(self, directory, candidate, original):
        a, b = Path(directory) / "candidate.tsv", Path(directory) / "original.tsv"
        self.table(a, candidate)
        self.table(b, original)
        return a, b

    def test_complete_and_partial_final_gene(self):
        with tempfile.TemporaryDirectory() as directory:
            candidate = self.rows(["g1", "g2", "g3"])
            a, b = self.pair(directory, candidate, candidate[:9])
            receipt = review.compare_file(a, b)
            self.assertEqual(receipt["matched_complete_genes"], 2)
            self.assertEqual(receipt["matched_rows"], 8)
            self.assertEqual(receipt["excluded_incomplete_final_original_gene"], "g3")
            self.assertEqual(receipt["original"]["captured_bytes"], b.stat().st_size)
            with self.assertRaises(ValueError):
                review.compare_file(a, b, original_complete=True)

    def test_pilot_is_complete_original_prefix(self):
        with tempfile.TemporaryDirectory() as directory:
            a, b = self.pair(directory, self.rows(["g1", "g2"]), self.rows(["g1", "g2", "g3"]))
            receipt = review.compare_file(a, b, pilot=True)
            self.assertEqual(receipt["matched_complete_genes"], 2)
            self.assertEqual(receipt["untested_original_rows_after_pilot"], 4)
            self.table(b, self.rows(["g1"]))
            with self.assertRaises(ValueError):
                review.compare_file(a, b, pilot=True)

    def test_missing_changed_reordered_or_duplicate_middle_rejected(self):
        candidate = self.rows(["g1", "g2", "g3"])
        cases = [candidate[:4] + candidate[8:], candidate[:4] + candidate[5:],
                 candidate[:4] + candidate[8:] + candidate[4:8], candidate[:4] + candidate[4:5] + candidate[4:]]
        changed = [dict(row) for row in candidate]
        changed[4]["standardized_squared_error"] = "0.123456788"
        cases.append(changed)
        for original in cases:
            with self.subTest(original=original), tempfile.TemporaryDirectory() as directory:
                a, b = self.pair(directory, candidate, original)
                with self.assertRaises(ValueError):
                    review.compare_file(a, b)

    def test_truncated_eof_row_requires_candidate_prefix(self):
        with tempfile.TemporaryDirectory() as directory:
            candidate = self.rows(["g1", "g2"])
            a, b = self.pair(directory, candidate, candidate)
            lines = b.read_bytes().splitlines(keepends=True)
            b.write_bytes(b"".join(lines[:6]) + lines[6][:18])
            receipt = review.compare_file(a, b)
            self.assertEqual(receipt["matched_complete_genes"], 1)
            self.assertEqual(receipt["excluded_incomplete_final_original_gene"], "g2")
            self.assertEqual(receipt["original"]["truncated_final_row_bytes"], 18)
            b.write_bytes(b"".join(lines[:6]) + b"bad row")
            with self.assertRaises(ValueError):
                review.compare_file(a, b)

    def test_initial_byte_bound_excludes_appended_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "growing.tsv"
            self.table(path, self.rows(["g1"]))
            snapshot = review.PrefixRows(path)
            with path.open("ab") as handle:
                handle.write(b"malformed appended data\n")
            self.assertEqual(len(list(snapshot)), 4)
            receipt = snapshot.receipt()
            self.assertLess(receipt["captured_bytes"], receipt["bytes_at_finish"])

    def test_tuning_exact_both_grids(self):
        rows = [dict(gene_id=gene, grid=grid, method="distance", ridge_alpha="1000.0", inner_mse=".02")
                for gene in ("g1", "g2") for grid in ("original_replay", "wide")]
        with tempfile.TemporaryDirectory() as directory:
            a, b = self.pair(directory, rows, rows)
            receipt = review.compare_file(a, b, original_complete=True, kind="tuning")
            self.assertEqual(receipt["matched_rows"], 4)


if __name__ == "__main__":
    unittest.main()

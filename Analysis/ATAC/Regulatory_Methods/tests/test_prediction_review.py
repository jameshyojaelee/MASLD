"""Small scientific invariants for the independent prediction review."""
import importlib.util
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location("review", Path(__file__).parents[1] / "scripts/review_prediction_results.py")
review = importlib.util.module_from_spec(spec)
spec.loader.exec_module(review)


class PredictionReviewTests(unittest.TestCase):
    def table(self, directory, rows):
        path = Path(directory) / "predictions.tsv"
        review.write(path, rows)
        return path

    def row(self, gene, donor, method, loss, fold=0):
        return dict(outer_fold=fold, gene_id=gene, endpoint="promoter_accessibility", method=method,
                    donor_id=donor, n_metacells=2, cohort="GSE244832", lineage="hepatocyte",
                    standardized_squared_error=loss)

    def test_matching_and_donor_averaging(self):
        # A has two common genes, B only one. The unmatched loss 100 must not enter.
        rows = [self.row("g1", d, m, x) for d, m, x in
                [("A", "intercept", 4), ("A", "distance", 2),
                 ("B", "intercept", 10), ("B", "distance", 1)]]
        rows += [self.row("g2", "A", "intercept", 8), self.row("g2", "A", "distance", 4),
                 self.row("g2", "B", "intercept", 100)]
        with tempfile.TemporaryDirectory() as directory:
            totals, genes, coverage, n = review.reduce_predictions(self.table(directory, rows),
                {"A": (0, 2), "B": (0, 2)}, "GSE244832", "hepatocyte", ["promoter_accessibility"])
        key = "promoter_accessibility", "distance", "intercept"
        self.assertEqual(totals[key]["A"][:2], [6, 2])
        self.assertEqual(totals[key]["B"][:2], [9, 1])
        self.assertEqual((6/2 + 9)/2, 6)  # Donor mean, not gene-row mean 5.
        self.assertEqual(len(genes[key]), 2)
        self.assertEqual(n, 7)

    def test_reject_duplicate_wrong_fold_and_disordered(self):
        a = self.row("g1", "A", "distance", 1)
        cases = [[a, a], [self.row("g1", "A", "distance", 1, fold=1)],
                 [self.row("g2", "A", "distance", 1), a]]
        for rows in cases:
            with self.subTest(rows=rows), tempfile.TemporaryDirectory() as directory:
                with self.assertRaises(ValueError):
                    review.reduce_predictions(self.table(directory, rows), {"A": (0, 2)},
                        "GSE244832", "hepatocyte", ["promoter_accessibility"])

    def test_exact_sign_null_and_numerical_zero(self):
        result = review.inference([1]*5, 42)
        self.assertEqual(result["pvalue"], 1/32)
        self.assertEqual((result["ci_low"], result["ci_high"]), (1, 1))
        self.assertEqual(review.inference([1e-17]*5, 42)["pvalue"], 1)
        self.assertEqual(review.inference([1]*4, 42)["status"], "untestable")

    def test_full_family_holm(self):
        self.assertEqual(len(list(review.planned())), 48)
        adjusted = review.holm([.001, .002] + [1]*46)
        self.assertAlmostEqual(adjusted[0], .048)
        self.assertAlmostEqual(adjusted[1], .094)

    def test_fixed_score_influence(self):
        result = review.score_influence([-1, 1, 9])
        self.assertEqual(result["fraction_positive_donor_gains"], 2/3)
        self.assertEqual(result["min_donor_gain"], -1)
        self.assertEqual(result["max_donor_gain"], 9)
        self.assertEqual(result["min_leave_one_donor_out_mean_gain"], 0)
        self.assertEqual(result["max_leave_one_donor_out_mean_gain"], 5)
        self.assertTrue(review.math.isnan(review.score_influence([1])["min_leave_one_donor_out_mean_gain"]))

    def test_partial_and_invalidated_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "completion.json").write_text('{"status":"running"}')
            with self.assertRaises(ValueError):
                review.review(root)
            (root / "completion.json").write_text('{"status":"completed"}')
            (root / "INVALIDATED.json").write_text('{}')
            with self.assertRaises(ValueError):
                review.review(root)


if __name__ == "__main__":
    unittest.main()

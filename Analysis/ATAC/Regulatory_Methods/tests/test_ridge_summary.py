"""Scientific invariants for the separate wider-ridge donor comparison."""
import importlib.util
import json
import sys
import tempfile
import unittest
from collections import Counter, defaultdict
from pathlib import Path

SCRIPTS = Path(__file__).parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location("ridge_summary", SCRIPTS / "summarize_ridge_sensitivity.py")
review = importlib.util.module_from_spec(spec)
spec.loader.exec_module(review)


class RidgeSummaryTests(unittest.TestCase):
    def row(self, gene, donor, grid, method, loss, alpha=100):
        return dict(outer_fold=0, gene_id=gene, endpoint="rna", grid=grid, method=method,
                    donor_id=donor, n_metacells=2 if donor == "A" else 3,
                    observed=1., predicted=.5, squared_error=loss,
                    standardized_squared_error=loss, ridge_alpha=alpha,
                    target_unit="log1p_raw_RNA_CPM", cohort="GSE296875", lineage="hepatocyte")

    def rows(self):
        result = []
        for gene in ("g1", "g2"):
            for grid in review.GRIDS:
                for method in ("promoter_only", "distance") if gene == "g1" else ("promoter_only",):
                    for donor in ("A", "B"):
                        loss = {("promoter_only", "A"): 4., ("promoter_only", "B"): 10.,
                                ("distance", "A"): 2., ("distance", "B"): 1.}[method, donor]
                        if grid == "wide":
                            loss -= .5
                        result.append(self.row(gene, donor, grid, method, loss))
        return result

    def reduce(self, rows):
        totals = defaultdict(lambda: defaultdict(lambda: [0., 0, 0., 0.]))
        genes, model_genes = defaultdict(set), defaultdict(set)
        alpha, coverage = Counter(), Counter()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "predictions.tsv"
            review.write(path, rows)
            review.reduce_fold(path, "hepatocyte", 0, {"A": (0, 2), "B": (0, 3)},
                               totals, genes, alpha, coverage, model_genes)
        return totals, genes, alpha, coverage

    def test_matching_donor_unit_and_unique_model_alpha(self):
        totals, genes, alpha, coverage = self.reduce(self.rows())
        # g2 lacks distal predictions and must not affect distal-vs-promoter loss.
        contrast = "wide", "distance"
        self.assertEqual(totals[contrast]["A"][:2], [2., 1])
        self.assertEqual(totals[contrast]["B"][:2], [9., 1])
        self.assertEqual(genes[contrast], {"g1"})
        self.assertEqual(alpha["wide", "promoter_only", 100.], 2)
        self.assertEqual(coverage["wide", "promoter_only"], 4)
        # Promoter grid change averages its two genes within each donor.
        self.assertEqual(totals["wide_vs_original", "promoter_only"]["A"][:2], [1., 2])
        row, donor = review.aggregate_row("hepatocyte", contrast, totals[contrast], genes[contrast], 42, 2000)
        self.assertEqual(row["candidate_standardized_mse"], 1.)
        self.assertEqual(row["baseline_standardized_mse"], 6.5)
        self.assertEqual(len(donor), 2)

    def test_coverage_fold_duplicate_and_tuning_rejected(self):
        cases = []
        rows = self.rows()
        cases.append(rows[:2] + rows[3:])  # One original/wide key missing.
        cases.append([dict(row, outer_fold=1) if i == 0 else row for i, row in enumerate(rows)])
        cases.append(rows[:1] + rows)  # Duplicate prediction key.
        cases.append([dict(row, ridge_alpha=10) if i == 0 else row for i, row in enumerate(rows)])
        cases.append([dict(row, ridge_alpha=123) if i == 0 else row for i, row in enumerate(rows)])
        cases.append([dict(row, observed=2) if i == 0 else row for i, row in enumerate(rows)])
        for rows in cases:
            with self.subTest(rows=rows[:2]), self.assertRaises(ValueError):
                self.reduce(rows)

    def test_holm_family_and_direction(self):
        rows = [dict(pvalue=.001, ci_low=.01), dict(pvalue=.005, ci_low=-.01)]
        rows.extend(dict(pvalue=1, ci_low=-1) for _ in range(4))
        review.adjusted(rows, "six_primary")
        self.assertEqual(rows[0]["holm_pvalue"], .006)
        self.assertTrue(rows[0]["predictive_improvement"])
        self.assertFalse(rows[1]["predictive_improvement"])
        self.assertTrue(all(row["family_size"] == 6 for row in rows))

    def test_reference_preserves_original_holm_family(self):
        originals, rows = [], []
        for lineage in review.LINEAGES:
            for method in review.METHODS:
                original = dict(cohort="GSE296875", lineage=lineage, endpoint="rna", method=method,
                                baseline="promoter_only", n_donors=39 if lineage == "hepatocyte" else 18,
                                n_common_predictions=100, n_common_genes=10, improvement=.02,
                                ci_low=.01, ci_high=.03, pvalue=.001, holm_pvalue=.048,
                                family="original_48", predictive_improvement=True)
                originals.append(original)
                rows.append({key: value for key, value in original.items()
                             if key not in {"holm_pvalue", "family", "predictive_improvement"}})
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reference.tsv"
            review.write(path, originals)
            review.check_reference(rows, path)
            self.assertTrue(all(row["family_size"] == 48 and row["holm_pvalue"] == .048 for row in rows))
            rows[0]["n_common_genes"] = 9
            with self.assertRaises(ValueError):
                review.check_reference(rows, path)

    def test_technical_or_changed_grid_receipt_rejected(self):
        valid = dict(status="completed", technical_check=False, max_genes=0, endpoint="rna",
                     original_grid=list(review.GRIDS["original_replay"]), wide_grid=list(review.GRIDS["wide"]))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            receipt = root / "completion.json"
            receipt.write_text(json.dumps(valid))
            self.assertEqual(review.completed(root)["status"], "completed")
            for patch in [dict(status="failed"), dict(technical_check=True), dict(max_genes=10),
                          dict(endpoint="promoter_accessibility"), dict(wide_grid=[.1, 1, 1000])]:
                receipt.write_text(json.dumps({**valid, **patch}))
                with self.assertRaises(ValueError):
                    review.completed(root)


if __name__ == "__main__":
    unittest.main()

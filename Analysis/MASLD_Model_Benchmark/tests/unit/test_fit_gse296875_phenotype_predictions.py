"""Prove the fold firewall and the imputation semantics of the fitting step."""

from __future__ import annotations

import ast
import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts/fit_gse296875_phenotype_predictions.py"
_spec = importlib.util.spec_from_file_location("fit_predictions", SCRIPT)
assert _spec is not None and _spec.loader is not None
fitter = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fitter)


class ProcessSeparationTests(unittest.TestCase):
    def test_the_fitting_process_cannot_import_the_metric_module(self) -> None:
        tree = ast.parse(SCRIPT.read_text())
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module.split(".")[0])
        self.assertNotIn("masld_bench", imported)

    def test_no_metric_is_computed_here(self) -> None:
        source = SCRIPT.read_text()
        for banned in (
            "average_precision",
            "spearman",
            "roc_auc",
            "auprc",
            "r2_score",
        ):
            self.assertNotIn(banned, source)


class NormalisationTests(unittest.TestCase):
    def test_log_cpm_is_per_unit_and_fits_no_parameter(self) -> None:
        counts = np.array([[1, 1, 2], [10, 10, 20]], dtype=np.uint64)
        transformed = fitter.log_cpm(counts)
        # Two units with identical composition and different depth must land on
        # the same normalised profile.
        np.testing.assert_allclose(transformed[0], transformed[1])

    def test_an_empty_unit_becomes_typed_missing_never_zero(self) -> None:
        """Donor 733 has no B cells at all, so that unit has no library."""

        counts = np.array([[1, 1, 2], [0, 0, 0]], dtype=np.uint64)
        transformed = fitter.log_cpm(counts)
        self.assertTrue(np.isfinite(transformed[0]).all())
        self.assertTrue(np.isnan(transformed[1]).all())
        self.assertFalse((transformed[1] == 0.0).any())

    def test_a_structurally_missing_unit_cannot_reach_a_fit(self) -> None:
        expression = np.zeros((6, 12), dtype=np.float64) + 3.0
        expression[0] = np.nan
        observed = np.ones(6, dtype=bool)
        with self.assertRaises(fitter.FitError):
            fitter.molecular_block(
                expression, observed, np.arange(5), np.arange(6)
            )

    def test_a_structurally_missing_unit_is_safe_once_masked(self) -> None:
        rng = np.random.default_rng(3)
        expression = rng.normal(size=(6, 12)) + 5.0
        expression[0] = np.nan
        observed = np.ones(6, dtype=bool)
        observed[0] = False  # masked, so it is imputed and never fitted on
        block = fitter.molecular_block(
            expression, observed, np.arange(5), np.arange(6)
        )
        self.assertTrue(np.isfinite(block).all())
        self.assertEqual(block[0, -1], 0.0)


class FoldFirewallTests(unittest.TestCase):
    def setUp(self) -> None:
        rng = np.random.default_rng(0)
        self.expression = rng.normal(size=(10, 40)) + 5.0
        self.observed = np.ones(10, dtype=bool)
        self.train = np.arange(8)
        self.all_rows = np.arange(10)

    def test_a_test_row_cannot_change_the_fitted_transform(self) -> None:
        first = fitter.molecular_block(
            self.expression, self.observed, self.train, self.all_rows
        )
        contaminated = self.expression.copy()
        contaminated[9] += 1_000.0  # a held-out donor becomes an extreme outlier
        second = fitter.molecular_block(
            contaminated, self.observed, self.train, self.all_rows
        )
        # Training rows keep their projection exactly; only the held-out row moves.
        np.testing.assert_allclose(first[:8], second[:8], atol=1e-9)
        self.assertFalse(np.allclose(first[9], second[9]))

    def test_masked_units_receive_the_training_mean_not_zero(self) -> None:
        observed = self.observed.copy()
        observed[9] = False
        block = fitter.molecular_block(
            self.expression, observed, self.train, self.all_rows
        )
        self.assertEqual(block[9, -1], 0.0)  # indicator says imputed
        self.assertEqual(block[8, -1], 1.0)  # indicator says observed
        # A training-mean fill lands at the centre of the fitted basis, which is
        # the origin after centring, and is not the projection of a zero row.
        np.testing.assert_allclose(block[9, :-1], 0.0, atol=1e-8)
        zero_expression = self.expression.copy()
        zero_expression[9] = 0.0
        zero_block = fitter.molecular_block(
            zero_expression, self.observed, self.train, self.all_rows
        )
        self.assertFalse(np.allclose(zero_block[9, :-1], 0.0, atol=1e-8))

    def test_masked_training_units_do_not_fit_the_transform(self) -> None:
        observed = self.observed.copy()
        observed[0] = False
        clean = fitter.molecular_block(
            self.expression, observed, self.train, self.all_rows
        )
        contaminated = self.expression.copy()
        contaminated[0] += 1_000.0
        dirty = fitter.molecular_block(
            contaminated, observed, self.train, self.all_rows
        )
        np.testing.assert_allclose(clean, dirty, atol=1e-9)


class GridTests(unittest.TestCase):
    def test_grids_match_the_frozen_specification(self) -> None:
        self.assertEqual(
            fitter.ALPHA_GRID, (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)
        )
        self.assertEqual(fitter.C_GRID, (0.001, 0.01, 0.1, 1.0, 10.0))
        self.assertEqual(fitter.MAX_COMPONENTS, 20)
        self.assertEqual(fitter.N_HIGHLY_VARIABLE, 2_000)
        self.assertEqual(fitter.DETECTION_FRACTION, 0.5)

    def test_scopes_are_the_frozen_partition_plus_secondaries(self) -> None:
        self.assertEqual(
            fitter.PRIMARY_LINEAGES,
            ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell"),
        )
        self.assertEqual(fitter.SECONDARY_LINEAGES, ("endothelial_cell", "b_cell"))
        self.assertEqual(fitter.SCOPES[0], "all_lineage")
        self.assertEqual(len(fitter.SCOPES), 8)

    def test_inner_loss_is_log_loss_for_the_binary_endpoint(self) -> None:
        confident_right = fitter.inner_loss("fibrosis", 1.0, 0.99)
        confident_wrong = fitter.inner_loss("fibrosis", 1.0, 0.01)
        self.assertLess(confident_right, confident_wrong)
        self.assertAlmostEqual(fitter.inner_loss("steatosis", 10.0, 7.0), 9.0)


if __name__ == "__main__":
    unittest.main()

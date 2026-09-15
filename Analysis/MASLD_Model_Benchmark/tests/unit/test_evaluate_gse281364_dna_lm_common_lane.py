from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.evaluate_gse281364_dna_lm_common_lane import (
    CommonLaneError,
    activity_delta,
    fit_mlp_outer,
    fit_ridge_outer,
    metric_values,
    _write_tsv,
    verify_frozen_tree,
)


class DNACommonLaneEvaluationTests(unittest.TestCase):
    def _frozen_fixture(self, root: Path) -> str:
        payload = root / "payload.txt"
        payload.write_text("original\n", encoding="utf-8")
        manifest = {
            "artifacts": [
                {
                    "path": "payload.txt",
                    "sha256": sha256(payload.read_bytes()).hexdigest(),
                    "size_bytes": payload.stat().st_size,
                }
            ],
            "metadata": {"fixture": True},
            "schema_version": "masld-bench-artifacts-v1",
        }
        manifest_path = root / "ARTIFACTS.json"
        manifest_path.write_text(
            json.dumps(manifest, separators=(",", ":"), sort_keys=True) + "\n",
            encoding="utf-8",
        )
        digest = sha256(manifest_path.read_bytes()).hexdigest()
        (root / "COMPLETE").write_text(
            json.dumps(
                {
                    "artifact_count": 1,
                    "manifest_sha256": digest,
                    "schema_version": "masld-bench-complete-v1",
                },
                separators=(",", ":"),
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        return digest

    def test_frozen_tree_verifier_rejects_tamper_and_missing_complete(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            digest = self._frozen_fixture(root)
            self.assertEqual(verify_frozen_tree(root, digest)["metadata"], {"fixture": True})
            (root / "payload.txt").write_text("tampered\n", encoding="utf-8")
            with self.assertRaises(CommonLaneError):
                verify_frozen_tree(root, digest)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            digest = self._frozen_fixture(root)
            (root / "COMPLETE").unlink()
            with self.assertRaises(CommonLaneError):
                verify_frozen_tree(root, digest)

    def test_activity_delta_is_signed_alt_minus_ref(self) -> None:
        self.assertAlmostEqual(activity_delta(9, 19, 9, 39, pseudocount=1.0), 1.0)
        self.assertAlmostEqual(activity_delta(10, 10, 10, 10), 0.0)

    def test_ridge_fit_cannot_read_held_out_outcomes(self) -> None:
        rng = np.random.default_rng(7)
        features = rng.normal(size=(50, 12))
        outcomes = features[:, 0] + rng.normal(scale=0.1, size=50)
        folds = np.repeat(np.arange(5), 10)
        first, receipt, _ = fit_ridge_outer(features, outcomes, folds, 2, alphas=(0.1, 1.0))
        modified = outcomes.copy()
        modified[folds == 2] += 1.0e6
        second, _, _ = fit_ridge_outer(features, modified, folds, 2, alphas=(0.1, 1.0))
        np.testing.assert_allclose(first, second)
        self.assertNotIn(2, receipt["outer_training_folds"])
        self.assertFalse(receipt["held_out_outcomes_read_during_fit"])

    def test_ridge_grid_includes_exact_intercept_only_limit(self) -> None:
        rng = np.random.default_rng(9)
        features = rng.normal(size=(50, 12))
        outcomes = rng.normal(size=50)
        folds = np.repeat(np.arange(5), 10)
        prediction, receipt, state = fit_ridge_outer(
            features, outcomes, folds, 3, alphas=(float("inf"),)
        )
        expected = outcomes[folds != 3].mean()
        np.testing.assert_allclose(prediction, expected)
        np.testing.assert_array_equal(state["coefficient"], np.zeros(12))
        self.assertTrue(np.isinf(receipt["selected_alpha"]))

    def test_mlp_fit_cannot_read_held_out_outcomes(self) -> None:
        rng = np.random.default_rng(11)
        features = rng.normal(size=(50, 8))
        outcomes = features[:, 0] + rng.normal(scale=0.1, size=50)
        folds = np.repeat(np.arange(5), 10)
        first, receipt, _ = fit_mlp_outer(
            features, outcomes, folds, 1, seed=13, hidden_width=4, max_epochs=4, patience=2
        )
        modified = outcomes.copy()
        modified[folds == 1] -= 1.0e6
        second, _, _ = fit_mlp_outer(
            features, modified, folds, 1, seed=13, hidden_width=4, max_epochs=4, patience=2
        )
        np.testing.assert_allclose(first, second)
        self.assertNotIn(1, receipt["outer_training_folds"])
        self.assertFalse(receipt["held_out_outcomes_read_during_fit"])

    def test_constant_prediction_has_no_correlation_or_calibration(self) -> None:
        result = metric_values(np.arange(10, dtype=float), np.zeros(10, dtype=float))
        self.assertIsNone(result["pearson"])
        self.assertIsNone(result["spearman"])
        self.assertIsNone(result["calibration_slope"])
        self.assertGreater(float(result["rmse"]), 0)

    def test_tsv_writer_unions_heterogeneous_selection_fields(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "selection.tsv"
            _write_tsv(path, [{"common": "a", "ridge": 1}, {"common": "b", "mlp": 2}])
            lines = path.read_text(encoding="utf-8").splitlines()
            self.assertEqual(lines[0], "common\tridge\tmlp")
            self.assertEqual(lines[1], "a\t1\t")
            self.assertEqual(lines[2], "b\t\t2")


if __name__ == "__main__":
    unittest.main()

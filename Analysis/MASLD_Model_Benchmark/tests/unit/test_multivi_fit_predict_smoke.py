from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "multivi_fit_predict_smoke.py"
SPEC = importlib.util.spec_from_file_location("multivi_fit_predict_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
fit_predict = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fit_predict)


class MultiVIFitPredictTests(unittest.TestCase):
    def test_explicit_missing_state_matches_zero_runtime_mask(self) -> None:
        observed = fit_predict.validate_modality_states([3, 1], ["observed"] * 2)
        missing = fit_predict.validate_modality_states(
            [0, 0], ["structurally_missing"] * 2
        )
        self.assertEqual(observed, (True, True))
        self.assertEqual(missing, (False, False))

    def test_missing_state_with_values_fails_closed(self) -> None:
        with self.assertRaisesRegex(
            fit_predict.MultiVIFitPredictError, "missing modality"
        ):
            fit_predict.validate_modality_states([1], ["structurally_missing"])

    def test_profiles_are_positive_compositions(self) -> None:
        profiles = fit_predict.normalize_profiles(
            np.asarray([[0.0, 2.0, 1.0], [5.0, 0.0, 0.0]])
        )
        self.assertTrue(np.all(profiles > 0))
        np.testing.assert_allclose(profiles.sum(axis=1), 1.0, rtol=0, atol=1e-12)

    def test_join_hash_is_namespaced_and_deterministic(self) -> None:
        first = fit_predict.join_hash("namespace", "row", "donor\0lineage\0peak")
        self.assertEqual(first, fit_predict.join_hash("namespace", "row", "donor\0lineage\0peak"))
        self.assertNotEqual(first, fit_predict.join_hash("namespace", "unit", "donor\0lineage\0peak"))


if __name__ == "__main__":
    unittest.main()

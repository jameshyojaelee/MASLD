from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts" / "evaluate_multivi_smoke.py"
SPEC = importlib.util.spec_from_file_location("multivi_evaluator_tested", MODULE)
assert SPEC is not None and SPEC.loader is not None
evaluator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(evaluator)


class MultiVIEvaluatorTests(unittest.TestCase):
    def test_identical_multinomial_profile_has_zero_deviance(self) -> None:
        observed = np.asarray([3.0, 1.0, 0.0, 2.0])
        self.assertAlmostEqual(
            evaluator.deviance_per_insertion(observed, observed + 1.0e-8),
            0.0,
            places=6,
        )

    def test_normalized_profile_is_positive(self) -> None:
        value = evaluator.normalize_profile([0.0, 2.0, 1.0])
        self.assertTrue(np.all(value > 0))
        self.assertAlmostEqual(float(value.sum()), 1.0)

    def test_deviance_rejects_zero_prediction_mass(self) -> None:
        with self.assertRaisesRegex(
            evaluator.MultiVIEvaluationError, "multinomial mass"
        ):
            evaluator.deviance_per_insertion([1.0, 0.0], [0.0, 0.0])

    def test_join_hash_is_namespaced(self) -> None:
        row = evaluator.join_hash("namespace", "row", "identifier")
        self.assertNotEqual(row, evaluator.join_hash("namespace", "unit", "identifier"))


if __name__ == "__main__":
    unittest.main()

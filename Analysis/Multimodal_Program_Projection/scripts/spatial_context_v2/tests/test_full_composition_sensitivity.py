#!/usr/bin/env python3
"""Lightweight fixtures for the all-cell-type composition sensitivity."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
PRODUCER = SCRIPT_DIR / "31_run_full_composition_sensitivity.py"
spec = importlib.util.spec_from_file_location("full_composition_producer_test", PRODUCER)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {PRODUCER}")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class FullCompositionSensitivitySemantics(unittest.TestCase):
    def test_design_standardization_matches_pinned_convention(self) -> None:
        design = module.standardize_design(
            [
                np.ones(4),
                np.array([1.0, 2.0, 3.0, 4.0]),
                np.repeat(7.0, 4),
            ]
        )
        np.testing.assert_array_equal(design[:, 0], np.ones(4))
        self.assertAlmostEqual(float(design[:, 1].mean()), 0.0)
        self.assertAlmostEqual(float(design[:, 1].std(ddof=1)), 1.0)
        np.testing.assert_array_equal(design[:, 2], np.zeros(4))

    def test_override_uses_all_factors_and_two_qc_covariates(self) -> None:
        index = pd.Index(["a", "b", "c"])
        obs = pd.DataFrame(index=index)
        base = np.column_stack(
            [np.ones(3), np.arange(18 * 3, dtype=float).reshape(3, 18)]
        )
        cached = {"fixture": {"obs": obs, "design": base}}
        engine = SimpleNamespace()
        module.install_design_override(engine, cached, "fixture")
        observed = engine.design_matrix(obs.copy(), np.zeros(3))
        self.assertEqual(observed.shape, (3, module.EXPECTED_FACTOR_COUNT + 3))
        np.testing.assert_array_equal(observed, base)
        zonation = engine.design_matrix(obs.copy(), np.zeros(3), zonation=np.arange(3.0))
        self.assertEqual(zonation.shape, (3, module.EXPECTED_FACTOR_COUNT + 4))

    def test_override_fails_on_observation_reordering(self) -> None:
        obs = pd.DataFrame(index=pd.Index(["a", "b", "c"]))
        cached = {"fixture": {"obs": obs, "design": np.ones((3, 19))}}
        engine = SimpleNamespace()
        module.install_design_override(engine, cached, "fixture")
        with self.assertRaisesRegex(RuntimeError, "observation order drifted"):
            engine.design_matrix(obs.iloc[::-1].copy(), np.zeros(3))


if __name__ == "__main__":
    unittest.main(verbosity=2)

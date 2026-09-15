#!/usr/bin/env python3
"""Unit checks for the outcome-blind SGD elastic-net convergence rescue."""

from __future__ import annotations

from pathlib import Path
import unittest
from unittest import mock

import numpy as np

from scripts import fit_predict_cell_baselines_study_50000 as baseline
from scripts import diagnose_sgd_elastic_net_rescue_study_50000 as diagnostic


ROOT = Path(__file__).parents[2]
WRAPPER = ROOT / "slurm/diagnose_sgd_elastic_net_rescue_study_50000.sbatch"


def prepared_fixture() -> list[dict[str, object]]:
    rng = np.random.default_rng(9)
    result = []
    for inner_fold in range(5):
        fitting_y = np.tile(np.arange(5), 8)
        validation_y = np.tile(np.arange(5), 2)
        fitting = np.zeros(50, dtype=bool)
        validation = np.zeros(50, dtype=bool)
        fitting[:40] = True
        validation[40:] = True
        result.append(
            {
                "inner_fold": inner_fold,
                "fitting": fitting,
                "validation": validation,
                "fitting_x": rng.normal(size=(40, 4)),
                "validation_x": rng.normal(size=(10, 4)),
                "fitting_weights": np.ones(40),
                "validation_weights": np.ones(10),
                "fitting_y": fitting_y,
                "validation_y": validation_y,
            }
        )
    return result


class SGDElasticNetRescueTests(unittest.TestCase):
    def test_weighted_sgd_elastic_net_converges_and_normalizes_probabilities(self) -> None:
        rng = np.random.default_rng(13)
        y = np.tile(np.arange(5), 80)
        x = rng.normal(size=(len(y), 6)) + y[:, None] * 0.25
        weights = np.linspace(0.5, 20.0, len(y))
        model, state = baseline.fit_sgd_elastic_net(
            x,
            y,
            weights,
            alpha=1.0e-3,
            l1_ratio=0.5,
            seed=17,
        )
        self.assertTrue(state["converged_before_ceiling"])
        probabilities = model.predict_proba(x[:7])
        np.testing.assert_allclose(probabilities.sum(axis=1), 1.0)

    def test_selection_fails_closed_when_candidates_do_not_converge(self) -> None:
        prepared = prepared_fixture()
        y = np.tile(np.arange(5), 10)
        fake_model = mock.Mock()
        failed = {
            "converged_before_ceiling": False,
            "iterations": baseline.ELASTIC_NET_MAX_ITER,
        }
        with mock.patch.object(
            baseline, "fit_sgd_elastic_net", return_value=(fake_model, failed)
        ):
            with self.assertRaisesRegex(
                baseline.BaselineScreenError, "too few converged"
            ):
                baseline.select_sgd_elastic_net(prepared, y, fold=0, seed=1)

    def test_frozen_diagnostic_hashes_and_wrapper_firewall(self) -> None:
        self.assertEqual(
            diagnostic.DIAGNOSTIC_504_SHA256,
            "bbccba9f03aaa9736a9ac9146e25d00eb7c9fa9b6b2a06c11ce62e4c6b26306e",
        )
        self.assertEqual(
            diagnostic.DIAGNOSTIC_505_SHA256,
            "8f00036292a60a1712e4da0f5dd8f9e9c4e8a47541e4af8a15778d145d2124c0",
        )
        text = WRAPPER.read_text(encoding="utf-8")
        header = "\n".join(
            line for line in text.splitlines() if line.startswith("#SBATCH")
        )
        self.assertIn("--partition=cpu", header)
        self.assertIn("--account=nslab", header)
        self.assertIn("--qos=nslab", header)
        self.assertIn("--cpus-per-task=16", header)
        self.assertIn("--mem=128G", header)
        self.assertNotIn("--array", header)
        self.assertNotIn("innovation", text.lower())
        self.assertNotIn("predict_test", text)
        self.assertIn(
            "8a7dfd2a15a635d5ef5b25fa8ee68355e6c6cfd047ec464cdc1dea83b7201298",
            text,
        )
        self.assertIn(
            "672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509",
            text,
        )
        self.assertIn(
            "10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680",
            text,
        )
        self.assertFalse(
            any(line.lstrip().startswith("sbatch ") for line in text.splitlines())
        )
        source = (
            ROOT / "scripts/diagnose_sgd_elastic_net_rescue_study_50000.py"
        ).read_text(encoding="utf-8")
        self.assertIn("np.flatnonzero(outer != arguments.outer_fold)", source)
        self.assertNotIn("np.flatnonzero(outer == arguments.outer_fold)", source)
        self.assertIn('"outer_test_partition_predicted": False', source)
        self.assertIn('"outer_test_metrics_calculated": False', source)


if __name__ == "__main__":
    unittest.main()

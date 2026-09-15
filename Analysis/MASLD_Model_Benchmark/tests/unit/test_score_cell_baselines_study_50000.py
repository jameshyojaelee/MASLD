#!/usr/bin/env python3
"""Unit checks for the independent 50,000-cell evaluation-only scorer."""

from __future__ import annotations

import ast
from pathlib import Path
import unittest

import numpy as np

from scripts import score_cell_baselines_study_50000 as scorer


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts/score_cell_baselines_study_50000.py"
WRAPPER = ROOT / "slurm/score_cell_baselines_study_50000.sbatch"


class CellBaselineStudy50000ScorerTests(unittest.TestCase):
    def test_perfect_calibration_and_scores(self) -> None:
        truth = np.arange(5, dtype=np.int8)
        probabilities = np.eye(5, dtype=np.float64)
        donors = np.asarray([f"d{index}" for index in range(5)])
        value = scorer.score_subset(truth, donors, probabilities)
        self.assertEqual(value["macro_f1"], 1.0)
        self.assertEqual(value["multiclass_brier"], 0.0)
        self.assertEqual(value["calibration"]["top_label_ece"], 0.0)

    def test_bootstrap_is_deterministic_and_preserves_study_counts(self) -> None:
        studies = np.asarray(["a", "a", "b", "b"])
        present = np.asarray(
            [
                [True, True, False, False, False],
                [False, False, True, True, True],
                [True, False, True, False, True],
                [False, True, False, True, False],
            ]
        )
        first = scorer.build_multiplicities(
            studies, present, n_resamples=100, seed=17
        )
        second = scorer.build_multiplicities(
            studies, present, n_resamples=100, seed=17
        )
        np.testing.assert_array_equal(first, second)
        np.testing.assert_array_equal(first[:, :2].sum(axis=1), np.full(100, 2))
        np.testing.assert_array_equal(first[:, 2:].sum(axis=1), np.full(100, 2))

    def test_bootstrap_accepted_draws_have_exact_multiplicities(self) -> None:
        studies = np.asarray(["a", "a"])
        present = np.ones((2, 5), dtype=bool)
        observed = scorer.build_multiplicities(
            studies, present, n_resamples=4, seed=17
        )
        expected = np.asarray([[0, 2], [2, 0], [1, 1], [1, 1]], dtype=np.int16)
        np.testing.assert_array_equal(observed, expected)

    def test_strongest_rule_is_deterministic(self) -> None:
        values = {
            model_id: {
                "donor_class_balanced_macro_f1": 0.5,
                "study_balanced_macro_f1": 0.4,
                "multiclass_brier": 0.3,
            }
            for model_id in scorer.MODEL_IDS
        }
        values["hvg_pca_logistic"]["donor_class_balanced_macro_f1"] = 0.6
        self.assertEqual(scorer.choose_strongest(values), "hvg_pca_logistic")

    def test_sufficient_statistics_rederive_perfect_endpoint(self) -> None:
        truth = np.arange(5, dtype=np.int8)
        probabilities = np.eye(5, dtype=np.float64)
        stats = scorer.sufficient_statistics(
            truth,
            truth,
            probabilities,
            np.arange(5, dtype=np.int16),
            donors=5,
        )
        value = scorer.endpoints_from_multiplicities(
            stats,
            np.ones((1, 5), dtype=np.int16),
            np.arange(5),
            np.ones(5, dtype=bool),
        )
        self.assertEqual(float(value["macro_f1"][0]), 1.0)
        self.assertEqual(float(value["brier"][0]), 0.0)
        self.assertEqual(float(value["top_label_ece"][0]), 0.0)
        np.testing.assert_array_equal(value["class_f1"], np.ones((1, 5)))

    def test_production_contract_is_evaluation_only(self) -> None:
        self.assertEqual(scorer.BOOTSTRAP_RESAMPLES, 10_000)
        self.assertEqual(scorer.SCREEN_SEEDS, (20260824, 20260825, 20260826))
        tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
        called_attributes = {
            node.func.attr
            for node in ast.walk(tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
        }
        self.assertNotIn("fit", called_attributes)
        self.assertNotIn("fit_transform", called_attributes)
        text = SCRIPT.read_text(encoding="utf-8")
        self.assertNotIn("sklearn", text)
        self.assertNotIn("anndata", text)
        self.assertIn('"models_fit": False', text)
        self.assertIn('"histology_read": False', text)
        self.assertIn('"sealed_outcomes_read": False', text)
        self.assertIn('"seed_metrics.tsv"', text)
        self.assertIn("ensemble prediction differs from frozen seed mean", text)

    def test_wrapper_is_generic_cpu_nslab_and_not_self_submitting(self) -> None:
        text = WRAPPER.read_text(encoding="utf-8")
        header = "\n".join(
            line for line in text.splitlines() if line.startswith("#SBATCH")
        )
        self.assertIn("--job-name=model-cpu-score-502", header)
        self.assertIn("--partition=cpu", header)
        self.assertIn("--account=nslab", header)
        self.assertIn("--qos=nslab", header)
        self.assertNotIn("--array", header)
        self.assertNotIn("innovation", text.lower())
        self.assertFalse(
            any(line.lstrip().startswith("sbatch ") for line in text.splitlines())
        )
        self.assertIn('PREDICTIONS_ROOT:?PREDICTIONS_ROOT is required', text)
        self.assertIn('PREDICTIONS_SHA256:?PREDICTIONS_SHA256 is required', text)


if __name__ == "__main__":
    unittest.main()

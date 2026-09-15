#!/usr/bin/env python3
"""Synthetic contract tests for the independent GSE274114 evaluator."""

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts import evaluate_gse274114_bundled_baselines as evaluator


ROOT = Path(__file__).resolve().parents[2]


class GSE274114IndependentEvaluatorTests(unittest.TestCase):
    def test_registered_metrics_have_expected_perfect_values(self) -> None:
        truth = np.asarray([0, 0, 1, 1], dtype=np.int64)
        probability = np.asarray([0.1, 0.2, 0.8, 0.9], dtype=np.float64)
        metrics = evaluator._registered_metrics(truth, probability)
        self.assertEqual(tuple(metrics), evaluator.REGISTERED_METRICS)
        self.assertEqual(metrics["participant_macro_f1"], 1.0)
        self.assertEqual(metrics["participant_auprc"], 1.0)
        self.assertEqual(metrics["participant_auroc"], 1.0)
        self.assertGreater(metrics["participant_brier"], 0.0)
        self.assertGreater(metrics["participant_log_loss"], 0.0)

    def test_bootstrap_is_registered_participant_stratified_and_deterministic(self) -> None:
        truth = np.asarray([0, 0, 0, 1, 1, 1], dtype=np.int64)
        probability = np.asarray([0.1, 0.4, 0.2, 0.6, 0.9, 0.8], dtype=np.float64)
        first = evaluator._participant_stratified_bootstrap(
            truth, probability, replicates=10_000
        )
        second = evaluator._participant_stratified_bootstrap(
            truth, probability, replicates=10_000
        )
        self.assertEqual(first, second)
        self.assertEqual(tuple(first), evaluator.REGISTERED_METRICS)
        with self.assertRaises(evaluator.GSE274114IndependentEvaluatorError):
            evaluator._participant_stratified_bootstrap(
                truth, probability, replicates=999
            )

    def test_task_contracts_are_platform_specific_and_prohibit_champion_claims(self) -> None:
        for task_id, task in evaluator.TASKS.items():
            with self.subTest(task_id=task_id):
                self.assertIn("within-", task["interpretation"])
                self.assertIn("champion", task["prohibited_claim"])
                self.assertNotIn("cross-platform utility", task["interpretation"])
                self.assertNotIn("MASLD diagnostic", task["interpretation"])

    def test_prediction_failure_occurs_before_missing_label_paths_are_opened(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "must_not_exist"
            with self.assertRaises(FileNotFoundError):
                evaluator.run_evaluation(
                    campaign_root=root / "missing_campaign",
                    campaign_artifacts_sha256="0" * 64,
                    labels_path=root / "missing_labels.tsv",
                    labels_sha256="1" * 64,
                    folds_path=root / "missing_folds.tsv",
                    folds_sha256="2" * 64,
                    task_spec_paths={
                        task_id: root / f"missing_{task_id}.toml"
                        for task_id in evaluator.TASKS
                    },
                    bootstrap_replicates=10_000,
                    output=output,
                )
            self.assertFalse(output.exists())

    def test_source_does_not_import_or_call_training_implementation(self) -> None:
        source = Path(evaluator.__file__).read_text(encoding="utf-8")
        self.assertNotIn("gse274114_baseline_firewall", source)
        self.assertNotIn("run_gse274114_baseline_task_campaign", source)
        self.assertNotIn("fit_predict", source)
        self.assertNotIn("_predict_baseline", source)
        self.assertNotIn("evaluate_model", source)
        self.assertNotIn("ttest", source)
        self.assertIn('"nominal_or_confirmatory_p_values_computed": False', source)
        self.assertIn('"pooled_cross_platform_scoring_run": False', source)
        self.assertIn('"task_champion_selected": False', source)

    def test_registered_taskspecs_pass_independent_contract_validation(self) -> None:
        paths = {
            "gse274114_hiseq_healthy_vs_hbv": ROOT
            / "config/evaluation/gse274114_hiseq_healthy_vs_hbv_task.toml",
            "gse274114_novaseq_mash_vs_mash_hbv": ROOT
            / "config/evaluation/gse274114_novaseq_mash_vs_mash_hbv_task.toml",
        }
        for task_id, path in paths.items():
            with self.subTest(task_id=task_id):
                evaluator._validate_task_spec(path, task_id)


if __name__ == "__main__":
    unittest.main()

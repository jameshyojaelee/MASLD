#!/usr/bin/env python3
"""Contract tests for the bundled GSE274114 prediction-only campaign."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
import unittest

from masld_bench.registry import load_task_spec
from scripts import gse274114_baseline_firewall as firewall
from scripts import run_gse274114_baseline_task_campaign as campaign


ROOT = Path(__file__).resolve().parents[2]


class GSE274114BaselineTaskCampaignTests(unittest.TestCase):
    def test_each_task_bundles_exactly_sixty_unique_prediction_runs(self) -> None:
        for task_id in firewall.TASKS:
            with self.subTest(task_id=task_id):
                runs = campaign.planned_runs(task_id)
                self.assertEqual(len(runs), 60)
                self.assertEqual(
                    len(
                        {
                            (row["model_id"], row["outer_fold"], row["seed"])
                            for row in runs
                        }
                    ),
                    60,
                )
                self.assertEqual(
                    Counter(row["model_kind"] for row in runs),
                    Counter(
                        {
                            "training_class_prior": 5,
                            "gene_rank_nearest_centroid": 5,
                            "hvg_pca_elastic_net": 25,
                            "hvg_pca_linear_svm": 25,
                        }
                    ),
                )

    def test_five_stochastic_seeds_are_genuinely_distinct(self) -> None:
        self.assertEqual(len(firewall.FIXED_STOCHASTIC_SEEDS), 5)
        self.assertEqual(len(set(firewall.FIXED_STOCHASTIC_SEEDS)), 5)
        self.assertNotIn(
            firewall.DETERMINISTIC_SEED, firewall.FIXED_STOCHASTIC_SEEDS
        )
        runs = campaign.planned_runs("gse274114_hiseq_healthy_vs_hbv")
        for model_kind in firewall.STOCHASTIC_MODELS:
            for outer_fold in range(5):
                self.assertEqual(
                    {
                        row["seed"]
                        for row in runs
                        if row["model_kind"] == model_kind
                        and row["outer_fold"] == outer_fold
                    },
                    set(firewall.FIXED_STOCHASTIC_SEEDS),
                )

    def test_planned_models_equal_each_frozen_taskspec_baseline_roster(self) -> None:
        paths = {
            "gse274114_hiseq_healthy_vs_hbv": ROOT
            / "config/evaluation/gse274114_hiseq_healthy_vs_hbv_task.toml",
            "gse274114_novaseq_mash_vs_mash_hbv": ROOT
            / "config/evaluation/gse274114_novaseq_mash_vs_mash_hbv_task.toml",
        }
        for task_id, path in paths.items():
            task = load_task_spec(path)
            observed = {row["model_id"] for row in campaign.planned_runs(task_id)}
            self.assertEqual(observed, set(task.baseline_model_ids))

    def test_campaign_source_has_no_evaluator_or_metric_call(self) -> None:
        source = Path(campaign.__file__).read_text(encoding="utf-8")
        self.assertNotIn("evaluate_model(", source)
        self.assertNotIn("sklearn.metrics", source)
        self.assertNotIn("joined_evaluator_rows", source)
        self.assertIn('"scoring_or_evaluation_run": False', source)
        self.assertIn('"held_fold_labels_in_views": False', source)


if __name__ == "__main__":
    unittest.main()

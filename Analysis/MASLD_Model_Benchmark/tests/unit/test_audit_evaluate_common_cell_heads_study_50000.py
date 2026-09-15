#!/usr/bin/env python3
"""Focused checks for the 50,000-cell common-head audit and evaluator."""

from __future__ import annotations

import ast
from pathlib import Path
from unittest import mock
import unittest

import numpy as np

from scripts import audit_aggregate_common_cell_heads_study_50000 as aggregator
from scripts import score_common_cell_heads_study_50000 as scorer


ROOT = Path(__file__).parents[2]
AGGREGATOR = ROOT / "scripts/audit_aggregate_common_cell_heads_study_50000.py"
SCORER = ROOT / "scripts/score_common_cell_heads_study_50000.py"
AUDIT_WRAPPER = ROOT / "slurm/audit_aggregate_common_cell_heads_study_50000.sbatch"
SCORE_WRAPPER = ROOT / "slurm/score_common_cell_heads_study_50000.sbatch"


class CommonCellHeadStudy50000Tests(unittest.TestCase):
    def test_probability_contract_accepts_argmax_and_rejects_label_field(self) -> None:
        row = {
            "predicted_class": "cholangiocyte",
            **{
                f"probability::{label}": "1.0" if index == 0 else "0.0"
                for index, label in enumerate(aggregator.ROSTER)
            },
        }
        aggregator.validate_probabilities([row])
        row["predicted_class"] = "immune"
        with self.assertRaisesRegex(
            aggregator.CommonHeadAuditError, "probability contract"
        ):
            aggregator.validate_probabilities([row])

    def test_donor_safe_split_rejects_cross_fold_donor(self) -> None:
        rows = [
            {"row_id": "r1", "donor_id": "d1", "dataset": "s1", "outer_fold": "0"},
            {"row_id": "r2", "donor_id": "d2", "dataset": "s2", "outer_fold": "1"},
        ]
        patches = (
            mock.patch.object(aggregator, "EXPECTED_ROWS", 2),
            mock.patch.object(aggregator, "EXPECTED_DONORS", 2),
            mock.patch.object(aggregator, "EXPECTED_STUDIES", 2),
            mock.patch.object(aggregator, "OUTER_FOLDS", (0, 1)),
        )
        with patches[0], patches[1], patches[2], patches[3]:
            aggregator.validate_donor_safe_split(rows)
            rows[1]["donor_id"] = "d1"
            with self.assertRaisesRegex(
                aggregator.CommonHeadAuditError, "donor- and study-safe"
            ):
                aggregator.validate_donor_safe_split(rows)

    def test_shard_census_is_exact(self) -> None:
        keys = aggregator.expected_shard_keys()
        self.assertEqual(len(keys), 30)
        self.assertEqual({key[0] for key in keys}, set(aggregator.HEADS))
        self.assertEqual({key[1] for key in keys}, set(aggregator.SCREEN_SEEDS))
        self.assertEqual({key[2] for key in keys}, set(aggregator.OUTER_FOLDS))

    def test_seed_ensemble_is_unweighted_probability_mean(self) -> None:
        values = [
            np.tile(np.asarray([[0.6, 0.1, 0.1, 0.1, 0.1]]), (2, 1)),
            np.tile(np.asarray([[0.3, 0.2, 0.2, 0.2, 0.1]]), (2, 1)),
            np.tile(np.asarray([[0.3, 0.1, 0.2, 0.1, 0.3]]), (2, 1)),
        ]
        with mock.patch.object(scorer, "EXPECTED_ROWS", 2):
            observed = scorer.mean_seed_probabilities(values)
        np.testing.assert_allclose(
            observed,
            np.tile(np.asarray([[0.4, 0.13333333333333333, 0.16666666666666666, 0.13333333333333333, 0.16666666666666666]]), (2, 1)),
        )

    def test_evaluator_uses_50k_metric_contract_not_training_or_smoke(self) -> None:
        tree = ast.parse(SCORER.read_text(encoding="utf-8"))
        imports = {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
        } | {
            node.module or ""
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom)
        }
        self.assertIn("scripts", imports)
        text = SCORER.read_text(encoding="utf-8")
        self.assertNotIn("fit_predict_common_cell_heads_study_50000", text)
        self.assertNotIn("cell_state_development", text)
        self.assertIn('"training_adapter_imported": False', text)
        self.assertIn('"head_selection_performed": False', text)
        self.assertIn('"comparison_to_baseline_performed": False', text)
        self.assertIn('"champion_promotion_performed": False', text)
        self.assertIn('"sealed_outcomes_read": False', text)
        self.assertNotIn('"transcriptformer_tf_sapiens"', text)
        self.assertIn("expected_model_id", text)
        self.assertIn("expected_exposure_status", text)

    def test_aggregator_does_not_read_development_labels_or_metrics(self) -> None:
        tree = ast.parse(AGGREGATOR.read_text(encoding="utf-8"))
        called_attributes = {
            node.func.attr
            for node in ast.walk(tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
        }
        text = AGGREGATOR.read_text(encoding="utf-8")
        self.assertNotIn("fit", called_attributes)
        self.assertNotIn("fit_transform", called_attributes)
        self.assertNotIn("selection.tsv", text)
        self.assertIn('"development_labels_read_by_aggregator": []', text)
        self.assertIn('"metrics_calculated": False', text)
        self.assertNotIn('"transcriptformer_tf_sapiens"', text)
        self.assertIn("expected_model_id", text)

    def test_model_identity_is_explicitly_pinned(self) -> None:
        self.assertEqual(
            aggregator.require_model_id("geneformer_v2_316m"),
            "geneformer_v2_316m",
        )
        with self.assertRaisesRegex(aggregator.CommonHeadAuditError, "model_id"):
            aggregator.require_model_id("../spoof")
        self.assertEqual(
            scorer.require_model_identity(
                "geneformer_v2_316m", "target_label_unexposed", True
            ),
            ("geneformer_v2_316m", "target_label_unexposed", True),
        )
        with self.assertRaisesRegex(scorer.CommonHeadScoreError, "contradicts"):
            scorer.require_model_identity("fixture_model", "unknown", True)

    def test_wrappers_are_cpu_nslab_and_not_self_submitting(self) -> None:
        for wrapper, job_name in (
            (AUDIT_WRAPPER, "model-cpu-audit-503"),
            (SCORE_WRAPPER, "model-cpu-score-504"),
        ):
            text = wrapper.read_text(encoding="utf-8")
            header = "\n".join(
                line for line in text.splitlines() if line.startswith("#SBATCH")
            )
            self.assertIn(f"--job-name={job_name}", header)
            self.assertIn("--partition=cpu", header)
            self.assertIn("--account=nslab", header)
            self.assertIn("--qos=nslab", header)
            self.assertNotIn("--array", header)
            self.assertNotIn("innovation", text.lower())
            self.assertFalse(
                any(line.lstrip().startswith("sbatch ") for line in text.splitlines())
            )
        score_text = SCORE_WRAPPER.read_text(encoding="utf-8")
        audit_text = AUDIT_WRAPPER.read_text(encoding="utf-8")
        self.assertIn('BUNDLE_ROOT:?BUNDLE_ROOT is required', audit_text)
        self.assertIn('BUNDLE_SHA256:?BUNDLE_SHA256 is required', audit_text)
        self.assertIn('EXPECTED_MODEL_ID:?EXPECTED_MODEL_ID is required', audit_text)
        self.assertIn('PREDICTIONS_ROOT:?PREDICTIONS_ROOT is required', score_text)
        self.assertIn('PREDICTIONS_SHA256:?PREDICTIONS_SHA256 is required', score_text)
        self.assertIn('EXPECTED_MODEL_ID:?EXPECTED_MODEL_ID is required', score_text)
        self.assertIn(
            'EXPECTED_EXPOSURE_STATUS:?EXPECTED_EXPOSURE_STATUS is required',
            score_text,
        )
        self.assertIn(
            'EXPECTED_SEALED_CHAMPION_ELIGIBLE:?EXPECTED_SEALED_CHAMPION_ELIGIBLE is required',
            score_text,
        )


if __name__ == "__main__":
    unittest.main()

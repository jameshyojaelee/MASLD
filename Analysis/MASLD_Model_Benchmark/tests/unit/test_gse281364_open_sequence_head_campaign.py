from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from scripts import evaluate_gse281364_open_sequence_heads as evaluate
from scripts import fit_gse281364_open_sequence_heads as fit
from scripts import freeze_gse281364_open_sequence_taskspec as freeze


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/gse281364_open_sequence_head_campaign.toml"


class OpenSequenceHeadCampaignTests(unittest.TestCase):
    def test_frozen_config_roster_and_firewalls(self) -> None:
        config = freeze.load_config(CONFIG)
        freeze.validate_config(config)
        candidates = {
            (row["model_id"], head)
            for row in config["candidate"]
            for head in row["head_ids"]
        }
        self.assertEqual(len(candidates), 12)
        self.assertEqual(tuple(config["task_spec"]["fixed_seeds"]), freeze.SEEDS)
        self.assertEqual(config["task_spec"]["biological_donors"], 0)
        self.assertFalse(config["stack_fit_authorized"])
        self.assertFalse(config["residual_correlation_authorized"])
        self.assertFalse(config["conditional_model_build_authorized"])
        self.assertFalse(config["sealed_asset_access_authorized"])
        self.assertFalse(config["completion_firewall"]["mandatory_baselines_complete"])
        self.assertEqual(
            tuple(config["completion_firewall"]["mandatory_baselines_missing"]),
            freeze.MISSING_BASELINES,
        )
        self.assertTrue(config["completion_firewall"]["shortlist_blocked"])
        self.assertTrue(config["completion_firewall"]["complementarity_blocked"])

    def test_taskspec_is_frozen_before_outcome_values_or_fit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "task_spec"
            receipt = freeze.freeze_spec(root=ROOT, config_path=CONFIG, output=output)
            spec = json.loads((output / "task_spec.json").read_text(encoding="utf-8"))
            self.assertEqual(receipt["status"], "pass_prespecified_taskspec")
            self.assertFalse(receipt["outcome_values_read"])
            self.assertFalse(receipt["model_fit"])
            self.assertEqual(spec["uncertainty"]["resampling_unit"], "long_range_block_id")
            self.assertEqual(spec["multiplicity"]["confirmatory_family"], "none_exposed_development_screen")
            self.assertFalse(spec["completion_firewall"]["mandatory_baselines_complete"])
            self.assertTrue(spec["completion_firewall"]["shortlist_blocked"])
            with self.assertRaises(freeze.TaskSpecError):
                freeze.freeze_spec(root=ROOT, config_path=CONFIG, output=output)

    def test_inner_projection_ignores_excluded_blocks(self) -> None:
        rng = np.random.default_rng(31)
        embeddings = rng.normal(size=(50, 4, 64))
        folds = np.tile(np.arange(5), 10)
        training = (folds != 0) & (folds != 1)
        changed = embeddings.copy()
        changed[~training] = rng.normal(loc=1.0e6, scale=1.0e3, size=changed[~training].shape)
        with patch.object(fit, "PROJECTION_WIDTH", 32):
            original = fit.fit_inner_projection(embeddings, training)
            altered = fit.fit_inner_projection(changed, training)
        self.assertEqual(set(original), set(altered))
        for name in original:
            np.testing.assert_array_equal(original[name], altered[name])

    def test_ridge_selection_and_fit_ignore_held_outcomes(self) -> None:
        rng = np.random.default_rng(73)
        folds = np.tile(np.arange(5), 12)
        inner_features = rng.normal(size=(60, 9))
        outer_features = rng.normal(size=(60, 9))
        outcomes = rng.normal(size=60)
        changed = outcomes.copy()
        changed[folds == 2] += 1.0e9
        first = fit.fit_ridge_outer(
            inner_features=inner_features,
            outer_features=outer_features,
            outcomes=outcomes,
            folds=folds,
            held_fold=2,
            alphas=(0.01, 1.0, 100.0),
        )
        second = fit.fit_ridge_outer(
            inner_features=inner_features,
            outer_features=outer_features,
            outcomes=changed,
            folds=folds,
            held_fold=2,
            alphas=(0.01, 1.0, 100.0),
        )
        np.testing.assert_array_equal(first[0], second[0])
        self.assertEqual(first[1], second[1])
        for name in first[2]:
            np.testing.assert_array_equal(first[2][name], second[2][name])

    def test_primary_metric_is_context_macro_fisher_z(self) -> None:
        left = np.asarray([1.0, 2.0, 3.0, 4.0])
        right = np.asarray([4.0, 1.0, 3.0, 2.0])
        first = evaluate.metric_values(left, left)["spearman"]
        second = evaluate.metric_values(left, right)["spearman"]
        observed = evaluate.fisher_macro([first, second])
        expected = float(
            np.tanh(
                np.mean(
                    np.arctanh(
                        np.clip(np.asarray([first, second]), -0.999999, 0.999999)
                    )
                )
            )
        )
        self.assertAlmostEqual(observed, expected)

    def test_block_bootstrap_is_paired_and_deterministic(self) -> None:
        observed = {
            context: np.linspace(0.0, 1.0, 20) + index
            for index, context in enumerate(freeze.CONTEXTS)
        }
        predictions = {
            ("candidate", "head"): {
                context: values + np.linspace(0.01, 0.2, 20)
                for context, values in observed.items()
            }
        }
        blocks = np.asarray([f"block-{index // 2}" for index in range(20)])
        first = evaluate.bootstrap_primary(
            observed=observed,
            predictions=predictions,
            block_ids=blocks,
            resamples=40,
            seed=20260825,
        )
        second = evaluate.bootstrap_primary(
            observed=observed,
            predictions=predictions,
            block_ids=blocks,
            resamples=40,
            seed=20260825,
        )
        np.testing.assert_array_equal(
            first[("candidate", "head")], second[("candidate", "head")]
        )
        self.assertEqual(first[("candidate", "head")].shape, (40,))

    def test_fit_and_evaluator_implementations_are_separate(self) -> None:
        fit_source = (ROOT / "scripts/fit_gse281364_open_sequence_heads.py").read_text(
            encoding="utf-8"
        )
        evaluator_source = (
            ROOT / "scripts/evaluate_gse281364_open_sequence_heads.py"
        ).read_text(encoding="utf-8")
        self.assertNotIn("evaluate_gse281364_open_sequence_heads", fit_source)
        self.assertNotIn("fit_gse281364_open_sequence_heads", evaluator_source)
        self.assertNotIn("pvalue", evaluator_source.lower())


if __name__ == "__main__":
    unittest.main()

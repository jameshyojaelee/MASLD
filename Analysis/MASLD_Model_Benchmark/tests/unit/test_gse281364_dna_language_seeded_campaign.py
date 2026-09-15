from __future__ import annotations

import copy
from pathlib import Path
import unittest

import numpy as np

from scripts.fit_gse281364_dna_language_seeded_heads import (
    DNASeededHeadFitError,
    fit_seeded_dual_projection_fold,
)
from scripts.fit_gse281364_enformer_sei_heads import fit_seeded_outer_fold
from scripts.project_gse281364_long_range_embeddings import MODELS as PROJECTION_MODELS
from scripts.reconcile_gse281364_dna_language_seeded_features import (
    CANDIDATES,
    MODELS,
    PREDICTION_FIELDS,
    load_config,
    validate_config,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/gse281364_dna_language_seeded_head_campaign.json"


class DNASeededCampaignTests(unittest.TestCase):
    @staticmethod
    def synthetic() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        rng = np.random.default_rng(17)
        blocks = np.repeat(np.asarray([f"block-{index:02d}" for index in range(50)]), 2)
        folds = np.repeat(np.arange(50) % 5, 2)
        inner = rng.normal(size=(100, 32))
        outer = inner + rng.normal(scale=0.01, size=inner.shape)
        outcome = inner[:, 5] - 0.3 * inner[:, 11] + rng.normal(scale=0.2, size=100)
        return inner, outer, outcome, folds, blocks

    def fit(self, inner: np.ndarray, outer: np.ndarray, outcome: np.ndarray, folds: np.ndarray, blocks: np.ndarray, seed: int = 1103):
        return fit_seeded_dual_projection_fold(
            inner_features=inner,
            outer_features=outer,
            outcomes=outcome,
            folds=folds,
            block_ids=blocks,
            held_fold=3,
            seed=seed,
            top_k_grid=(4, 16, 32),
            alphas=(0.01, 1.0, 100.0),
        )

    def test_exact_roster_and_projection_support(self) -> None:
        config = load_config(CONFIG)
        self.assertEqual(MODELS, PROJECTION_MODELS)
        self.assertEqual(len(CANDIDATES), 11)
        self.assertEqual([row["model_id"] for row in config["models"]], list(MODELS))
        self.assertTrue(config["models"][-1]["restricted"])
        self.assertFalse(config["models"][-1]["open_champion_eligible_after_external_gates"])
        self.assertEqual(
            config["fit"]["shared_control_prediction_policy"],
            "reuse_bit_exact_frozen_enformer_sei_oof",
        )
        self.assertNotIn("observed", PREDICTION_FIELDS)

    def test_held_outcomes_cannot_change_fit(self) -> None:
        inner, outer, outcome, folds, blocks = self.synthetic()
        prediction, receipt, state = self.fit(inner, outer, outcome, folds, blocks)
        poisoned = outcome.copy()
        poisoned[folds == 3] = 1.0e15
        other_prediction, other_receipt, other_state = self.fit(inner, outer, poisoned, folds, blocks)
        np.testing.assert_array_equal(prediction, other_prediction)
        self.assertEqual(receipt, other_receipt)
        for key in state:
            np.testing.assert_array_equal(state[key], other_state[key])

    def test_held_features_cannot_change_fit_state(self) -> None:
        inner, outer, outcome, folds, blocks = self.synthetic()
        _, receipt, state = self.fit(inner, outer, outcome, folds, blocks)
        poisoned_inner = inner.copy()
        poisoned_outer = outer.copy()
        poisoned_inner[folds == 3] = 1.0e15
        poisoned_outer[folds == 3] = 1.0e15
        _, other_receipt, other_state = self.fit(poisoned_inner, poisoned_outer, outcome, folds, blocks)
        self.assertEqual(receipt, other_receipt)
        for key in state:
            np.testing.assert_array_equal(state[key], other_state[key])

    def test_seeds_are_distinct_block_bootstraps(self) -> None:
        inner, outer, outcome, folds, blocks = self.synthetic()
        first, _, _ = self.fit(inner, outer, outcome, folds, blocks, seed=1103)
        second, _, _ = self.fit(inner, outer, outcome, folds, blocks, seed=2909)
        self.assertFalse(np.array_equal(first, second))

    def test_shared_feature_head_matches_frozen_enformer_sei_draws(self) -> None:
        inner, _, outcome, folds, blocks = self.synthetic()
        prediction, receipt, state = fit_seeded_dual_projection_fold(
            inner_features=inner,
            outer_features=inner,
            outcomes=outcome,
            folds=folds,
            block_ids=blocks,
            held_fold=3,
            seed=1103,
            top_k_grid=(4, 16, 32),
            alphas=(0.01, 1.0, 100.0),
        )
        expected_prediction, expected_receipt, expected_state = fit_seeded_outer_fold(
            features=inner,
            outcomes=outcome,
            folds=folds,
            block_ids=blocks,
            held_fold=3,
            seed=1103,
            top_k_grid=(4, 16, 32),
            alphas=(0.01, 1.0, 100.0),
        )
        np.testing.assert_array_equal(prediction, expected_prediction)
        self.assertEqual(receipt["selected_feature_count"], expected_receipt["selected_feature_count"])
        self.assertEqual(receipt["selected_alpha"], expected_receipt["selected_alpha"])
        self.assertEqual(receipt["selected_inner_rmse"], expected_receipt["selected_inner_rmse"])
        for key in state:
            np.testing.assert_array_equal(state[key], expected_state[key])

    def test_cross_fold_block_fails_closed(self) -> None:
        inner, outer, outcome, folds, blocks = self.synthetic()
        changed = blocks.copy()
        changed[np.flatnonzero(folds == 3)[0]] = "block-00"
        with self.assertRaisesRegex(DNASeededHeadFitError, "overlaps training"):
            self.fit(inner, outer, outcome, folds, changed)

    def test_sealed_authorization_mutation_fails_closed(self) -> None:
        config = load_config(CONFIG)
        changed = copy.deepcopy(config)
        changed["authorization"]["sealed_asset_access"] = True
        with self.assertRaisesRegex(Exception, "campaign contract"):
            validate_config(changed)


if __name__ == "__main__":
    unittest.main()

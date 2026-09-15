from __future__ import annotations

import copy
from pathlib import Path
import unittest

import numpy as np

from scripts.fit_gse281364_enformer_sei_heads import (
    EnformerSeiHeadFitError,
    fit_seeded_outer_fold,
)
from scripts.freeze_gse281364_enformer_sei_head_taskspec import (
    CANDIDATES,
    PREDICTION_FIELDS,
    load_config,
    validate_config,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/gse281364_enformer_sei_head_campaign.json"


class EnformerSeiHeadCampaignTests(unittest.TestCase):
    @staticmethod
    def synthetic() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        rng = np.random.default_rng(42)
        blocks = np.repeat(np.asarray([f"block-{index:02d}" for index in range(50)]), 2)
        folds = np.repeat(np.arange(50) % 5, 2)
        features = rng.normal(size=(100, 24))
        outcomes = 0.8 * features[:, 3] - 0.4 * features[:, 7] + rng.normal(scale=0.2, size=100)
        return features, outcomes, folds, blocks

    def fit(self, features: np.ndarray, outcomes: np.ndarray, folds: np.ndarray, blocks: np.ndarray, *, seed: int = 1103):
        return fit_seeded_outer_fold(
            features=features,
            outcomes=outcomes,
            folds=folds,
            block_ids=blocks,
            held_fold=2,
            seed=seed,
            top_k_grid=(4, 12, 24),
            alphas=(0.01, 1.0, 100.0),
        )

    def test_frozen_config_and_roster(self) -> None:
        config = load_config(CONFIG)
        validate_config(config)
        self.assertEqual(len(CANDIDATES), 6)
        self.assertNotIn("observed", PREDICTION_FIELDS)
        self.assertFalse(config["completion_firewall"]["shortlist_authorized"])

    def test_held_outcomes_cannot_change_fit_or_prediction(self) -> None:
        features, outcomes, folds, blocks = self.synthetic()
        prediction, receipt, state = self.fit(features, outcomes, folds, blocks)
        poisoned = outcomes.copy()
        poisoned[folds == 2] = 1.0e12
        other_prediction, other_receipt, other_state = self.fit(features, poisoned, folds, blocks)
        np.testing.assert_array_equal(prediction, other_prediction)
        self.assertEqual(receipt, other_receipt)
        for key in state:
            np.testing.assert_array_equal(state[key], other_state[key])

    def test_held_features_do_not_change_preprocessing_or_selection(self) -> None:
        features, outcomes, folds, blocks = self.synthetic()
        _, receipt, state = self.fit(features, outcomes, folds, blocks)
        poisoned = features.copy()
        poisoned[folds == 2] = 1.0e12
        _, other_receipt, other_state = self.fit(poisoned, outcomes, folds, blocks)
        self.assertEqual(receipt, other_receipt)
        for key in state:
            np.testing.assert_array_equal(state[key], other_state[key])

    def test_five_seed_bootstrap_is_not_schema_repetition(self) -> None:
        features, outcomes, folds, blocks = self.synthetic()
        first, _, _ = self.fit(features, outcomes, folds, blocks, seed=1103)
        second, _, _ = self.fit(features, outcomes, folds, blocks, seed=2909)
        self.assertFalse(np.array_equal(first, second))

    def test_cross_fold_block_fails_closed(self) -> None:
        features, outcomes, folds, blocks = self.synthetic()
        bad_blocks = blocks.copy()
        bad_blocks[folds == 2][0:1] = "block-00"
        # Boolean advanced indexing returns a copy, so change the explicit held index.
        bad_blocks[np.flatnonzero(folds == 2)[0]] = "block-00"
        with self.assertRaisesRegex(EnformerSeiHeadFitError, "held-block isolation"):
            self.fit(features, outcomes, folds, bad_blocks)

    def test_authorization_mutation_fails_closed(self) -> None:
        config = load_config(CONFIG)
        changed = copy.deepcopy(config)
        changed["authorization"]["sealed_asset_access"] = True
        with self.assertRaisesRegex(Exception, "authorization boundary"):
            validate_config(changed)


if __name__ == "__main__":
    unittest.main()

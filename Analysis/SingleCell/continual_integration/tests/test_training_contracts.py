from __future__ import annotations

import unittest

import anndata as ad
import numpy as np
import pandas as pd

from masld_cl.training import (
    TrainingError, _reference_fisher_indices, _set_model_labels,
    _validate_held_study_rosters, _validate_scanvi_label_codes,
    _validate_update_trainability,
)


class TestTrainingContracts(unittest.TestCase):
    def test_scanvi_unlabeled_category_is_final_registry_code(self):
        obs = pd.DataFrame({
            "audit_cell_type": ["A", "pDCs", "A", "pDCs"],
            "strict_reference": [True, True, False, False],
        })
        value = ad.AnnData(np.ones((4, 2)), obs=obs)
        _set_model_labels(value, query_unknown=True, unlabeled="Unknown")
        self.assertEqual(
            list(value.obs["model_label"].cat.categories),
            ["A", "pDCs", "Unknown"],
        )
        _validate_scanvi_label_codes(value, unlabeled="Unknown")
        self.assertLess(
            int(value.obs.loc[value.obs["model_label"] != "Unknown", "model_label"].cat.codes.max()),
            2,
        )

    def test_scanvi_nonfinal_unlabeled_category_fails_closed(self):
        obs = pd.DataFrame({
            "model_label": pd.Categorical(
                ["A", "pDCs"], categories=["A", "Unknown", "pDCs"]
            )
        })
        value = ad.AnnData(np.ones((2, 2)), obs=obs)
        with self.assertRaisesRegex(TrainingError, "final registry code"):
            _validate_scanvi_label_codes(value, unlabeled="Unknown")

    def test_scanvi_unfrozen_update_allows_only_fixed_class_prior(self):
        _validate_update_trainability(
            {"encoder.weight": True, "classifier.weight": True, "y_prior": False},
            model_kind="all_lineage", method="continual_learning",
        )
        with self.assertRaises(TrainingError):
            _validate_update_trainability(
                {"encoder.weight": False, "classifier.weight": True, "y_prior": False},
                model_kind="all_lineage", method="continual_learning",
            )
        with self.assertRaises(TrainingError):
            _validate_update_trainability(
                {"encoder.weight": True, "classifier.weight": True, "y_prior": True},
                model_kind="all_lineage", method="continual_learning",
            )

    def test_architecture_surgery_requires_mixed_optimizable_parameters(self):
        _validate_update_trainability(
            {"encoder.weight": False, "classifier.weight": True, "y_prior": False},
            model_kind="all_lineage", method="architecture_surgery",
        )
        with self.assertRaises(TrainingError):
            _validate_update_trainability(
                {"encoder.weight": True, "classifier.weight": True, "y_prior": False},
                model_kind="all_lineage", method="architecture_surgery",
            )

    def test_continual_fisher_is_exact_replay_subset(self):
        obs = pd.DataFrame({
            "donor_id": ["d1"] * 6 + ["d2"] * 6,
            "audit_cell_type": ["A"] * 12,
            "preparation_method": ["p"] * 12,
        })
        value = ad.AnnData(np.ones((12, 2)), obs=obs)
        config = {"sampling": {
            "lineage_fisher_min_cells": 5,
            "primary_replay_fraction": 0.2,
        }}
        replay = np.asarray([0, 2, 6, 8], dtype=np.int64)
        pool = np.arange(12, dtype=np.int64)
        fisher = _reference_fisher_indices(
            value, "continual_learning", replay, pool, config, 17
        )
        np.testing.assert_array_equal(fisher, replay)
        self.assertTrue(set(fisher).issubset(set(replay)))

    def test_low_source_count_group_is_removed_from_fisher(self):
        obs = pd.DataFrame({
            "donor_id": ["d1"] * 5 + ["d2"] * 2,
            "audit_cell_type": ["A"] * 7,
            "preparation_method": ["p"] * 7,
        })
        value = ad.AnnData(np.ones((7, 2)), obs=obs)
        config = {"sampling": {
            "lineage_fisher_min_cells": 5,
            "primary_replay_fraction": 0.2,
        }}
        fisher = _reference_fisher_indices(
            value, "continual_learning", np.asarray([0, 5]),
            np.arange(7), config, 17,
        )
        np.testing.assert_array_equal(fisher, np.asarray([0]))

    def test_held_study_requires_disjoint_query_and_control_rosters(self):
        config = {"evaluation": {
            "powered_query_studies": ["GSE202379", "GSE244832"]
        }}
        held = _validate_held_study_rosters(
            config, ["GSE244832"], ["GSE244832"], ["GSE202379"]
        )
        self.assertEqual(held, {"GSE202379"})
        with self.assertRaisesRegex(TrainingError, "control Fisher"):
            _validate_held_study_rosters(
                config, ["GSE244832"], None, ["GSE202379"]
            )
        with self.assertRaisesRegex(TrainingError, "disjoint"):
            _validate_held_study_rosters(
                config, ["GSE202379"], ["GSE244832"], ["GSE202379"]
            )


if __name__ == "__main__":
    unittest.main()

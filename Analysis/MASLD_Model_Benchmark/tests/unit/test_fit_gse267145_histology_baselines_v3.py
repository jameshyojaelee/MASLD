from __future__ import annotations

import json
import numpy as np
from pathlib import Path
import unittest
from unittest.mock import patch

import scripts.fit_gse267145_histology_baselines_v3 as v3


ROOT = Path(__file__).resolve().parents[2]


class _Cache:
    def get(self, split_id, fitting, evaluation, feature_request):
        return {
            "fit": np.ones((len(fitting), 2), dtype=float),
            "evaluation": np.ones((len(evaluation), 2), dtype=float),
            "max_components": 2,
        }


class _ValidModel:
    coef_ = np.ones((3, 2), dtype=float)

    def predict(self, values):
        return np.zeros(len(values), dtype=int)

    def predict_proba(self, values):
        return np.repeat(np.asarray([[0.5, 0.3, 0.2]]), len(values), axis=0)


class GSE267145HistologyFitterV3Tests(unittest.TestCase):
    def setUp(self) -> None:
        self.assignments = np.repeat(np.arange(4), 3)
        self.groups = np.tile(np.arange(3), 4)
        self.arguments = {
            "cache": _Cache(),
            "representation": {"feature_request": 2, "pca_components": 2},
            "outer_training_indices": np.arange(12),
            "outer_test_indices": np.arange(3),
            "inner_assignment": self.assignments,
            "groups": self.groups,
            "c_values": list(v3.FROZEN_C_VALUES),
            "l1_ratios": list(v3.FROZEN_L1_RATIOS),
            "seed": 10,
        }

    def test_exact_convergence_only_grid_exhaustion_uses_training_priors(self) -> None:
        audits = []
        with patch.object(
            v3.v2,
            "_logistic",
            side_effect=v3.v2.HistologyBaselineError(
                "elastic-net logistic candidate did not converge"
            ),
        ):
            result = v3.tune_group_classifier_v3(
                **self.arguments,
                audit_callback=lambda endpoint, payload: audits.append(
                    (endpoint, payload)
                ),
            )
        np.testing.assert_allclose(result["oof"], np.full((12, 3), 1.0 / 3.0))
        np.testing.assert_allclose(result["test"], np.full((3, 3), 1.0 / 3.0))
        selected = result["selected"]
        self.assertTrue(selected["fallback_triggered"])
        self.assertEqual(selected["fallback_id"], v3.FALLBACK_ID)
        self.assertEqual(selected["candidate_count"], 12)
        self.assertEqual(selected["valid_candidate_count"], 0)
        self.assertEqual(selected["candidate_failure_count"], 12)
        self.assertFalse(selected["outer_test_features_used"])
        self.assertFalse(selected["outer_test_outcomes_used"])
        self.assertFalse(selected["randomness_used"])
        self.assertEqual(audits[0][0], "fibrosis_group3")
        self.assertEqual(len(audits[0][1]["candidates"]), 12)

    def test_one_valid_candidate_preserves_original_grid_selection(self) -> None:
        def fit_or_fail(x, y, *, c_value, l1_ratio, seed):
            if c_value == 0.01 and l1_ratio == 0.0:
                return _ValidModel()
            raise v3.v2.HistologyBaselineError(
                "elastic-net logistic candidate did not converge"
            )

        with patch.object(v3.v2, "_logistic", side_effect=fit_or_fail):
            result = v3.tune_group_classifier_v3(**self.arguments)
        selected = result["selected"]
        self.assertFalse(selected["fallback_triggered"])
        self.assertEqual(selected["parameters"], {"c": 0.01, "l1_ratio": 0.0})
        self.assertEqual(selected["valid_candidate_count"], 1)
        self.assertEqual(selected["candidate_failure_count"], 11)
        np.testing.assert_allclose(result["test"], [[0.5, 0.3, 0.2]] * 3)

    def test_nonfrozen_grid_remains_fail_closed(self) -> None:
        arguments = dict(self.arguments)
        arguments["c_values"] = [0.01]
        arguments["l1_ratios"] = [0.0]
        with patch.object(
            v3.v2,
            "_logistic",
            side_effect=v3.v2.HistologyBaselineError(
                "elastic-net logistic candidate did not converge"
            ),
        ):
            with self.assertRaisesRegex(
                v3.v2.HistologyBaselineError, "exact 12-candidate grid"
            ):
                v3.tune_group_classifier_v3(**arguments)

    def test_nonconvergence_failure_remains_fail_closed(self) -> None:
        with patch.object(
            v3.v2,
            "_logistic",
            side_effect=ValueError("synthetic shape failure"),
        ):
            with self.assertRaisesRegex(
                v3.v2.HistologyBaselineError, "convergence-only failures"
            ):
                v3.tune_group_classifier_v3(**self.arguments)

    def test_contract_preserves_v2_axes_and_unscored_firewall(self) -> None:
        contract = json.loads(v3.CONTRACT_PATH.read_text(encoding="utf-8"))
        preserved = contract["preserved_contract"]
        self.assertEqual(preserved["outer_folds"], [0, 1, 2, 3, 4])
        self.assertEqual(preserved["model_seeds"], [1701, 1709, 1721, 1723, 1733])
        self.assertEqual(preserved["primary_endpoint"], "stage3")
        self.assertFalse(preserved["primary_endpoint_code_path_changed"])
        self.assertFalse(preserved["other_secondary_endpoint_code_paths_changed"])
        self.assertFalse(preserved["outer_test_outcomes_available_to_fit_worker"])
        self.assertFalse(contract["scoring_authorized"])

    def test_wrapper_patches_only_group_selection_and_unit_receipt_paths(self) -> None:
        source = Path(v3.__file__).read_text(encoding="utf-8")
        self.assertNotIn("v2.tune_stage_classifier =", source)
        self.assertNotIn("v2.tune_regression =", source)
        self.assertNotIn("v2.tune_cumulative_fibrosis =", source)
        self.assertNotIn("v2._fit_outer_seed =", source)
        self.assertIn("v2.tune_group_classifier = tune_group_classifier_v3", source)
        self.assertIn("v2._selection_receipt = _selection_receipt_v3", source)
        self.assertIn("v2.fit_preflight = fit_preflight_v3", source)


if __name__ == "__main__":
    unittest.main()

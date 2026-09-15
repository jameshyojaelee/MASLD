from __future__ import annotations

import json
import numpy as np
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import scripts.fit_gse267145_histology_baselines as fitter


class GSE267145HistologyFitterTests(unittest.TestCase):
    def setUp(self) -> None:
        fitter.FIT_COUNTS.clear()
        fitter.ELASTIC_LOGISTIC_ITERATIONS.clear()

    def test_one_standard_error_prefers_lower_complexity(self) -> None:
        selected = fitter.one_standard_error_select(
            [
                {"candidate_id": "complex", "fold_scores": [0.8, 0.7, 0.8, 0.7], "complexity": [10]},
                {"candidate_id": "simple", "fold_scores": [0.74, 0.74, 0.74, 0.74], "complexity": [1]},
            ],
            maximize=True,
        )
        self.assertEqual(selected["candidate_id"], "simple")

    def test_one_invalid_candidate_is_excluded_and_audited(self) -> None:
        selected = fitter.select_with_failure_audit(
            [
                {
                    "candidate_id": "failed",
                    "fold_scores": [],
                    "complexity": [0],
                    "failure_reason": "HistologyBaselineError:did not converge",
                },
                {
                    "candidate_id": "valid",
                    "fold_scores": [0.5, 0.6, 0.5, 0.6],
                    "complexity": [1],
                    "failure_reason": None,
                },
            ],
            maximize=True,
        )
        self.assertEqual(selected["candidate_id"], "valid")
        self.assertEqual(selected["candidate_failure_count"], 1)
        self.assertEqual(selected["valid_candidate_count"], 1)

    def test_all_invalid_candidates_fail_closed(self) -> None:
        with self.assertRaises(fitter.HistologyBaselineError):
            fitter.select_with_failure_audit(
                [
                    {
                        "candidate_id": "failed",
                        "fold_scores": [],
                        "complexity": [0],
                        "failure_reason": "HistologyBaselineError:did not converge",
                    }
                ],
                maximize=False,
            )

    def test_tuner_all_invalid_reports_contract_error_not_complexity_error(self) -> None:
        class Cache:
            def get(self, split_id, fitting, evaluation, feature_request):
                return {
                    "fit": np.ones((len(fitting), 2)),
                    "evaluation": np.ones((len(evaluation), 2)),
                    "max_components": 2,
                }

        assignments = np.repeat(np.arange(4), 3)
        groups = np.tile(np.arange(3), 4)
        with patch.object(
            fitter,
            "_logistic",
            side_effect=fitter.HistologyBaselineError("candidate did not converge"),
        ):
            with self.assertRaisesRegex(
                fitter.HistologyBaselineError, "every hyperparameter candidate failed"
            ):
                fitter.tune_group_classifier(
                    cache=Cache(),
                    representation={"feature_request": 2, "pca_components": 2},
                    outer_training_indices=np.arange(12),
                    outer_test_indices=np.arange(3),
                    inner_assignment=assignments,
                    groups=groups,
                    c_values=[0.1],
                    l1_ratios=[0.0],
                    seed=10,
                )

    def test_cache_key_binds_actual_fit_and_evaluation_indices(self) -> None:
        rng = np.random.default_rng(17)
        matrix = rng.poisson(3, size=(16, 12)).astype(float)
        cache = fitter.RepresentationCache(matrix, [f"g{i}" for i in range(12)], [8], [3], 29)
        first = cache.get("inner0", np.arange(0, 10), np.arange(10, 16), 8)
        repeated = cache.get("inner0", np.arange(0, 10), np.arange(10, 16), 8)
        changed = cache.get("inner0", np.arange(1, 11), np.r_[0, np.arange(11, 16)], 8)
        self.assertIs(first, repeated)
        self.assertIsNot(first, changed)
        self.assertNotEqual(first["fitting_index_sha256"], changed["fitting_index_sha256"])
        self.assertEqual(fitter.FIT_COUNTS["pca_fits"], 2)

    def test_distance_probabilities_are_normalized(self) -> None:
        x = np.asarray([[0, 0], [0, 1], [4, 4], [4, 5], [8, 0], [8, 1]], dtype=float)
        labels = np.asarray([0, 0, 1, 1, 2, 2])
        for kind in ("nearest_centroid", "knn"):
            probabilities = fitter._distance_probabilities(
                kind, x, labels, x, neighbor_count=3, temperature=1.0
            )
            self.assertEqual(probabilities.shape, (6, 3))
            np.testing.assert_allclose(probabilities.sum(axis=1), 1.0)

    def test_elastic_logistic_uses_corrected_iteration_ceiling(self) -> None:
        from sklearn.linear_model import LogisticRegression

        x = np.asarray(
            [[-2.0, 0.0], [-1.0, 0.1], [-0.5, -0.1], [0.5, 0.1], [1.0, -0.1], [2.0, 0.0]]
        )
        labels = np.asarray([0, 0, 0, 1, 1, 1])
        model = fitter._logistic(x, labels, c_value=1.0, l1_ratio=0.5, seed=17)
        previous = LogisticRegression(
            C=1.0,
            penalty="elasticnet",
            solver="saga",
            l1_ratio=0.5,
            max_iter=20_000,
            tol=1e-4,
            fit_intercept=True,
            class_weight=fitter._class_weights(labels),
            random_state=17,
        ).fit(x, labels)
        self.assertEqual(fitter.ELASTIC_LOGISTIC_MAX_ITER, 500_000)
        self.assertEqual(model.max_iter, fitter.ELASTIC_LOGISTIC_MAX_ITER)
        np.testing.assert_allclose(model.coef_, previous.coef_, rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose(
            model.predict_proba(x), previous.predict_proba(x), rtol=1e-10, atol=1e-12
        )
        diagnostics = fitter.elastic_logistic_solver_diagnostics()
        self.assertEqual(diagnostics["previous_max_iter"], 20_000)
        self.assertEqual(diagnostics["intermediate_diagnostic_max_iter"], 100_000)
        self.assertEqual(diagnostics["max_iter"], 500_000)
        self.assertEqual(diagnostics["fitted_candidates"], 1)

    def test_fibrosis_pav_is_monotone(self) -> None:
        projected = fitter._pav_decreasing(np.asarray([[0.2, 0.8, 0.4], [0.9, 0.7, 0.1]]))
        self.assertTrue(np.all(projected[:, :-1] >= projected[:, 1:]))

    def test_model_roster_contains_distance_controls(self) -> None:
        self.assertEqual(len(fitter.MODEL_IDS), 11)
        self.assertIn("rna_hvg_pca_nearest_centroid", fitter.MODEL_IDS)
        self.assertIn("h3_variance_pca_knn", fitter.MODEL_IDS)

    def test_selected_cumulative_predictions_reuse_exact_candidate_seed_path(self) -> None:
        class Cache:
            def get(self, split_id, fitting, evaluation, feature_request):
                return {
                    "fit": np.column_stack((np.arange(len(fitting)), np.ones(len(fitting)))),
                    "evaluation": np.column_stack((np.arange(len(evaluation)), np.ones(len(evaluation)))),
                    "max_components": 2,
                }

        class Model:
            coef_ = np.ones((1, 2))

            def predict_proba(self, values):
                probability = np.full(len(values), 0.4)
                return np.column_stack((1.0 - probability, probability))

        observed_seeds: list[int] = []

        def fake_logistic(x, y, *, c_value, l1_ratio, seed):
            observed_seeds.append(seed)
            return Model()

        audits = []
        with patch.object(fitter, "_logistic", side_effect=fake_logistic):
            result = fitter.tune_cumulative_fibrosis(
                cache=Cache(),
                representation={"feature_request": 2, "pca_components": 2},
                outer_training_indices=np.arange(16),
                outer_test_indices=np.arange(4),
                inner_assignment=np.repeat(np.arange(4), 4),
                fibrosis=np.tile(np.arange(4), 4).astype(float),
                c_values=[0.1],
                l1_ratios=[0.0],
                seed=100,
                audit_callback=lambda endpoint, payload: audits.append((endpoint, payload)),
            )
        expected_inner = [100 + fold * 10 + threshold for fold in range(4) for threshold in (1, 2, 3)]
        expected_outer = [1100 + threshold for threshold in (1, 2, 3)]
        self.assertCountEqual(observed_seeds, expected_inner + expected_outer)
        self.assertEqual(len(observed_seeds), 15)
        self.assertTrue(result["selected"]["seed_lineage"]["selected_predictions_reuse_exact_candidate_fits"])
        self.assertEqual(audits[0][0], "fibrosis_cumulative")
        self.assertTrue(audits[0][1]["candidates"][0]["outer_training_refit_valid"])

    def test_failure_receipt_persists_without_predictions(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            audit_root = output / "outer_0/seed_1701"
            audit_root.mkdir(parents=True)
            (audit_root / "candidate_audit--model--stage3.json").write_text("{}\n")
            path = fitter.write_fit_failure_receipt(
                output, 0, 1701, fitter.HistologyBaselineError("synthetic failure")
            )
            receipt = json.loads(path.read_text())
            self.assertEqual(receipt["candidate_audit_files_persisted"], 1)
            self.assertEqual(receipt["prediction_files_written"], 0)
            self.assertFalse((output / "predictions").exists())


if __name__ == "__main__":
    unittest.main()

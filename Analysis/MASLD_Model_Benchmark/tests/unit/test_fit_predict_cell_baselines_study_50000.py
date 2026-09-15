#!/usr/bin/env python3
"""Unit checks for study-held-out task-native cell-state predictions."""

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np
from scipy import sparse

from scripts import fit_predict_cell_baselines_study_50000 as baseline


ROOT = Path(__file__).parents[2]
WRAPPER = ROOT / "slurm/run_cell_baselines_study_50000.sbatch"


class CellBaselineStudy50000Tests(unittest.TestCase):
    def test_float_storage_accepts_only_integer_valued_raw_counts(self) -> None:
        baseline.validate_integer_like_counts(
            sparse.csr_matrix(np.asarray([[0.0, 1.0], [2.0, 3.0]]))
        )
        with self.assertRaises(baseline.BaselineScreenError):
            baseline.validate_integer_like_counts(
                sparse.csr_matrix(np.asarray([[0.0, 1.5]]))
            )

    def test_centroid_and_knn_probabilities_are_normalized(self) -> None:
        train_x = np.asarray(
            [[class_id * 5.0 + offset, class_id] for class_id in range(5) for offset in (0.0, 0.2)],
            dtype=np.float64,
        )
        train_y = np.repeat(np.arange(5), 2)
        weights = np.ones(len(train_y), dtype=np.float64)
        centroid, state = baseline.nearest_centroid_probabilities(
            train_x, train_y, weights, train_x, classes=5
        )
        knn = baseline.weighted_knn_probabilities(
            train_x, train_y, weights, train_x, classes=5, neighbors=1, n_jobs=1
        )
        self.assertEqual(centroid.shape, (10, 5))
        self.assertGreater(float(state["temperature"]), 0.0)
        np.testing.assert_allclose(centroid.sum(axis=1), 1.0)
        np.testing.assert_allclose(knn.sum(axis=1), 1.0)
        np.testing.assert_array_equal(np.argmax(knn, axis=1), train_y)

    def test_classifier_projection_standardization_is_training_only(self) -> None:
        train = np.asarray([[1.0, 4.0, 9.0], [3.0, 4.0, 13.0]], dtype=np.float32)
        test = np.asarray([[101.0, 400.0, 900.0]], dtype=np.float32)
        train_z, test_z, mean, scale = baseline.standardize_projection(train, test)
        np.testing.assert_allclose(mean, [2.0, 4.0, 11.0])
        np.testing.assert_allclose(scale, [1.0, 1.0, 2.0])
        np.testing.assert_allclose(train_z.mean(axis=0), 0.0, atol=1.0e-7)
        np.testing.assert_allclose(test_z, [[99.0, 396.0, 444.5]])
        self.assertEqual(train_z.dtype, np.float64)
        self.assertEqual(test_z.dtype, np.float64)

    def test_frozen_parameters_and_wrapper_contract(self) -> None:
        self.assertEqual(baseline.N_TOP_HVG, 2_000)
        self.assertEqual(baseline.N_PCA, 50)
        self.assertEqual(baseline.KNN_NEIGHBORS, 15)
        self.assertEqual(baseline.SEEDS, (20260824, 20260825, 20260826))
        self.assertEqual(baseline.LOGISTIC_MAX_ITER, 2_000)
        self.assertEqual(baseline.ELASTIC_NET_MAX_ITER, 5_000)
        self.assertEqual(baseline.ELASTIC_NET_TOL, 1.0e-4)
        self.assertEqual(baseline.ELASTIC_NET_N_ITER_NO_CHANGE, 20)
        self.assertEqual(baseline.ELASTIC_NET_ALPHAS, (1.0e-5, 1.0e-4, 1.0e-3))
        self.assertEqual(baseline.ELASTIC_NET_L1_RATIOS, (0.15, 0.5, 0.85))
        self.assertEqual(len(baseline.MODEL_IDS), 5)
        text = WRAPPER.read_text(encoding="utf-8")
        header = "\n".join(line for line in text.splitlines() if line.startswith("#SBATCH"))
        self.assertIn("--job-name=model-cpu-train-501", header)
        self.assertIn("--partition=cpu", header)
        self.assertIn("--cpus-per-task=16", header)
        self.assertIn("--mem=128G", header)
        self.assertIn("--qos=nslab", header)
        self.assertNotIn("innovation", text.lower())
        self.assertNotIn("--array", header)
        self.assertFalse(
            any(line.lstrip().startswith("sbatch ") for line in text.splitlines())
        )
        self.assertIn("672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509", text)
        self.assertIn("10927e162581375d866adbcf2b3bb7bd0dc4fda9fd6dfe033198df7d95be2680", text)
        self.assertIn("8a7dfd2a15a635d5ef5b25fa8ee68355e6c6cfd047ec464cdc1dea83b7201298", text)
        self.assertIn("freeze_python_environment.py", text)
        self.assertNotIn(".agent/scanpy_env", text)

    def test_study_firewall_and_inner_preprocessing_are_explicit(self) -> None:
        text = (ROOT / "scripts/fit_predict_cell_baselines_study_50000.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("set(datasets[train]) & set(datasets[test])", text)
        self.assertIn("raw_train[fitting]", text)
        self.assertIn("fitting_weights", text)
        self.assertIn("oof_seen[validation] += 1", text)
        self.assertIn("inner_svm_preprocessing_fit_on_inner_training_only", text)
        self.assertIn(
            "inner_svm_projection_standardization_fit_on_inner_training_only", text
        )
        self.assertIn("prepare_inner_projections", text)
        self.assertIn("select_sgd_elastic_net", text)
        self.assertIn("donor_class_balanced_multiclass_log_loss", text)
        self.assertNotIn('solver="saga"', text)

    def test_fold_models_complete_inner_cross_fit(self) -> None:
        rng = np.random.default_rng(17)
        donor_roster = np.asarray([f"donor-{index:02d}" for index in range(25)])
        train_donors = np.repeat(donor_roster, 5)
        train_y = np.tile(np.arange(5, dtype=np.int64), 25)
        raw_train = sparse.csr_matrix(
            rng.poisson(2.0, size=(len(train_y), baseline.N_TOP_HVG)).astype(np.float32)
        )
        train_x = rng.normal(size=(len(train_y), 50)).astype(np.float32)
        test_x = rng.normal(size=(10, 50)).astype(np.float32)
        assignment = {
            donor: index % 5 for index, donor in enumerate(donor_roster.tolist())
        }
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            values = baseline._fit_fold_models(
                raw_train=raw_train,
                gene_ids=[
                    f"ENSG{index:011d}" for index in range(baseline.N_TOP_HVG)
                ],
                train_x=train_x,
                train_y=train_y,
                train_donors=train_donors,
                train_row_ids=np.asarray([f"row-{index}" for index in range(len(train_y))]),
                test_x=test_x,
                inner_assignment=assignment,
                fold=0,
                seed=baseline.SEEDS[0],
                output=output,
            )
            self.assertTrue((output / "classifier_projection_mean.npy").is_file())
            self.assertTrue((output / "classifier_projection_scale.npy").is_file())
        self.assertEqual(set(values), set(baseline.MODEL_IDS))
        for probabilities in values.values():
            self.assertEqual(probabilities.shape, (10, 5))
            np.testing.assert_allclose(probabilities.sum(axis=1), 1.0)

    def test_outer_fit_follows_inner_rank_until_training_convergence(self) -> None:
        x = np.zeros((10, 2), dtype=np.float64)
        y = np.tile(np.arange(5), 2)
        weights = np.ones(10, dtype=np.float64)
        audit = [
            {
                "alpha": 1.0e-4,
                "l1_ratio": 0.5,
                "valid_all_inner_folds": True,
                "mean_inner_log_loss": 0.8,
            },
            {
                "alpha": 1.0e-3,
                "l1_ratio": 0.85,
                "valid_all_inner_folds": True,
                "mean_inner_log_loss": 0.9,
            },
        ]
        failed = {"converged_before_ceiling": False, "iterations": 5000}
        passed = {"converged_before_ceiling": True, "iterations": 37}
        models = [mock.Mock(), mock.Mock()]
        with mock.patch.object(
            baseline,
            "fit_sgd_elastic_net",
            side_effect=[(models[0], failed), (models[1], passed)],
        ) as fitted:
            model, state, selected, attempts = (
                baseline.fit_ranked_outer_sgd_elastic_net(
                    x, y, weights, audit, fold=3, seed=20260826
                )
            )
        self.assertIs(model, models[1])
        self.assertIs(state, passed)
        self.assertEqual(selected, {"alpha": 1.0e-3, "l1_ratio": 0.85})
        self.assertEqual([item["inner_validation_rank"] for item in attempts], [1, 2])
        self.assertEqual(fitted.call_count, 2)
        self.assertEqual(
            [call.kwargs["seed"] for call in fitted.call_args_list],
            [20260829, 20260829],
        )


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.evaluate_sequence_accessibility import (
    SequenceEvaluationError,
    bind_prediction_model_id,
    fisher_mean,
    load_prediction_manifest,
    multinomial_deviance_per_insertion,
    smooth_distribution,
    spearman_or_zero,
    validate_evaluation_role,
    validate_evaluator_input_metadata,
    within_window_deviance_per_insertion,
)


class SequenceAccessibilityEvaluatorTests(unittest.TestCase):
    def test_only_validation_role_is_authorized_before_finalist_lock(self) -> None:
        self.assertEqual(validate_evaluation_role("valid"), "valid")
        with self.assertRaises(SequenceEvaluationError):
            validate_evaluation_role("test")

    def test_model_id_can_bind_from_verified_root_when_subview_omits_it(self) -> None:
        self.assertEqual(bind_prediction_model_id({}, "chrombpnet"), "chrombpnet")

    def test_model_id_rejects_conflicting_subview_declaration(self) -> None:
        with self.assertRaises(SequenceEvaluationError):
            bind_prediction_model_id({"model_id": "other"}, "chrombpnet")

    def test_distribution_smoothing_is_positive_and_normalized(self) -> None:
        values = smooth_distribution([0.0, 2.0, 0.0])
        self.assertTrue(np.all(values > 0))
        self.assertAlmostEqual(float(values.sum()), 1.0)

    def test_perfect_multinomial_prediction_has_zero_deviance(self) -> None:
        truth = np.asarray([1.0, 2.0, 3.0])
        value = multinomial_deviance_per_insertion(truth, truth)
        self.assertLess(value, 1.0e-10)

    def test_within_window_deviance_ignores_between_window_scale(self) -> None:
        truth = np.asarray([[1.0, 3.0], [4.0, 2.0]])
        predicted = truth * np.asarray([[10.0], [0.2]])
        self.assertLess(within_window_deviance_per_insertion(truth, predicted), 1.0e-10)

    def test_zero_mass_prediction_uses_finite_uniform_fallback(self) -> None:
        truth = np.asarray([[1.0, 3.0], [4.0, 2.0]])
        predicted = np.zeros_like(truth)
        value = within_window_deviance_per_insertion(truth, predicted)
        self.assertTrue(np.isfinite(value))
        self.assertGreater(value, 0.0)

    def test_spearman_constant_baseline_is_zero(self) -> None:
        self.assertEqual(spearman_or_zero([1, 2, 3], [1, 1, 1]), 0.0)

    def test_fisher_mean_is_symmetric(self) -> None:
        self.assertAlmostEqual(fisher_mean([-0.5, 0.5]), 0.0)

    def test_view_manifest_separates_candidate_from_root_model(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "predictions.tsv"
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(
                    (
                        "candidate_id",
                        "root_model_id",
                        "root",
                        "artifacts_sha256",
                        "prediction_subdir",
                        "seed",
                    )
                )
                writer.writerow(
                    (
                        "chrombpnet_full__s17",
                        "chrombpnet",
                        "/tmp/model",
                        "a" * 64,
                        "predictions/full_model",
                        "17",
                    )
                )
                writer.writerow(
                    (
                        "chrombpnet_nobias__s17",
                        "chrombpnet",
                        "/tmp/model",
                        "a" * 64,
                        "predictions/nobias_model",
                        "17",
                    )
                )
            rows = load_prediction_manifest(path)
        self.assertEqual(len(rows), 2)
        self.assertEqual({row["root_model_id"] for row in rows}, {"chrombpnet"})

    def test_evaluator_input_metadata_rejects_wrong_split(self) -> None:
        from scripts.evaluate_sequence_accessibility import sha256_file

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            payload = {
                "metadata": {
                    "artifact_class": "sequence_training_mean_baseline",
                    "dataset_id": "gse296875",
                    "split_id": "donor1_genomic1",
                    "lineage_id": "hepatocyte",
                    "biological_unit": "donor",
                    "status": "passed",
                    "benchmark_metrics_calculated": False,
                    "model_id": "training_pseudobulk_mean",
                    "held_donor_outcomes_exposed": False,
                }
            }
            (root / "ARTIFACTS.json").write_text(json.dumps(payload))
            with self.assertRaises(SequenceEvaluationError):
                validate_evaluator_input_metadata(
                    root=root,
                    expected_artifacts_sha256=sha256_file(root / "ARTIFACTS.json"),
                    artifact_class="sequence_training_mean_baseline",
                    dataset_id="gse296875",
                    split_id="donor0_genomic0",
                    lineage_id="hepatocyte",
                )


if __name__ == "__main__":
    unittest.main()

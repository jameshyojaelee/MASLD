from __future__ import annotations

import unittest

from scripts.audit_gse267145_histology_rna_finalist_selection import (
    RNAFinalistSelectionError,
    rederive_one_standard_error,
    select_rna_finalist,
)


MODELS = (
    "rna_hvg_pca_elastic_net",
    "rna_hvg_pca_linear_svm",
    "rna_hvg_pca_nearest_centroid",
    "rna_hvg_pca_knn",
)


class RNAFinalistSelectionTests(unittest.TestCase):
    def test_primary_and_calibration_select_linear_svm(self) -> None:
        rows = [
            {"model_id": MODELS[0], "stage3_macro_f1": 0.56, "stage3_multiclass_brier": 0.49},
            {"model_id": MODELS[1], "stage3_macro_f1": 0.68, "stage3_multiclass_brier": 0.42},
            {"model_id": MODELS[2], "stage3_macro_f1": 0.65, "stage3_multiclass_brier": 0.66},
            {"model_id": MODELS[3], "stage3_macro_f1": 0.58, "stage3_multiclass_brier": 0.56},
        ]
        result = select_rna_finalist(rows, MODELS)
        self.assertEqual(result["selected_model_id"], MODELS[1])
        self.assertTrue(result["unique_primary_winner"])
        self.assertTrue(result["unique_calibration_winner"])

    def test_calibration_cannot_rescue_nonprimary_winner(self) -> None:
        rows = [
            {"model_id": model, "stage3_macro_f1": value, "stage3_multiclass_brier": brier}
            for model, value, brier in zip(
                MODELS,
                (0.56, 0.68, 0.65, 0.58),
                (0.49, 0.50, 0.40, 0.56),
                strict=True,
            )
        ]
        with self.assertRaises(RNAFinalistSelectionError):
            select_rna_finalist(rows, MODELS)

    def test_one_standard_error_rederived_from_four_inner_folds(self) -> None:
        candidates = [
            {
                "candidate_id": "best",
                "parameters": {"pca": 20},
                "fold_scores": [0.70, 0.60, 0.80, 0.70],
                "complexity": [20, 20],
                "failure_reason": None,
                "outer_training_refit_valid": True,
            },
            {
                "candidate_id": "simple",
                "parameters": {"pca": 5},
                "fold_scores": [0.68, 0.66, 0.67, 0.67],
                "complexity": [5, 5],
                "failure_reason": None,
                "outer_training_refit_valid": True,
            },
        ]
        audit = {
            "candidates": candidates,
            "selected": {
                "candidate_id": "simple",
                "best_mean_score": 0.7,
                "best_standard_error": 0.04082482904638631,
                "one_standard_error_threshold": 0.6591751709536137,
                "one_standard_error_candidate_count": 2,
            },
        }
        result = rederive_one_standard_error(audit)
        self.assertEqual(result["selected_candidate_id"], "simple")
        self.assertEqual(result["one_standard_error_candidate_count"], 2)


if __name__ == "__main__":
    unittest.main()

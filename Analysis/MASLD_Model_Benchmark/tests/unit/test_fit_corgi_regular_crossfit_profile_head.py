from __future__ import annotations

import unittest

import numpy as np

from scripts.fit_corgi_regular_crossfit_profile_head import (
    CorgiHeadFitError,
    geometric_profile,
    select_alpha,
    validate_contract_receipt,
)


class CorgiCrossfitProfileHeadTests(unittest.TestCase):
    def test_alpha_zero_is_training_mean_and_alpha_one_is_corgi(self) -> None:
        corgi = np.array([9.0, 3.0, 1.0])
        baseline = np.array([1.0, 3.0, 9.0])
        np.testing.assert_allclose(
            geometric_profile(corgi, baseline, 0.0),
            np.array([1.0, 3.0, 9.0]) / 13.0,
            atol=1.0e-6,
        )
        np.testing.assert_allclose(
            geometric_profile(corgi, baseline, 1.0),
            np.array([9.0, 3.0, 1.0]) / 13.0,
            atol=1.0e-6,
        )

    def test_one_standard_error_prefers_smaller_alpha(self) -> None:
        losses = {
            0.0: [1.03, 1.01, 1.02, 1.04],
            0.25: [1.00, 1.01, 1.01, 1.02],
            0.5: [0.99, 1.01, 1.00, 1.02],
            0.75: [1.00, 1.02, 1.01, 1.03],
            1.0: [1.04, 1.02, 1.03, 1.05],
        }
        selected, rows = select_alpha(losses)
        self.assertEqual(selected, 0.25)
        self.assertEqual(sum(int(row["selected"]) for row in rows), 1)

    def test_invalid_alpha_is_rejected(self) -> None:
        with self.assertRaises(CorgiHeadFitError):
            geometric_profile([1, 2, 3], [3, 2, 1], 0.33)

    def test_contract_rejects_fit_fold_containing_evaluation_fold(self) -> None:
        value = self._receipt()
        value["folds"][0]["fit_valid_folds"] = [0, 2, 3, 4]
        with self.assertRaises(CorgiHeadFitError):
            validate_contract_receipt(value)

    def test_contract_accepts_five_crossed_folds(self) -> None:
        self.assertEqual(len(validate_contract_receipt(self._receipt())), 5)

    @staticmethod
    def _receipt() -> dict[str, object]:
        folds = []
        for evaluation in range(5):
            fit_valid = [fold for fold in range(5) if fold != evaluation]
            folds.append(
                {
                    "evaluation_fold": evaluation,
                    "evaluation_base_outer_fold": (evaluation - 1) % 5,
                    "fit_valid_folds": fit_valid,
                    "fit_base_outer_folds": [(fold - 1) % 5 for fold in fit_valid],
                }
            )
        return {
            "schema_version": "masld-bench-corgi-crossfit-profile-head-preflight-v1",
            "status": "pass_outcome_free_crossfit_head_contract",
            "model_id": "corgi_regular",
            "task_id": "rna_conditioned_atac",
            "dataset_id": "gse296875",
            "lineage_id": "hepatocyte",
            "adaptation_rung": "new_assay_head",
            "donor_and_genomic_fold_crossfit": True,
            "source_signal_arrays_read": False,
            "development_outcomes_read": False,
            "test_or_sealed_features_or_outcomes_read": False,
            "model_fit_performed": False,
            "evaluation_performed": False,
            "promotion_gate_evaluated": False,
            "champion_or_external_claim_allowed": False,
            "global_census_modified": False,
            "folds": folds,
        }


if __name__ == "__main__":
    unittest.main()

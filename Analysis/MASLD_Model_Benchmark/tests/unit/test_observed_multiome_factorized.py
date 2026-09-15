from __future__ import annotations

import inspect
from pathlib import Path
import tempfile
import unittest

import numpy as np

from masld_bench.adapters.observed_multiome_factorized import (
    MODEL_IDS,
    FactorizedBaselineError,
    export,
    fit,
    load,
    predict,
)


class ObservedMultiomeFactorizedTests(unittest.TestCase):
    def setUp(self) -> None:
        generator = np.random.default_rng(19)
        self.rows = [f"row-{index}" for index in range(12)]
        self.strata = [f"lineage-{index % 2}" for index in range(12)]
        self.targets = [f"train-target-{index}" for index in range(7)]
        self.held = [f"held-target-{index}" for index in range(3)]
        self.rna = generator.normal(size=(12, 3))
        self.atac = generator.poisson(3.0, size=(12, 4)).astype(float)
        self.sequence = generator.normal(size=(7, 2))
        self.held_sequence = generator.normal(size=(3, 2))
        self.covariates = generator.normal(size=(7, 2))
        self.held_covariates = generator.normal(size=(3, 2))
        self.counts = generator.poisson(2.0, size=(12, 7)).astype(float)

    def _state(self, model_id: str):
        return fit(
            model_id,
            row_ids=self.rows,
            strata=self.strata,
            target_ids=self.targets,
            rna=self.rna,
            observed_atac=self.atac,
            sequence_features=self.sequence,
            target_covariates=self.covariates,
            target_counts=self.counts,
            seed=23,
        )

    def test_all_models_predict_disjoint_held_targets_without_target_coefficients(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            for model_id in MODEL_IDS:
                state = self._state(model_id)
                output = predict(
                    state,
                    row_ids=self.rows,
                    strata=self.strata,
                    target_ids=self.held,
                    rna=self.rna,
                    observed_atac=self.atac,
                    sequence_features=self.held_sequence,
                    target_covariates=self.held_covariates,
                )
                self.assertEqual(output.shape, (12, 3))
                self.assertEqual(state.coefficients.ndim, 1)
                export(state, Path(directory) / model_id)
                with np.load(Path(directory) / model_id / "state.npz", allow_pickle=False) as arrays:
                    self.assertEqual(arrays["coefficients"].ndim, 1)

    def test_export_load_round_trip_preserves_predictions_and_target_firewall(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            for model_id in MODEL_IDS:
                state = self._state(model_id)
                path = Path(directory) / model_id
                export(state, path)
                restored = load(path)
                expected = predict(
                    state,
                    row_ids=self.rows,
                    strata=self.strata,
                    target_ids=self.held,
                    rna=self.rna,
                    observed_atac=self.atac,
                    sequence_features=self.held_sequence,
                    target_covariates=self.held_covariates,
                )
                observed = predict(
                    restored,
                    row_ids=self.rows,
                    strata=self.strata,
                    target_ids=self.held,
                    rna=self.rna,
                    observed_atac=self.atac,
                    sequence_features=self.held_sequence,
                    target_covariates=self.held_covariates,
                )
                np.testing.assert_array_equal(expected, observed)
                with self.assertRaises(FactorizedBaselineError):
                    predict(
                        restored,
                        row_ids=self.rows,
                        strata=self.strata,
                        target_ids=[self.targets[0]],
                        rna=self.rna,
                        observed_atac=self.atac,
                        sequence_features=self.sequence[:1],
                        target_covariates=self.covariates[:1],
                    )

    def test_prediction_rejects_training_target_overlap(self) -> None:
        state = self._state("rna_only")
        with self.assertRaises(FactorizedBaselineError):
            predict(
                state,
                row_ids=self.rows,
                strata=self.strata,
                target_ids=[self.targets[0], self.held[0]],
                rna=self.rna,
                observed_atac=self.atac,
                sequence_features=np.vstack((self.sequence[0], self.held_sequence[0])),
                target_covariates=np.vstack((self.covariates[0], self.held_covariates[0])),
            )

    def test_prediction_has_no_outcome_argument(self) -> None:
        arguments = inspect.signature(predict).parameters
        self.assertNotIn("target", arguments)
        self.assertNotIn("target_counts", arguments)

    def test_shuffled_prediction_is_row_order_invariant(self) -> None:
        state = self._state("shuffled_modality")
        original = predict(
            state,
            row_ids=self.rows,
            strata=self.strata,
            target_ids=self.held,
            rna=self.rna,
            observed_atac=self.atac,
            sequence_features=self.held_sequence,
            target_covariates=self.held_covariates,
        )
        order = np.asarray([5, 0, 7, 1, 9, 3, 10, 2, 11, 4, 8, 6])
        reordered = predict(
            state,
            row_ids=[self.rows[index] for index in order],
            strata=[self.strata[index] for index in order],
            target_ids=self.held,
            rna=self.rna[order],
            observed_atac=self.atac[order],
            sequence_features=self.held_sequence,
            target_covariates=self.held_covariates,
        )
        np.testing.assert_array_equal(original, reordered[np.argsort(order)])

    def test_required_atac_is_not_zero_filled(self) -> None:
        state = self._state("observed_atac_only")
        with self.assertRaises(FactorizedBaselineError):
            predict(
                state,
                row_ids=self.rows,
                strata=self.strata,
                target_ids=self.held,
                rna=self.rna,
                observed_atac=None,
                sequence_features=self.held_sequence,
                target_covariates=self.held_covariates,
            )


if __name__ == "__main__":
    unittest.main()

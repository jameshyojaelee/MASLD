from __future__ import annotations

import inspect
from pathlib import Path
import tempfile
import unittest

import numpy as np

from masld_bench.adapters.observed_multiome_baselines import (
    MODEL_IDS,
    ObservedMultiomeBaselineError,
    export,
    fit,
    predict,
)


class ObservedMultiomeBaselineTests(unittest.TestCase):
    def setUp(self) -> None:
        generator = np.random.default_rng(7)
        self.rows = [f"row-{index}" for index in range(12)]
        self.strata = [f"lineage-{index % 2}" for index in range(12)]
        self.rna = generator.normal(size=(12, 4))
        self.sequence = generator.normal(size=(12, 3))
        self.atac = generator.poisson(3.0, size=(12, 5)).astype(float)
        self.target = generator.poisson(2.0, size=(12, 2)).astype(float)

    def _state(self, model_id: str):
        return fit(
            model_id,
            row_ids=self.rows,
            strata=self.strata,
            rna=self.rna,
            sequence=self.sequence,
            observed_atac=self.atac,
            target_counts=self.target,
            seed=11,
        )

    def test_all_five_controls_fit_predict_and_export_without_pickle(self) -> None:
        self.assertEqual(len(MODEL_IDS), 5)
        with tempfile.TemporaryDirectory() as directory:
            for model_id in MODEL_IDS:
                state = self._state(model_id)
                output = predict(
                    state,
                    row_ids=self.rows,
                    strata=self.strata,
                    rna=self.rna,
                    sequence=self.sequence,
                    observed_atac=self.atac,
                )
                self.assertEqual(output.shape, self.target.shape)
                self.assertTrue(np.isfinite(output).all())
                self.assertTrue(np.all(output >= 0))
                export(state, Path(directory) / model_id)
                with np.load(Path(directory) / model_id / "state.npz", allow_pickle=False) as arrays:
                    self.assertIn("coefficients", arrays.files)

    def test_predict_has_no_target_argument(self) -> None:
        self.assertNotIn("target", inspect.signature(predict).parameters)
        self.assertNotIn("target_counts", inspect.signature(predict).parameters)

    def test_removed_modalities_cannot_change_predictions(self) -> None:
        for model_id in ("masked_modality", "rna_only"):
            state = self._state(model_id)
            left = predict(
                state,
                row_ids=self.rows,
                strata=self.strata,
                rna=self.rna,
                sequence=self.sequence,
                observed_atac=self.atac,
            )
            right = predict(
                state,
                row_ids=self.rows,
                strata=self.strata,
                rna=self.rna,
                sequence=self.sequence,
                observed_atac=self.atac + 1000,
            )
            np.testing.assert_array_equal(left, right)

    def test_required_missing_atac_is_rejected_not_zero_filled(self) -> None:
        state = self._state("observed_atac_only")
        with self.assertRaises(ObservedMultiomeBaselineError):
            predict(
                state,
                row_ids=self.rows,
                strata=self.strata,
                rna=self.rna,
                sequence=self.sequence,
                observed_atac=None,
            )

    def test_shuffled_query_is_row_keyed_and_order_invariant(self) -> None:
        state = self._state("shuffled_modality")
        original = predict(
            state,
            row_ids=self.rows,
            strata=self.strata,
            rna=self.rna,
            sequence=self.sequence,
            observed_atac=self.atac,
        )
        order = np.asarray([5, 0, 7, 1, 9, 3, 10, 2, 11, 4, 8, 6])
        reordered = predict(
            state,
            row_ids=[self.rows[index] for index in order],
            strata=[self.strata[index] for index in order],
            rna=self.rna[order],
            sequence=self.sequence[order],
            observed_atac=self.atac[order],
        )
        inverse = np.argsort(order)
        np.testing.assert_array_equal(original, reordered[inverse])


if __name__ == "__main__":
    unittest.main()

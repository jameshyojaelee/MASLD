from __future__ import annotations

import inspect
from pathlib import Path
import tempfile
import unittest

import numpy as np

from masld_bench.gse244832_task_native_baselines import (
    MODEL_IDS,
    TaskNativeBaselineError,
    export_state,
    fit,
    load_state,
    predict,
)


class GSE244832TaskNativeBaselinesTests(unittest.TestCase):
    def setUp(self) -> None:
        rng = np.random.default_rng(20260825)
        self.lineages = np.tile(np.arange(4, dtype=np.int64), 8)
        lineage_signal = self.lineages[:, None] * np.linspace(0.0, 0.8, 48)[None, :]
        self.context = rng.poisson(np.exp(-0.2 + lineage_signal)).astype(np.float64)
        self.target = rng.poisson(np.exp(-0.1 + lineage_signal[:, ::-1])).astype(np.float64)
        self.context[:, 0] += 1
        self.target[:, 0] += 1
        self.query = self.context[:8].copy()
        self.query_lineages = self.lineages[:8].copy()

    def test_predict_has_no_target_argument(self) -> None:
        self.assertNotIn("target_counts", inspect.signature(predict).parameters)

    def test_all_three_models_are_finite_repeatable_and_composition_invariant(self) -> None:
        for model_id in MODEL_IDS:
            state = fit(
                model_id,
                context_counts=self.context,
                target_counts=self.target,
                lineage_indices=self.lineages,
                lsi_components=4,
            )
            first = predict(
                state,
                context_counts=self.query,
                lineage_indices=self.query_lineages,
            )
            repeat = predict(
                state,
                context_counts=self.query,
                lineage_indices=self.query_lineages,
            )
            individual = np.vstack(
                [
                    predict(
                        state,
                        context_counts=self.query[index : index + 1],
                        lineage_indices=self.query_lineages[index : index + 1],
                    )
                    for index in range(self.query.shape[0])
                ]
            )
            self.assertEqual(first.shape, (8, 48))
            self.assertTrue(np.isfinite(first).all())
            self.assertTrue((first >= 0).all())
            np.testing.assert_array_equal(first, repeat)
            np.testing.assert_allclose(first, individual, rtol=1.0e-12, atol=1.0e-12)

    def test_poisson_glm_is_closed_form_lineage_rate(self) -> None:
        state = fit(
            "observed_atac_glm",
            context_counts=self.context,
            target_counts=self.target,
            lineage_indices=self.lineages,
        )
        prediction = predict(
            state,
            context_counts=self.query,
            lineage_indices=self.query_lineages,
        )
        for row, lineage in enumerate(self.query_lineages):
            selected = self.lineages == lineage
            expected_rate = self.target[selected].sum(axis=0) / self.context[selected].sum()
            np.testing.assert_allclose(
                prediction[row], self.query[row].sum() * expected_rate, rtol=1.0e-12
            )

    def test_export_load_round_trip(self) -> None:
        state = fit(
            "lsi",
            context_counts=self.context,
            target_counts=self.target,
            lineage_indices=self.lineages,
            lsi_components=4,
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "state"
            export_state(state, path)
            restored = load_state(path)
            observed = predict(
                restored,
                context_counts=self.query,
                lineage_indices=self.query_lineages,
            )
            expected = predict(
                state,
                context_counts=self.query,
                lineage_indices=self.query_lineages,
            )
            np.testing.assert_array_equal(observed, expected)

    def test_observed_zero_query_is_rejected_not_encoded_as_missing_zero(self) -> None:
        state = fit(
            "observed_atac_glm",
            context_counts=self.context,
            target_counts=self.target,
            lineage_indices=self.lineages,
        )
        with self.assertRaises(TaskNativeBaselineError):
            predict(
                state,
                context_counts=np.zeros((1, self.context.shape[1])),
                lineage_indices=np.asarray([0]),
            )


if __name__ == "__main__":
    unittest.main()

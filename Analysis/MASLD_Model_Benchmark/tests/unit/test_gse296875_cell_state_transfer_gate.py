from __future__ import annotations

import tempfile
from pathlib import Path
import unittest

import numpy as np

from scripts import evaluate_gse296875_cell_state_transfer as evaluation
from scripts import lock_gse296875_cell_state_transfer_predictions as lock
from scripts import predict_gse296875_cell_state_baseline_outer_ensemble as baseline
from scripts import predict_gse296875_transcriptformer_outer_ensemble as transcriptformer


class GSE296875CellStateTransferGateTest(unittest.TestCase):
    def test_probability_normalization_rejects_nonfinite(self) -> None:
        values = np.full((3, 5), 0.2)
        self.assertTrue(np.allclose(baseline._probabilities(values, 3), values))
        values[0, 0] = np.nan
        with self.assertRaisesRegex(baseline.BaselineTransferError, "probability tensor"):
            baseline._probabilities(values, 3)

    def test_linear_and_two_layer_head_states_are_both_supported(self) -> None:
        rng = np.random.default_rng(13)
        features = rng.normal(size=(3, 2_048)).astype(np.float32)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            linear = root / "linear.npz"
            np.savez(
                linear,
                feature_mean=np.zeros(2_048, dtype=np.float32),
                feature_scale=np.ones(2_048, dtype=np.float32),
                **{
                    "state::weight": np.zeros((5, 2_048), dtype=np.float32),
                    "state::bias": np.arange(5, dtype=np.float32),
                },
            )
            linear_logits = transcriptformer._head_logits(linear, features, "linear")
            self.assertEqual(linear_logits.shape, (3, 5))
            self.assertTrue(np.allclose(linear_logits, np.arange(5)[None, :]))
            mlp = root / "mlp.npz"
            np.savez(
                mlp,
                feature_mean=np.zeros(2_048, dtype=np.float32),
                feature_scale=np.ones(2_048, dtype=np.float32),
                **{
                    "state::0.weight": np.zeros((256, 2_048), dtype=np.float32),
                    "state::0.bias": np.zeros(256, dtype=np.float32),
                    "state::3.weight": np.zeros((5, 256), dtype=np.float32),
                    "state::3.bias": np.arange(5, dtype=np.float32),
                },
            )
            mlp_logits = transcriptformer._head_logits(mlp, features, "two_layer_mlp")
            self.assertEqual(mlp_logits.shape, (3, 5))
            self.assertTrue(np.allclose(mlp_logits, np.arange(5)[None, :]))

    def test_prediction_lock_rejects_evaluator_columns(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "predictions.tsv"
            path.write_text(
                "row_id\tdonor_id\tpredicted_class\t"
                + "\t".join(f"probability::{label}" for label in lock.ROSTER)
                + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(lock.PredictionLockError, "schema differs"):
                lock._validate_prediction(path, [])

    def test_transfer_roster_keeps_both_common_heads(self) -> None:
        self.assertEqual(
            evaluation.TF_MODELS,
            (
                "transcriptformer_tf_sapiens__linear",
                "transcriptformer_tf_sapiens__two_layer_mlp",
            ),
        )

    def test_transcriptformer_binds_study_split_not_activation_folds(self) -> None:
        text = Path(transcriptformer.__file__).read_text(encoding="utf-8")
        self.assertIn("source_activation_outer", text)
        self.assertIn("activation_split_mismatch_rows != 42_584", text)
        self.assertIn("resource_atlas_study_outer_5fold_v1", text)

    def test_prediction_scripts_do_not_name_evaluator_label_artifact(self) -> None:
        root = Path(__file__).parents[2] / "scripts"
        for name in (
            "predict_gse296875_cell_state_baseline_outer_ensemble.py",
            "predict_gse296875_transcriptformer_outer_ensemble.py",
        ):
            text = (root / name).read_text(encoding="utf-8")
            self.assertNotIn("labels.tsv.gz", text)
            self.assertNotIn("/evaluator/", text)


if __name__ == "__main__":
    unittest.main()

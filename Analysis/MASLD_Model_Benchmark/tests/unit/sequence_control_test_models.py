from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts import sequence_control_architectures as architectures
from scripts import sequence_control_predict_ccre as prediction
from scripts import sequence_control_prepare_training as preparation
from scripts import sequence_control_run_training as training


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config" / "sequence_control_architectures.json"


@dataclass
class Arguments:
    seed: int = 20260824
    learning_rate: float = 0.001


def parameters(model_id: str) -> dict[str, str]:
    raw = json.loads(CONFIG.read_text(encoding="utf-8"))
    model = raw["architectures"][model_id]
    values = {
        "control_model_id": model_id,
        "control_blocks": str(model["blocks"]),
        "control_dropout": str(model["dropout"]),
        "control_stem_kernel_size": str(model["stem_kernel_size"]),
        "counts_loss_weight": "10.0",
        "inputlen": "2114",
        "outputlen": "1000",
    }
    if model_id == "sequence_cnn_control":
        values.update(
            {
                "control_channels": str(model["channels"]),
                "control_kernel_size": str(model["kernel_size"]),
            }
        )
    else:
        values.update(
            {
                "control_attention_dropout": str(model["attention_dropout"]),
                "control_ffn_width": str(model["ffn_width"]),
                "control_heads": str(model["heads"]),
                "control_profile_kernel_size": str(model["profile_kernel_size"]),
                "control_token_stride": str(model["token_stride"]),
                "control_width": str(model["width"]),
            }
        )
    return values


class SequenceControlModelTests(unittest.TestCase):
    def test_models_have_matched_geometry_and_parameter_budgets(self) -> None:
        import tensorflow as tf

        raw = json.loads(CONFIG.read_text(encoding="utf-8"))
        budget = raw["parameter_budget"]
        counts = {}
        for model_id in architectures.MODEL_IDS:
            tf.keras.backend.clear_session()
            model = architectures.getModelGivenModelOptionsAndWeightInits(
                Arguments(), parameters(model_id)
            )
            self.assertEqual(model.name, model_id)
            self.assertEqual(model.input_shape, (None, 2114, 4))
            self.assertEqual(
                [tuple(shape) for shape in model.output_shape],
                [(None, 1000), (None, 1)],
            )
            counts[model_id] = model.count_params()
            self.assertGreaterEqual(
                model.count_params(), budget["minimum_trainable_parameters"]
            )
            self.assertLessEqual(
                model.count_params(), budget["maximum_trainable_parameters"]
            )
        ratio = max(counts.values()) / min(counts.values())
        self.assertLessEqual(ratio, budget["maximum_pairwise_ratio"])

    def test_initial_predictions_are_seed_deterministic(self) -> None:
        import tensorflow as tf

        inputs = np.eye(4, dtype=np.float32)[
            np.random.default_rng(7).integers(0, 4, size=(1, 2114))
        ]
        for model_id in architectures.MODEL_IDS:
            outputs = []
            for _ in range(2):
                tf.keras.backend.clear_session()
                model = architectures.getModelGivenModelOptionsAndWeightInits(
                    Arguments(), parameters(model_id)
                )
                outputs.append(model.predict(inputs, verbose=0))
            for first, second in zip(outputs[0], outputs[1]):
                np.testing.assert_array_equal(first, second)

    def test_prediction_transforms_and_reverse_complement(self) -> None:
        values = prediction.one_hot_dna(["ACGT", "TTAG"])
        restored = prediction.reverse_complement_one_hot(
            prediction.reverse_complement_one_hot(values)
        )
        np.testing.assert_array_equal(restored, values)
        probabilities = prediction.softmax(
            np.asarray([[0.0, 1.0, -1.0], [1000.0, 999.0, 998.0]])
        )
        np.testing.assert_allclose(
            probabilities.sum(axis=1), np.ones(2), rtol=0, atol=1.0e-7
        )
        np.testing.assert_allclose(
            prediction.mass_from_log1p(np.asarray([0.0, np.log(3.0)])),
            np.asarray([0.0, 2.0]),
            rtol=0,
            atol=1.0e-12,
        )

    def test_config_rejects_held_atac_exposure(self) -> None:
        raw = json.loads(CONFIG.read_text(encoding="utf-8"))
        raw["shared_contract"]["held_donor_atac_available_to_model"] = True
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "config.json"
            path.write_text(json.dumps(raw), encoding="utf-8")
            with self.assertRaises(training.SequenceControlTrainingError):
                training._load_config(path, "sequence_cnn_control")

    def test_geometry_and_ambiguous_bases_fail_closed(self) -> None:
        invalid = parameters("sequence_cnn_control")
        invalid["inputlen"] = "2112"
        with self.assertRaises(architectures.SequenceControlArchitectureError):
            architectures.getModelGivenModelOptionsAndWeightInits(
                Arguments(), invalid
            )
        with self.assertRaises(prediction.SequenceControlPredictionError):
            prediction.one_hot_dna(["ACNT"])

    def test_train_only_preparation_matches_native_fixture(self) -> None:
        fixture_root = (
            ROOT
            / "executions"
            / "chrombpnet-train-only-hyperparams-probe-21064020"
        )
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "prepared"
            summary = preparation.prepare(
                bigwig_path=fixture_root / "fixture" / "fixture.tn5.bw",
                peaks_path=fixture_root / "fixture" / "peaks.bed",
                nonpeaks_path=fixture_root / "fixture" / "nonpeaks.bed",
                fold_path=fixture_root / "fixture" / "actual.fold.json",
                output=output,
                seed=20260824,
                negative_sampling_ratio=0.1,
                outlier_threshold=0.9999,
                input_length=2114,
                output_length=1000,
                max_jitter=500,
            )
            self.assertEqual(summary["counts_loss_weight"], 210.0)
            self.assertFalse(
                summary["validation_outcomes_used_for_hyperparameter_fit"]
            )
            self.assertFalse(summary["test_outcomes_used_for_hyperparameter_fit"])
            self.assertFalse(summary["bias_model_loaded"])
            self.assertEqual(
                (output / "sequence_control.filtered.peaks.bed").read_bytes(),
                (
                    fixture_root
                    / "prepared"
                    / "train_only_filtered.peaks.bed"
                ).read_bytes(),
            )
            self.assertEqual(
                (output / "sequence_control.filtered.nonpeaks.bed").read_bytes(),
                (
                    fixture_root
                    / "prepared"
                    / "train_only_filtered.nonpeaks.bed"
                ).read_bytes(),
            )


if __name__ == "__main__":
    unittest.main()

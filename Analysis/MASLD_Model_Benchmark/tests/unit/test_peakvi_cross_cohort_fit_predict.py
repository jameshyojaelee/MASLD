from __future__ import annotations

import tempfile
from pathlib import Path
import unittest

import numpy as np

from scripts.peakvi_cross_cohort_fit_predict import (
    LINEAGES,
    PeakVICrossCohortError,
    donor_lineage_predictions,
    feature_rows,
    make_prediction_manifest,
)


class _FakeTensor:
    def __init__(self, value: np.ndarray) -> None:
        self.value = value

    def cuda(self) -> "_FakeTensor":
        return self

    def __getitem__(self, item: object) -> "_FakeTensor":
        return _FakeTensor(self.value[item])


class _FakeOutput:
    def __init__(self, value: np.ndarray) -> None:
        self.value = value

    def cpu(self) -> "_FakeOutput":
        return self

    def numpy(self) -> np.ndarray:
        return self.value


class _FakeTorch:
    class inference_mode:
        def __enter__(self) -> None:
            return None

        def __exit__(self, *_args: object) -> None:
            return None

    @staticmethod
    def from_numpy(value: np.ndarray) -> _FakeTensor:
        return _FakeTensor(value)

    @staticmethod
    def sigmoid(value: _FakeTensor) -> _FakeOutput:
        return _FakeOutput(1.0 / (1.0 + np.exp(-value.value)))


class _FakeDecoder:
    out_features = 16_000

    def eval(self) -> None:
        return None

    def __call__(self, value: _FakeTensor) -> _FakeTensor:
        return _FakeTensor(np.repeat(value.value[:, :1], self.out_features, axis=1))


class PeakVICrossCohortFitPredictTests(unittest.TestCase):
    def test_donor_lineage_aggregation_preserves_unit_axis(self) -> None:
        rows = []
        values = []
        for donor in range(1, 13):
            for lineage in LINEAGES:
                for replicate in range(2):
                    rows.append({"donor_id": f"Z{donor:02d}", "lineage_id": lineage})
                    values.append([float(donor + replicate)])
        result = donor_lineage_predictions(
            _FakeTorch(),
            _FakeDecoder(),
            np.asarray(values, dtype=np.float32),
            rows,
            batch_size=7,
        )
        self.assertEqual(result.shape, (12, 4, 16_000))
        self.assertTrue(np.all(result > 0))
        self.assertTrue(np.all(result < 1))

    def test_feature_axis_rejects_non_kilobase_window(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "features.tsv"
            with path.open("w", encoding="utf-8") as handle:
                handle.write("feature_index\twindow_id\tcontig\tstart\tend\n")
                for index in range(16_000):
                    end = index * 1000 + (999 if index == 9 else 1000)
                    handle.write(f"{index}\tw{index}\tchr1\t{index * 1000}\t{end}\n")
            with self.assertRaises(PeakVICrossCohortError):
                feature_rows(path)

    def test_manifest_keeps_native_and_projected_outputs_separate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            model = root / "model.json"
            model.write_text("{}\n", encoding="utf-8")
            predictions = root / "predictions.npy"
            missing = root / "missing.npy"
            np.save(predictions, np.ones((12, 4, 16_000), dtype=np.float32))
            np.save(missing, np.zeros((12, 4, 16_000), dtype=np.uint8))
            manifest = make_prediction_manifest(
                config={
                    "exchange_axis": {"artifacts_sha256": "a" * 64},
                    "seeds": [1103, 2207, 3301],
                },
                config_sha256="b" * 64,
                model_manifest=model,
                rotation_id="valid_context_test_target",
                predictions=predictions,
                missing_state=missing,
                output=root,
            )
            self.assertEqual(
                manifest["native_family_contract"]["native_output_family"],
                "joint_representation",
            )
            self.assertEqual(manifest["output_family"], "masked_accessibility_count")
            self.assertFalse(manifest["evidence_contract"]["target_role_atac_consumed"])
            self.assertFalse(manifest["metrics_calculated"])


if __name__ == "__main__":
    unittest.main()

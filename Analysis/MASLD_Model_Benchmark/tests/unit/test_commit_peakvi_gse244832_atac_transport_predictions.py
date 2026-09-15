from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from masld_bench.artifacts import freeze_tree
from scripts.commit_gse244832_atac_transport_predictions import (
    ATACPredictionCommitError,
    _digest,
)
from scripts.commit_peakvi_gse244832_atac_transport_predictions import (
    SEEDS,
    validate_peakvi_source_reuse,
)


class CommitPeakVIGSE244832ATACTransportPredictionsTests(unittest.TestCase):
    def _fixture(self, root: Path) -> tuple[Path, str, Path]:
        rotation = "valid_context_test_target"
        source_root = root / "source"
        source_rotation = source_root / "predictions" / rotation
        source_seed_records = []
        local_seed_records = []
        for seed in SEEDS:
            seed_root = source_rotation / f"seed-{seed}"
            seed_root.mkdir(parents=True)
            peakvi = seed_root / "peakvi.state_dict.pt"
            decoder = seed_root / "decoder.state_dict.pt"
            peakvi.write_bytes(f"peakvi-{seed}".encode())
            decoder.write_bytes(f"decoder-{seed}".encode())
            peakvi_sha = _digest(peakvi)
            decoder_sha = _digest(decoder)
            source_seed_records.append(
                {
                    "seed": seed,
                    "peakvi_state_file_sha256": peakvi_sha,
                    "decoder_state_file_sha256": decoder_sha,
                    "peakvi_parameter_sha256": f"peakvi-parameter-{seed}",
                    "decoder_parameter_sha256": f"decoder-parameter-{seed}",
                }
            )
            local_seed_records.append(
                {
                    "seed": seed,
                    "source_peakvi_state_path": str(peakvi.resolve()),
                    "source_peakvi_state_sha256": peakvi_sha,
                    "source_decoder_state_path": str(decoder.resolve()),
                    "source_decoder_state_sha256": decoder_sha,
                    "source_peakvi_parameter_sha256": f"peakvi-parameter-{seed}",
                    "source_decoder_parameter_sha256": f"decoder-parameter-{seed}",
                    "strict_state_restore_passed": True,
                    "target_encoding_repeat_bit_identical": True,
                    "fit_or_refit_performed": False,
                }
            )
        (source_rotation / "fit_receipt.json").write_text(
            json.dumps(
                {
                    "status": "passed",
                    "rotation_id": rotation,
                    "target_query_adaptation_performed": False,
                    "metrics_calculated": False,
                    "seeds": source_seed_records,
                }
            )
            + "\n",
            encoding="utf-8",
        )
        source_sha = freeze_tree(
            source_root,
            {
                "artifact_class": "peakvi_cross_cohort_fit_predict_and_blind_commit",
                "model_id": "peakvi_source_fit_fixed_window_decoder",
                "source_dataset_id": "gse296875",
                "target_dataset_id": "gse281367",
                "development_outcomes_read": False,
                "metrics_calculated": False,
                "status": "passed",
            },
        )

        prediction_root = root / "prediction"
        prediction_root.mkdir()
        model_path = prediction_root / "model_execution_manifest.json"
        model_path.write_text(
            json.dumps(
                {
                    "schema_version": "masld-bench-peakvi-gse244832-model-execution-v1",
                    "model_id": "peakvi",
                    "rotation_id": rotation,
                    "source_dataset_id": "gse296875",
                    "target_dataset_id": "gse244832",
                    "source_fit_job_id": "21101210",
                    "source_fit_result_path": str(source_root.resolve()),
                    "source_fit_result_artifacts_sha256": source_sha,
                    "fit_or_refit_performed": False,
                    "target_dataset_fit_or_adaptation": False,
                    "registered_output_family": "masked_accessibility_count",
                    "prediction_scale": (
                        "mean_predicted_single_nucleus_binary_accessibility_probability"
                    ),
                    "prediction_is_observed_fragment_count": False,
                    "seeds": local_seed_records,
                }
            )
            + "\n",
            encoding="utf-8",
        )
        (prediction_root / "prediction_bundle.json").write_text(
            json.dumps(
                {
                    "rotation_id": rotation,
                    "model_manifest": {
                        "path": "model_execution_manifest.json",
                        "sha256": _digest(model_path),
                    },
                }
            )
            + "\n",
            encoding="utf-8",
        )
        return source_root, source_sha, prediction_root

    def test_exact_frozen_source_state_reuse_passes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source, source_sha, prediction = self._fixture(Path(directory))
            receipt = validate_peakvi_source_reuse(prediction, source, source_sha)
            self.assertFalse(receipt["fit_or_refit_performed"])
            self.assertEqual(len(receipt["state_hashes"]), 3)

    def test_target_refit_declaration_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source, source_sha, prediction = self._fixture(Path(directory))
            model_path = prediction / "model_execution_manifest.json"
            model = json.loads(model_path.read_text(encoding="utf-8"))
            model["fit_or_refit_performed"] = True
            model_path.write_text(json.dumps(model) + "\n", encoding="utf-8")
            bundle_path = prediction / "prediction_bundle.json"
            bundle = json.loads(bundle_path.read_text(encoding="utf-8"))
            bundle["model_manifest"]["sha256"] = _digest(model_path)
            bundle_path.write_text(json.dumps(bundle) + "\n", encoding="utf-8")
            with self.assertRaises(ATACPredictionCommitError):
                validate_peakvi_source_reuse(prediction, source, source_sha)

    def test_source_artifact_hash_mismatch_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source, _, prediction = self._fixture(Path(directory))
            with self.assertRaises(ATACPredictionCommitError):
                validate_peakvi_source_reuse(prediction, source, "0" * 64)


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from masld_bench.gse244832_task_native_baselines import fit, predict
from scripts.commit_gse244832_atac_transport_predictions import validate_array_contract
from scripts.run_gse244832_task_native_blind_baselines import main, write_prediction


class RunGSE244832TaskNativeBlindBaselinesTests(unittest.TestCase):
    def test_cli_dispatches_config_as_config_path(self) -> None:
        with patch(
            "scripts.run_gse244832_task_native_blind_baselines.run",
            return_value={"status": "fixture"},
        ) as mocked, patch(
            "sys.argv",
            [
                "run",
                "--root",
                "/tmp/root",
                "--config",
                "/tmp/config.json",
                "--source-input",
                "/tmp/source",
                "--source-artifacts-sha256",
                "1" * 64,
                "--output",
                "/tmp/output",
            ],
        ):
            main()
        self.assertEqual(mocked.call_args.kwargs["config_path"], Path("/tmp/config.json"))
        self.assertNotIn("config", mocked.call_args.kwargs)

    def test_missing_mask_is_applied_before_prediction_and_written_as_nan(self) -> None:
        rng = np.random.default_rng(20260825)
        source_lineages = np.tile(np.arange(4), 8)
        source_context = rng.poisson(1.0, size=(32, 40)).astype(np.float64)
        source_target = rng.poisson(1.0, size=(32, 40)).astype(np.float64)
        source_context[:, 0] += 1
        source_target[:, 0] += 1
        state = fit(
            "observed_atac_glm",
            context_counts=source_context,
            target_counts=source_target,
            lineage_indices=source_lineages,
        )
        unit_states = np.zeros((18, 4), dtype=np.uint8)
        unit_states.reshape(-1)[48:52] = 1
        unit_states.reshape(-1)[52:] = 3
        eligible = unit_states.reshape(-1) == 0
        query = rng.poisson(1.0, size=(48, 40)).astype(np.float64)
        query[:, 0] += 1
        query_lineages = np.tile(np.arange(4), 12)
        small = predict(
            state,
            context_counts=query,
            lineage_indices=query_lineages,
        )
        expanded = np.tile(small, (1, 400))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "prediction"
            config = {
                "execution_registration": {"artifacts_sha256": "1" * 64},
                "exchange_axis": {"artifacts_sha256": "2" * 64},
                "query_authorities": {
                    "valid": {"artifacts_sha256": "3" * 64},
                    "test": {"artifacts_sha256": "4" * 64},
                },
            }
            paths = {"query_valid": root / "query_valid", "query_test": root / "query_test"}
            artifact_sha = write_prediction(
                output=output,
                model_id="observed_atac_glm",
                rotation_id="valid_context_test_target",
                state=state,
                predicted_eligible=expanded,
                eligible=eligible,
                unit_states=unit_states,
                config=config,
                paths=paths,
                source=root / "source",
                source_artifacts_sha256="5" * 64,
            )
            self.assertEqual(len(artifact_sha), 64)
            predictions = np.load(output / "predictions.float32.npy", allow_pickle=False)
            missing = np.load(output / "missing_state.uint8.npy", allow_pickle=False)
            receipt = validate_array_contract(
                predictions,
                missing,
                unit_states,
                output_family="masked_accessibility_count",
            )
            self.assertEqual(receipt["eligible_units"], 48)
            self.assertTrue(np.isnan(predictions.reshape(72, 16_000)[~eligible]).all())
            manifest = json.loads(
                (output / "model_execution_manifest.json").read_text(encoding="utf-8")
            )
            self.assertTrue(manifest["mask_applied_before_query_transform"])
            self.assertFalse(manifest["ineligible_query_units_read_by_transform"])


if __name__ == "__main__":
    unittest.main()

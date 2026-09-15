from __future__ import annotations

import unittest
from pathlib import Path
from unittest.mock import patch

from scripts.verify_gse244832_task_native_blind_baselines import (
    BlindBaselineVerificationError,
    _index,
    main,
)


class VerifyGSE244832TaskNativeBlindBaselinesTests(unittest.TestCase):
    def test_cli_dispatches_config_as_config_path(self) -> None:
        with patch(
            "scripts.verify_gse244832_task_native_blind_baselines.verify",
            return_value={"status": "fixture"},
        ) as mocked, patch(
            "sys.argv",
            [
                "verify",
                "--root",
                "/tmp/root",
                "--config",
                "/tmp/config.json",
                "--source-input",
                "/tmp/source",
                "--source-artifacts-sha256",
                "1" * 64,
                "--run-input",
                "/tmp/run",
                "--run-artifacts-sha256",
                "2" * 64,
                "--output",
                "/tmp/output",
            ],
        ):
            main()
        self.assertEqual(mocked.call_args.kwargs["config_path"], Path("/tmp/config.json"))
        self.assertNotIn("config", mocked.call_args.kwargs)

    def test_child_index_rejects_duplicate_model_rotation(self) -> None:
        records = [
            {"model_id": "lsi", "rotation_id": "valid_context_test_target"},
            {"model_id": "lsi", "rotation_id": "valid_context_test_target"},
        ]
        with self.assertRaises(BlindBaselineVerificationError):
            _index(records)

    def test_child_index_keys_model_and_rotation(self) -> None:
        records = [
            {"model_id": "lsi", "rotation_id": "valid_context_test_target"},
            {"model_id": "lsi", "rotation_id": "test_context_valid_target"},
        ]
        observed = _index(records)
        self.assertEqual(len(observed), 2)


if __name__ == "__main__":
    unittest.main()

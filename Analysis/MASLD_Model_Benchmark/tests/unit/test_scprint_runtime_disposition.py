from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tomllib
import unittest


ROOT = Path(__file__).resolve().parents[2]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ScprintRuntimeDispositionTests(unittest.TestCase):
    def test_terminal_disposition_is_bound_to_failed_strict_restore(self) -> None:
        disposition = json.loads(
            (ROOT / "config/artifacts/models/scprint/runtime_disposition.json").read_text()
        )
        self.assertEqual(disposition["model_id"], "scprint_v1_5_medium")
        self.assertEqual(disposition["disposition"], "terminal_runtime_incompatibility")
        self.assertFalse(disposition["strict_state_restoration_passed"])
        self.assertFalse(disposition["outcomes_or_labels_read"])

        artifact = ROOT / disposition["artifact"]["path"]
        self.assertEqual(
            digest(artifact / "ARTIFACTS.json"),
            disposition["artifact"]["artifacts_sha256"],
        )
        self.assertEqual(
            digest(artifact / "runtime_preflight.json"),
            disposition["artifact"]["runtime_preflight_sha256"],
        )
        manifest = json.loads((artifact / "ARTIFACTS.json").read_text())
        report = json.loads((artifact / "runtime_preflight.json").read_text())
        self.assertEqual(manifest["status"], "failed_strict_restore")
        self.assertEqual(report["status"], "failed_strict_restore")
        self.assertFalse(report["outcomes_or_labels_read"])
        self.assertIn("compress_class_dim", report["architecture"]["restoration_error"])

    def test_model_registry_is_blocked_without_checkpoint_substitution(self) -> None:
        registry = tomllib.loads(
            (ROOT / "config/models/cell_foundation.toml").read_text()
        )
        model = next(
            item
            for item in registry["models"]
            if item["model_id"] == "scprint_v1_5_medium"
        )
        self.assertEqual(model["status"], "deferred")
        self.assertTrue(model["admission_blocking"])
        self.assertEqual(
            model["checkpoint_sha256"],
            "a4cf0753270d4ff451a5dbddadd4c59a53e10e4c3e96a8af0db407ad893c36c5",
        )
        self.assertIn("Do not drop compress_class_dim", model["blockers"][0])


if __name__ == "__main__":
    unittest.main()

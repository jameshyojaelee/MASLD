#!/usr/bin/env python3
"""Unit checks for the Enformer and Sei fail-closed lane audit."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

from masld_bench.artifacts import freeze_tree


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts/audit_enformer_sei_native_lanes.py"
SPEC = importlib.util.spec_from_file_location("enformer_sei_lane_audit", MODULE)
assert SPEC and SPEC.loader
lane_audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(lane_audit)


class EnformerSeiLaneAuditTest(unittest.TestCase):
    @staticmethod
    def _frozen_fixture(root: Path) -> str:
        (root / "payload.json").write_text('{"status":"pass"}\n', encoding="utf-8")
        freeze_tree(root, {"fixture": True})
        return hashlib.sha256((root / "ARTIFACTS.json").read_bytes()).hexdigest()

    def test_require_artifact_accepts_recursive_frozen_tree(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            expected = self._frozen_fixture(root)
            self.assertEqual(
                lane_audit.require_artifact(root, expected)["metadata"],
                {"fixture": True},
            )

    def test_require_artifact_rejects_tampered_payload(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            expected = self._frozen_fixture(root)
            (root / "payload.json").write_text('{"status":"changed"}\n', encoding="utf-8")
            with self.assertRaises(lane_audit.LaneAuditError):
                lane_audit.require_artifact(root, expected)

    def test_require_authority_rejects_digest_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "authority.json").write_text('{"status":"pass"}\n', encoding="utf-8")
            with self.assertRaises(lane_audit.LaneAuditError):
                lane_audit.require_authority(root, "authority.json", "0" * 64)

    def test_model_registry_rejects_duplicate_ids(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "config/models"
            path.mkdir(parents=True)
            (path / "regulatory_sequence.toml").write_text(
                "[[models]]\nmodel_id='sei'\n[[models]]\nmodel_id='sei'\n",
                encoding="utf-8",
            )
            with self.assertRaises(lane_audit.LaneAuditError):
                lane_audit.model_registry(root)

    def test_require_fails_closed(self) -> None:
        with self.assertRaisesRegex(lane_audit.LaneAuditError, "expected failure"):
            lane_audit.require(False, "expected failure")


if __name__ == "__main__":
    unittest.main()

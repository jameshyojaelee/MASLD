#!/usr/bin/env python3
"""Unit checks for the fail-closed EpiBERT/EPCOTv2 lane audit."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

from masld_bench.artifacts import freeze_tree


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts/audit_epibert_epcotv2_lane.py"
SPEC = importlib.util.spec_from_file_location("lane_audit", MODULE)
assert SPEC and SPEC.loader
lane_audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(lane_audit)


class LaneAuditTest(unittest.TestCase):
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

    def test_require_artifact_rejects_missing_or_modified_complete(self) -> None:
        for missing in (True, False):
            with self.subTest(missing=missing), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                expected = self._frozen_fixture(root)
                complete = root / "COMPLETE"
                if missing:
                    complete.unlink()
                else:
                    value = json.loads(complete.read_text(encoding="utf-8"))
                    value["artifact_count"] += 1
                    complete.write_text(json.dumps(value) + "\n", encoding="utf-8")
                with self.assertRaises(lane_audit.LaneAuditError):
                    lane_audit.require_artifact(root, expected)

    def test_require_artifact_rejects_missing_manifest_member(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            expected = self._frozen_fixture(root)
            (root / "payload.json").unlink()
            with self.assertRaises(lane_audit.LaneAuditError):
                lane_audit.require_artifact(root, expected)

    def test_license_absence_requires_exact_revision(self) -> None:
        metadata = {
            "sha": "abc",
            "cardData": {"title": "model"},
            "siblings": [{"rfilename": "README.md"}],
        }
        self.assertTrue(lane_audit.api_license_absent(metadata, "abc"))
        metadata["cardData"]["license"] = "mit"
        self.assertFalse(lane_audit.api_license_absent(metadata, "abc"))
        with self.assertRaises(lane_audit.LaneAuditError):
            lane_audit.api_license_absent(metadata, "def")

    def test_license_file_prevents_absence_finding(self) -> None:
        metadata = {
            "sha": "abc",
            "cardData": {},
            "siblings": [{"rfilename": "LICENSE.txt"}],
        }
        self.assertFalse(lane_audit.api_license_absent(metadata, "abc"))


if __name__ == "__main__":
    unittest.main()

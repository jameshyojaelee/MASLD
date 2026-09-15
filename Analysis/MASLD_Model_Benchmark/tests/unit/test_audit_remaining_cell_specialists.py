#!/usr/bin/env python3
"""Unit checks for remaining cell-specialist reconciliation."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

from masld_bench.artifacts import freeze_tree


ROOT = Path(__file__).parents[2]
MODULE = ROOT / "scripts/audit_remaining_cell_specialists.py"
SPEC = importlib.util.spec_from_file_location("remaining_cell_specialists", MODULE)
assert SPEC and SPEC.loader
cell_audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(cell_audit)


class RemainingCellSpecialistAuditTest(unittest.TestCase):
    @staticmethod
    def _census() -> dict:
        return json.loads((ROOT / cell_audit.CENSUS_RELATIVE).read_text(encoding="utf-8"))

    def test_validate_census_accepts_single_direct_wrapper(self) -> None:
        expected, rows = cell_audit.validate_census(self._census())
        self.assertEqual(set(rows), expected)
        self.assertTrue(rows["scimilarity_v1_1"]["direct_preregistered_wrapper_ready"])

    def test_validate_census_rejects_second_ready_model(self) -> None:
        census = copy.deepcopy(self._census())
        census["models"][0]["direct_preregistered_wrapper_ready"] = True
        with self.assertRaisesRegex(
            cell_audit.CellSpecialistAuditError, "readiness differs"
        ):
            cell_audit.validate_census(census)

    def test_validate_census_rejects_scheduling_authority(self) -> None:
        census = copy.deepcopy(self._census())
        census["models"][0]["schedulable"] = True
        with self.assertRaisesRegex(
            cell_audit.CellSpecialistAuditError, "authorized scheduling"
        ):
            cell_audit.validate_census(census)

    def test_validate_census_rejects_query_label_refinement(self) -> None:
        census = copy.deepcopy(self._census())
        census["common_lane_contract"]["query_label_refinement_prohibited"] = False
        with self.assertRaisesRegex(
            cell_audit.CellSpecialistAuditError, "leakage boundary"
        ):
            cell_audit.validate_census(census)

    def test_model_rows_rejects_duplicate_ids(self) -> None:
        value = {"models": [{"model_id": "scbert"}, {"model_id": "scbert"}]}
        with self.assertRaises(cell_audit.CellSpecialistAuditError):
            cell_audit.model_rows(value)

    def test_require_artifact_rejects_tampered_small_tree(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            artifact = root / "artifact"
            artifact.mkdir()
            (artifact / "payload.json").write_text(
                '{"status":"pass"}\n', encoding="utf-8"
            )
            freeze_tree(artifact, {"fixture": True})
            expected = hashlib.sha256(
                (artifact / "ARTIFACTS.json").read_bytes()
            ).hexdigest()
            original = set(cell_audit.DEEP_VERIFY)
            cell_audit.DEEP_VERIFY.add("artifact")
            try:
                (artifact / "payload.json").write_text(
                    '{"status":"changed"}\n', encoding="utf-8"
                )
                with self.assertRaises(cell_audit.CellSpecialistAuditError):
                    cell_audit.require_artifact(root, "artifact", expected)
            finally:
                cell_audit.DEEP_VERIFY.clear()
                cell_audit.DEEP_VERIFY.update(original)

    def test_legacy_checksum_tree_rejects_tampered_payload(self) -> None:
        relative = "legacy"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            artifact = root / relative
            artifact.mkdir()
            payloads = {
                "ARTIFACTS.json": b'{"status":"failed"}\n',
                "runtime_preflight.json": b'{"status":"failed"}\n',
            }
            for name, value in payloads.items():
                (artifact / name).write_bytes(value)
            checksums = {
                name: hashlib.sha256(value).hexdigest()
                for name, value in payloads.items()
            }
            (artifact / "SHA256SUMS").write_text(
                "".join(f"{checksum}  {name}\n" for name, checksum in checksums.items()),
                encoding="utf-8",
            )
            original = dict(cell_audit.LEGACY_CHECKSUM_VERIFY)
            cell_audit.LEGACY_CHECKSUM_VERIFY[relative] = checksums
            try:
                expected = checksums["ARTIFACTS.json"]
                self.assertEqual(
                    cell_audit.require_artifact(root, relative, expected)["status"],
                    "failed",
                )
                (artifact / "runtime_preflight.json").write_text(
                    '{"status":"changed"}\n', encoding="utf-8"
                )
                with self.assertRaises(cell_audit.CellSpecialistAuditError):
                    cell_audit.require_artifact(root, relative, expected)
            finally:
                cell_audit.LEGACY_CHECKSUM_VERIFY.clear()
                cell_audit.LEGACY_CHECKSUM_VERIFY.update(original)


if __name__ == "__main__":
    unittest.main()

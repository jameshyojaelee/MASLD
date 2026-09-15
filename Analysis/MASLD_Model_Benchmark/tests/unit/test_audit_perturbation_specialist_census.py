#!/usr/bin/env python3
"""Unit checks for the perturbation specialist census audit."""

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
MODULE = ROOT / "scripts/audit_perturbation_specialist_census.py"
SPEC = importlib.util.spec_from_file_location("perturbation_specialist_census", MODULE)
assert SPEC and SPEC.loader
census_audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(census_audit)


class PerturbationSpecialistCensusAuditTest(unittest.TestCase):
    @staticmethod
    def _frozen_fixture(root: Path) -> str:
        (root / "payload.json").write_text('{"status":"pass"}\n', encoding="utf-8")
        freeze_tree(root, {"fixture": True})
        return hashlib.sha256((root / "ARTIFACTS.json").read_bytes()).hexdigest()

    @staticmethod
    def _census() -> dict:
        return json.loads((ROOT / census_audit.CENSUS_RELATIVE).read_text(encoding="utf-8"))

    def test_require_artifact_accepts_recursive_frozen_tree(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            expected = self._frozen_fixture(root)
            self.assertEqual(
                census_audit.require_artifact(root.parent, root.name, expected)["metadata"],
                {"fixture": True},
            )

    def test_require_artifact_rejects_tampered_payload(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            expected = self._frozen_fixture(root)
            (root / "payload.json").write_text('{"status":"changed"}\n', encoding="utf-8")
            with self.assertRaises(census_audit.CensusAuditError):
                census_audit.require_artifact(root.parent, root.name, expected)

    def test_validate_census_accepts_frozen_contract(self) -> None:
        expected, rows = census_audit.validate_census(self._census())
        self.assertEqual(set(rows), expected)

    def test_validate_census_rejects_silent_scheduling(self) -> None:
        census = copy.deepcopy(self._census())
        census["models"][0]["schedulable"] = True
        with self.assertRaisesRegex(
            census_audit.CensusAuditError, "became schedulable"
        ):
            census_audit.validate_census(census)

    def test_validate_census_rejects_cells_as_replicates(self) -> None:
        census = copy.deepcopy(self._census())
        census["global_contract"]["not_independent_replicates"].remove("cells")
        with self.assertRaisesRegex(
            census_audit.CensusAuditError, "pseudoreplication"
        ):
            census_audit.validate_census(census)

    def test_model_rows_rejects_duplicate_ids(self) -> None:
        config = {"models": [{"model_id": "gears"}, {"model_id": "gears"}]}
        with self.assertRaises(census_audit.CensusAuditError):
            census_audit.model_rows(config)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Checks for the content-free GSE268273 evaluator-receipt incident."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from masld_bench.artifacts import verify_frozen_tree
from scripts import validate_gse268273_evaluator_receipt_incident as incident


ROOT = Path(__file__).parents[2]
SOURCE = (
    ROOT
    / "config/artifacts/incidents/gse268273_evaluator_receipt_exposure_20260824.json"
)


class GSE268273EvaluatorReceiptIncidentTests(unittest.TestCase):
    def test_source_has_only_content_free_incident_fields(self) -> None:
        value = json.loads(SOURCE.read_text(encoding="utf-8"))
        self.assertEqual(value, incident.EXPECTED)
        self.assertEqual(set(value), set(incident.EXPECTED))
        self.assertTrue(value["outcomes.tsv_not_opened"])
        self.assertTrue(value["labels_not_joined"])
        self.assertEqual(value["disposition"], "evaluator_receipt_exposed")

    def test_freezes_fail_closed_incident(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "incident"
            incident.run(source=SOURCE, output=output)
            manifest = verify_frozen_tree(output)
            self.assertFalse(
                manifest["metadata"]["clean_or_sealed_champion_eligible"]
            )
            self.assertTrue(
                manifest["metadata"][
                    "independent_unexposed_rederivation_or_audit_required"
                ]
            )
            self.assertEqual(
                json.loads((output / "incident_receipt.json").read_text()),
                incident.EXPECTED,
            )

    def test_rejects_any_additional_content(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.json"
            value = {**incident.EXPECTED, "extra": "forbidden"}
            source.write_text(json.dumps(value), encoding="utf-8")
            with self.assertRaisesRegex(incident.IncidentFirewallError, "fields"):
                incident.run(source=source, output=root / "rejected")


if __name__ == "__main__":
    unittest.main()

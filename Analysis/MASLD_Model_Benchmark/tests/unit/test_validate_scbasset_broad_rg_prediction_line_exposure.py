from __future__ import annotations

from pathlib import Path
import unittest

from scripts.validate_scbasset_broad_rg_prediction_line_exposure import (
    validate_incident,
)


ROOT = Path(__file__).resolve().parents[2]
INCIDENT = (
    ROOT
    / "config/artifacts/incidents/scbasset_broad_rg_prediction_line_exposure_20260825.json"
)


class ScBassetBroadRgPredictionLineExposureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.receipt = validate_incident(INCIDENT)

    def test_unrelated_prediction_line_exposure_is_explicit(self) -> None:
        self.assertTrue(self.receipt["unrelated_prediction_tsv_lines_emitted"])

    def test_scbasset_outcomes_metrics_and_sealed_data_are_excluded(self) -> None:
        self.assertFalse(self.receipt["scbasset_prediction_values_emitted"])
        self.assertFalse(self.receipt["scbasset_metrics_emitted"])
        self.assertFalse(self.receipt["raw_outcomes_emitted"])
        self.assertFalse(self.receipt["sealed_data_emitted"])

    def test_incident_values_are_unused_and_validator_does_not_reopen_predictions(self) -> None:
        self.assertFalse(self.receipt["incident_values_used"])
        self.assertFalse(self.receipt["referenced_prediction_files_reopened_by_validator"])


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

from scripts.audit_gse296875_observed_multiome_complete_rectangle_production_admission import ProductionAdmissionError, validate


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "config/campaigns/gse296875_observed_multiome_factorized_complete_rectangle_production_revision2_20260825.json"


class ObservedMultiomeCompleteRectangleProductionAdmissionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))

    def test_exact_complete_rectangle_is_authorized_without_array_deserialization(self) -> None:
        receipt = validate(ROOT, self.manifest)
        self.assertTrue(receipt["production_aggregation_authorized"])
        self.assertEqual(receipt["seed_bundles_verified"], 5)
        self.assertEqual(receipt["surface_seed_runs_verified"], 125)
        self.assertFalse(receipt["production_prediction_arrays_deserialized"])
        self.assertFalse(receipt["production_evaluator_arrays_deserialized"])

    def test_claim_boundary_drift_is_rejected(self) -> None:
        changed = deepcopy(self.manifest)
        changed["claim_boundary"]["external_claim_allowed"] = True
        with self.assertRaises(ProductionAdmissionError):
            validate(ROOT, changed)


if __name__ == "__main__":
    unittest.main()

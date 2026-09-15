from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

from scripts.audit_gse296875_observed_multiome_factorized_seed_bundle_admission import SeedBundleAdmissionError, validate


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse296875_observed_multiome_factorized_seed_bundle_admission_20260825.json"


class ObservedMultiomeFactorizedSeedBundleAdmissionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_exact_execution_surface_is_authorized_without_values(self) -> None:
        receipt = validate(ROOT, self.config)
        self.assertTrue(receipt["seed_bundle_execution_authorized"])
        self.assertEqual(receipt["seed_bundle_jobs_authorized"], 5)
        self.assertEqual(receipt["surface_seed_runs_authorized"], 125)
        self.assertFalse(receipt["evaluator_values_read"])
        self.assertFalse(receipt["prediction_values_read"])
        self.assertFalse(receipt["partial_ranking_authorized"])

    def test_command_drift_is_rejected(self) -> None:
        changed = deepcopy(self.config)
        changed["execution"]["submission_commands"][0] += " --extra"
        with self.assertRaises(SeedBundleAdmissionError):
            validate(ROOT, changed)

    def test_resource_drift_is_rejected(self) -> None:
        changed = deepcopy(self.config)
        changed["resources"]["cpus_per_job"] = 9
        with self.assertRaises(SeedBundleAdmissionError):
            validate(ROOT, changed)


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import unittest

from scripts.verify_gse296875_observed_multiome_model_inputs import (
    ModelInputVerificationError,
    validate_config,
    verify,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse296875_observed_multiome_model_input_verification_20260825.json"


class VerifyObservedMultiomeModelInputsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_real_artifact_passes_independent_verification(self) -> None:
        receipt = verify(ROOT, deepcopy(self.config))
        self.assertTrue(receipt["promotion_gate_passed"])
        self.assertFalse(receipt["target_atac_values_present_in_model_artifacts"])
        self.assertEqual(receipt["folds_verified"], 5)

    def test_artifact_hash_drift_fails_closed(self) -> None:
        mutated = deepcopy(self.config)
        mutated["model_input_artifact"]["artifacts_sha256"] = "0" * 64
        with self.assertRaises(ModelInputVerificationError):
            validate_config(ROOT, mutated)

    def test_firewall_mutation_fails_closed(self) -> None:
        mutated = deepcopy(self.config)
        mutated["firewall"]["model_fit"] = True
        with self.assertRaises(ModelInputVerificationError):
            validate_config(ROOT, mutated)


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

from copy import deepcopy
import inspect
import json
from pathlib import Path
import unittest

from scripts import verify_gse296875_observed_multiome_evaluator_outcomes as verifier


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse296875_observed_multiome_evaluator_verification_20260825.json"


class VerifyObservedMultiomeEvaluatorOutcomesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_real_evaluator_artifact_passes(self) -> None:
        receipt = verifier.verify(ROOT, deepcopy(self.config))
        self.assertTrue(receipt["promotion_gate_passed"])
        self.assertTrue(receipt["target_atac_count_values_read_by_evaluator_verifier"])
        self.assertFalse(receipt["model_input_biological_count_values_read"])
        self.assertEqual(len(receipt["folds_verified"]), 5)

    def test_firewall_and_hash_mutations_fail_closed(self) -> None:
        firewall = deepcopy(self.config)
        firewall["firewall"]["prediction_bundle_read"] = True
        with self.assertRaises(verifier.EvaluatorVerificationError):
            verifier.validate_config(ROOT, firewall)
        artifact = deepcopy(self.config)
        artifact["evaluator_artifact"]["artifacts_sha256"] = "0" * 64
        with self.assertRaises(verifier.EvaluatorVerificationError):
            verifier.validate_config(ROOT, artifact)

    def test_source_does_not_open_model_input_counts(self) -> None:
        source = inspect.getsource(verifier)
        self.assertNotIn('model["rna/', source)
        self.assertNotIn('model["observed_atac_input/counts_csr', source)
        self.assertNotIn('model["target_atac', source)
        self.assertNotIn("prediction", source.split("def verify", 1)[1].split("def main", 1)[0].replace('"prediction_bundle_read"', ""))


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

from copy import deepcopy
import inspect
import json
from pathlib import Path
import unittest

from scripts import verify_gse296875_observed_multiome_training_labels as verifier


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/campaigns/gse296875_observed_multiome_training_label_verification_20260825.json"


class VerifyObservedMultiomeTrainingLabelsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads(CONFIG.read_text(encoding="utf-8"))

    def test_real_training_rectangle_passes(self) -> None:
        receipt = verifier.verify(ROOT, deepcopy(self.config))
        self.assertTrue(receipt["promotion_gate_passed"])
        self.assertEqual(receipt["surface_count"], 25)
        self.assertFalse(receipt["held_donor_rows_present"])

    def test_firewall_and_hash_mutations_fail_closed(self) -> None:
        firewall = deepcopy(self.config)
        firewall["firewall"]["evaluator_artifact_read"] = True
        with self.assertRaises(verifier.TrainingLabelVerificationError):
            verifier.validate_config(ROOT, firewall)
        artifact = deepcopy(self.config)
        artifact["training_label_artifact"]["artifacts_sha256"] = "0" * 64
        with self.assertRaises(verifier.TrainingLabelVerificationError):
            verifier.validate_config(ROOT, artifact)

    def test_source_does_not_open_model_or_evaluator_counts(self) -> None:
        source = inspect.getsource(verifier)
        self.assertNotIn('model["rna/', source)
        self.assertNotIn('model["observed_atac_input/counts_csr', source)
        self.assertNotIn('model["targets_without_values/counts_csr', source)
        self.assertNotIn("evaluator_outcomes.h5", source)


if __name__ == "__main__":
    unittest.main()

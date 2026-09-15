from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest

from scripts.audit_scooby_epicardioids_runtime_schema_gate import (
    ScoobySchemaGateError,
    install_peft_inference_placeholder,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/scooby_epicardioids_runtime_schema_gate.json"


class ScoobyEpicardioidsRuntimeSchemaGateTests(unittest.TestCase):
    def test_observed_atac_topology_is_not_rna_conditioned_atac(self) -> None:
        config = json.loads(CONFIG.read_text(encoding="utf-8"))
        task = config["task_scope"]
        self.assertTrue(task["observed_atac_in_context"])
        self.assertEqual(task["rna_conditioned_atac"], "ineligible_observed_ATAC_contributed_to_context")
        self.assertFalse(task["liver_query_encoder_available"])
        self.assertFalse(task["sealed_inference_allowed"])
        self.assertFalse(task["universal_or_champion_claim_allowed"])

    def test_schema_gate_cannot_load_checkpoint_or_full_sequence(self) -> None:
        policy = json.loads(CONFIG.read_text(encoding="utf-8"))["gate_policy"]
        self.assertFalse(policy["checkpoint_payload_downloaded_in_this_gate"])
        self.assertFalse(policy["checkpoint_deserialized_in_this_gate"])
        self.assertFalse(policy["full_sequence_forward_allowed"])
        self.assertTrue(policy["random_weight_decoder_probe_allowed"])

    def test_peft_placeholder_blocks_training_calls(self) -> None:
        prior = sys.modules.pop("peft", None)
        try:
            install_peft_inference_placeholder()
            import peft

            with self.assertRaises(ScoobySchemaGateError):
                peft.LoraConfig()
            with self.assertRaises(ScoobySchemaGateError):
                peft.get_peft_model(object())
        finally:
            sys.modules.pop("peft", None)
            if prior is not None:
                sys.modules["peft"] = prior


if __name__ == "__main__":
    unittest.main()

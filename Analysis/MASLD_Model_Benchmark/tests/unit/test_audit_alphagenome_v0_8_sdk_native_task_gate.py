from __future__ import annotations

import json
from pathlib import Path
import unittest

from scripts.audit_alphagenome_v0_8_sdk_native_task_gate import (
    OUTPUT_ENUM,
    enum_members,
    message_body,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config/alphagenome_v0_8_sdk_native_task_gate.json"


class AlphaGenomeSdkNativeTaskGateTests(unittest.TestCase):
    def test_contract_keeps_every_execution_and_release_gate_closed(self) -> None:
        config = json.loads(CONFIG.read_text(encoding="utf-8"))
        disposition = config["disposition"]
        self.assertTrue(disposition["sdk_source_and_offline_runtime_probe_allowed"])
        self.assertFalse(disposition["api_key_inspection_allowed"])
        self.assertFalse(disposition["api_connection_allowed"])
        self.assertFalse(disposition["checkpoint_download_allowed"])
        self.assertFalse(disposition["model_forward_allowed"])
        self.assertFalse(disposition["outcome_access_allowed"])
        self.assertFalse(disposition["sealed_access_allowed"])
        self.assertFalse(disposition["open_champion_eligible"])
        self.assertFalse(disposition["conditional_open_model_training_allowed"])

    def test_task_scope_is_sequence_native_not_observed_multiome(self) -> None:
        task = json.loads(CONFIG.read_text(encoding="utf-8"))["task_scope"]
        self.assertFalse(task["donor_specific_prediction"])
        self.assertEqual(task["rna_conditioned_atac"], "unsupported_static_sequence_baseline_only")
        self.assertEqual(task["observed_multiome"], "unsupported_as_observed_assay_model")
        self.assertEqual(task["cell_state_mapping"], "unsupported")
        self.assertEqual(task["histology_or_clinical_prediction"], "unsupported")

    def test_proto_helpers_identify_version_asymmetry(self) -> None:
        proto = """
enum OutputType {
  OUTPUT_TYPE_UNSPECIFIED = 0;
  OUTPUT_TYPE_ATAC = 1;
  OUTPUT_TYPE_RNA_SEQ = 2;
}
message PredictSequenceRequest {
  string sequence = 1;
  string model_version = 5;
}
message PredictSequenceResponse {
  bytes output = 1;
}
message MetadataRequest {
  int64 organism = 1;
}
"""
        self.assertEqual(enum_members(proto, "OutputType", prefix="OUTPUT_TYPE_"), ("ATAC", "RNA_SEQ"))
        self.assertIn("model_version", message_body(proto, "PredictSequenceRequest"))
        self.assertNotIn("model_version", message_body(proto, "PredictSequenceResponse"))
        self.assertNotIn("model_version", message_body(proto, "MetadataRequest"))

    def test_native_output_contract_is_complete_and_ordered(self) -> None:
        config = json.loads(CONFIG.read_text(encoding="utf-8"))
        self.assertEqual(tuple(config["native_interface"]["output_types"]), OUTPUT_ENUM)
        self.assertEqual(config["native_interface"]["logical_version_for_any_future_probe"], "ALL_FOLDS")
        self.assertFalse(config["native_interface"]["default_model_version_allowed"])


if __name__ == "__main__":
    unittest.main()

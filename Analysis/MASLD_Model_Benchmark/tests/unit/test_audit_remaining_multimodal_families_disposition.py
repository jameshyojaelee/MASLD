from __future__ import annotations

from pathlib import Path
import unittest

from scripts.audit_remaining_multimodal_families_disposition import audit_disposition


ROOT = Path(__file__).resolve().parents[2]
AUTHORITY = (
    ROOT
    / "config/artifacts/models/remaining_multimodal_families/disposition_20260824.json"
)


class RemainingMultimodalFamiliesDispositionTests(unittest.TestCase):
    def test_all_native_output_classes_remain_fail_closed(self) -> None:
        receipt = audit_disposition(ROOT, AUTHORITY)
        self.assertEqual(
            receipt["status"],
            "pass_metadata_audit_all_new_execution_gates_closed",
        )
        self.assertEqual(receipt["verified_metadata_files"], 72)
        self.assertEqual(receipt["model_identity_count"], 13)
        self.assertEqual(receipt["runtime_manifest_count"], 5)
        self.assertEqual(receipt["execution_manifest_count"], 20)
        self.assertEqual(receipt["historical_source_deviation_count"], 3)
        self.assertEqual(receipt["shared_registry_drift_count"], 3)
        self.assertEqual(receipt["admitted_model_count"], 0)
        self.assertIs(receipt["training_authorized"], False)
        self.assertIs(receipt["prediction_authorized"], False)
        self.assertIs(receipt["scoring_authorized"], False)
        self.assertIs(receipt["sbatch_submission_authorized"], False)
        self.assertIs(receipt["sealed_outcomes_opened"], False)
        self.assertIs(receipt["development_outcomes_opened"], False)
        self.assertIs(receipt["model_payloads_opened"], False)
        self.assertIs(receipt["data_matrices_opened"], False)


if __name__ == "__main__":
    unittest.main()

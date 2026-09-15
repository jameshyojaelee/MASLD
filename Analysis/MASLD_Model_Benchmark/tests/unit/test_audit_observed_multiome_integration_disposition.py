from __future__ import annotations

from pathlib import Path
import unittest

from scripts.audit_observed_multiome_integration_disposition import audit_disposition


ROOT = Path(__file__).resolve().parents[2]
AUTHORITY = (
    ROOT
    / "config/artifacts/models/observed_multiome_integration/disposition_20260824.json"
)


class ObservedMultiomeIntegrationDispositionTests(unittest.TestCase):
    def test_all_four_model_lanes_remain_fail_closed(self) -> None:
        receipt = audit_disposition(ROOT, AUTHORITY)
        self.assertEqual(
            receipt["status"], "pass_metadata_audit_all_model_lanes_blocked"
        )
        self.assertEqual(receipt["verified_metadata_files"], 34)
        self.assertEqual(receipt["runtime_manifest_count"], 4)
        self.assertEqual(receipt["execution_manifest_count"], 11)
        self.assertEqual(receipt["historical_source_drift_count"], 1)
        self.assertEqual(receipt["shared_registry_drift_count"], 2)
        self.assertEqual(receipt["admitted_model_count"], 0)
        self.assertIs(receipt["task_spec_bound"], False)
        self.assertIs(receipt["external_seal_bound"], False)
        self.assertIs(receipt["training_authorized"], False)
        self.assertIs(receipt["scoring_authorized"], False)
        self.assertIs(receipt["gpu_sbatch_authorized"], False)
        self.assertIs(receipt["sealed_outcomes_opened"], False)
        self.assertIs(receipt["development_outcomes_opened"], False)
        self.assertIs(receipt["model_artifacts_opened"], False)
        self.assertIs(receipt["data_matrices_opened"], False)


if __name__ == "__main__":
    unittest.main()

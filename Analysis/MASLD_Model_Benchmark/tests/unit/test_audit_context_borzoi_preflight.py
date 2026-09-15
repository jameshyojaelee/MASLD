from __future__ import annotations

from pathlib import Path
import unittest

from scripts.audit_context_borzoi_preflight import audit_preflight


ROOT = Path(__file__).resolve().parents[2]
AUTHORITY = (
    ROOT
    / "config/artifacts/models/context_borzoi/outcome_blind_architecture_preflight_20260824.json"
)


class ContextBorzoiPreflightAuditTests(unittest.TestCase):
    def test_preflight_is_outcome_blind_and_both_activation_gates_are_closed(self) -> None:
        receipt = audit_preflight(ROOT, AUTHORITY)
        self.assertEqual(receipt["status"], "pass_outcome_blind_preflight_blocked")
        self.assertIs(receipt["complementarity_gate_passed"], False)
        self.assertIs(receipt["component_admission_passed"], False)
        self.assertIs(receipt["architecture_executable"], False)
        self.assertIs(receipt["training_authorized"], False)
        self.assertIs(receipt["gpu_sbatch_authorized"], False)
        self.assertEqual(receipt["ablation_count"], 5)
        self.assertIs(receipt["missing_as_zero_allowed"], False)
        self.assertIs(receipt["checkpoint_opened"], False)
        self.assertIs(receipt["development_outcomes_opened"], False)
        self.assertIs(receipt["sealed_data_opened"], False)
        self.assertIs(receipt["evaluator_outputs_opened"], False)
        self.assertIs(receipt["borzoi_weight_admitted"], False)
        self.assertIs(receipt["corgi_weight_admitted_as_component"], False)


if __name__ == "__main__":
    unittest.main()

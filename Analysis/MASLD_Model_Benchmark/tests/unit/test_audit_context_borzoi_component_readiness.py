from __future__ import annotations

from pathlib import Path
import unittest

from scripts.audit_context_borzoi_component_readiness import audit_readiness


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = (
    ROOT
    / "config/artifacts/models/context_borzoi/component_interface_readiness_20260825.json"
)


class ContextBorzoiComponentReadinessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.receipt = audit_readiness(ROOT, CONTRACT)

    def test_interface_is_reconciled_but_every_build_gate_remains_closed(self) -> None:
        receipt = self.receipt
        self.assertEqual(receipt["status"], "pass_interface_contract_build_blocked")
        self.assertTrue(receipt["component_interface_schema_ready"])
        self.assertFalse(receipt["complementarity_trigger_passed"])
        self.assertFalse(receipt["component_selection_locked"])
        self.assertFalse(receipt["architecture_built"])
        self.assertFalse(receipt["training_authorized"])
        self.assertFalse(receipt["prediction_authorized"])
        self.assertFalse(receipt["gpu_bundle_registration_authorized"])
        self.assertFalse(receipt["gpu_job_submission_authorized"])

    def test_provisional_geometry_and_family_native_input_boundary_are_explicit(self) -> None:
        receipt = self.receipt
        self.assertEqual(receipt["provisional_input_shape"], [None, 4, 524288])
        self.assertEqual(receipt["provisional_latent_shape"], [None, 1920, 6144])
        self.assertFalse(receipt["context_missing_as_zero_allowed"])
        self.assertFalse(receipt["query_atac_allowed_for_rna_conditioned_task"])

    def test_component_evidence_does_not_select_or_reuse_a_parent(self) -> None:
        receipt = self.receipt
        self.assertFalse(receipt["borzoi_selected"])
        self.assertFalse(receipt["corgi_selected_or_reused"])
        self.assertFalse(receipt["scooby_selected_or_reused"])
        self.assertFalse(receipt["checkpoint_opened_by_audit"])

    def test_outcome_firewall_is_closed(self) -> None:
        receipt = self.receipt
        self.assertFalse(receipt["biological_data_read"])
        self.assertFalse(receipt["development_outcomes_read"])
        self.assertFalse(receipt["sealed_features_read"])
        self.assertFalse(receipt["sealed_labels_or_outcomes_read"])


if __name__ == "__main__":
    unittest.main()

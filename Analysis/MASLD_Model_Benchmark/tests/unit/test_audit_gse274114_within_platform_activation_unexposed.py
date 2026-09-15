#!/usr/bin/env python3
"""Outcome-blind tests for the GSE274114 within-platform activation audit."""

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from scripts import audit_gse274114_within_platform_activation_unexposed as audit


ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "executions/model-data-080-21082249-gse274114"


class GSE274114WithinPlatformActivationAuditTests(unittest.TestCase):
    def test_pinned_fixture_controls_are_selectively_audited(self) -> None:
        observed = audit.audit_fixture_controls_without_protected_bytes(FIXTURE)
        self.assertEqual(observed["evaluator_only_files_declared_not_opened"], 2)
        self.assertEqual(observed["molecular_quant_files_declared_not_opened"], 39)

    def test_right_reference_topology_contracts_pass(self) -> None:
        observed = audit.validate_authorities(FIXTURE)
        self.assertEqual(observed["participants"], 39)
        self.assertEqual(observed["technical_runs"], 71)
        self.assertFalse(observed["participant_metadata_inputs_available"])

    def test_two_donor_safe_task_specs_prohibit_global_claims(self) -> None:
        observed = audit.validate_task_specs(
            ROOT / "config/datasets/gse274114_mash_hbv.toml",
            ROOT / "config/evaluation/gse274114_hiseq_healthy_vs_hbv_task.toml",
            ROOT / "config/evaluation/gse274114_novaseq_mash_vs_mash_hbv_task.toml",
            ROOT / "config/evaluation/gse274114_within_instrument_etiology_promotion_gate.json",
        )
        self.assertEqual(observed["hiseq_participants"], 20)
        self.assertEqual(observed["novaseq_participants"], 19)
        self.assertFalse(observed["global_four_class_performance_allowed"])
        self.assertFalse(observed["global_ood_accuracy_allowed"])
        self.assertFalse(observed["champion_eligible"])

    def test_auditor_never_uses_full_tree_verifier(self) -> None:
        source = Path(audit.__file__).read_text(encoding="utf-8")
        self.assertNotIn("verify_frozen_tree", source)
        self.assertNotIn("evaluator_only/source_labels.tsv\").read", source)
        self.assertNotIn("evaluator_only/participant_folds.tsv\").read", source)


if __name__ == "__main__":
    unittest.main()

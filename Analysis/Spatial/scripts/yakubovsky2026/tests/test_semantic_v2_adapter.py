#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "09_build_plan13_adapter_semantic_v2.py"
SPEC = importlib.util.spec_from_file_location("yakubovsky_semantic_v2_adapter", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def effect_fixture() -> dict[str, str]:
    return {
        "program_uid": "p1",
        "membership_sha256": "membership",
        "estimate": "0.1",
        "pvalue": "0.4",
        "padj": "0.8",
        "testable": "TRUE",
        "direction_observed": "positive",
        "robustness_pass": "FALSE",
        "evidence_state": "tested_negative",
    }


def native_fixture() -> dict[str, str]:
    return {
        "program_uid": "p1",
        "membership_sha256": "membership",
        "estimate": "0.1",
        "pvalue": "0.4",
        "qvalue": "0.8",
        "combined_z": "-1.25",
        "combined_direction": "negative",
        "median_donor_slope_direction": "positive",
        "combined_and_median_direction_agree": "FALSE",
        "cochran_q_descriptive": "12.5",
        "cochran_q_df": "2",
        "cochran_q_p_descriptive": "0.0019",
    }


class SemanticV2AdapterTests(unittest.TestCase):
    def test_nonsignificant_without_adequacy_rule_is_indeterminate(self):
        self.assertEqual(
            MODULE.semantic_evidence_state("tested_negative", False, ""),
            "indeterminate",
        )
        self.assertEqual(
            MODULE.semantic_evidence_state(
                "tested_negative", False, "prespecified_adequate_negative_rule_v1"
            ),
            "tested_negative",
        )

    def test_semantic_row_separates_descriptive_and_inferential_directions(self):
        row = MODULE.semantic_effect_row(
            effect_fixture(),
            native_fixture(),
            producer="Analysis/Spatial/scripts/yakubovsky2026/09_build_plan13_adapter_semantic_v2.py",
            producer_sha256="producer_hash",
            source_manifest_sha256="source_hash",
        )
        self.assertEqual(row["evidence_state"], "indeterminate")
        self.assertEqual(row["descriptive_effect_direction"], "positive")
        self.assertEqual(row["inferential_test_direction"], "negative")
        self.assertFalse(row["direction_agreement"])
        self.assertEqual(row["heterogeneity_statistic"], "12.5")
        self.assertEqual(row["heterogeneity_df"], "2")
        self.assertEqual(row["heterogeneity_pvalue"], "0.0019")
        self.assertEqual(row["matched_null_sd"], "")
        self.assertEqual(row["negative_call_rule_id"], "")
        self.assertEqual(row["producer_sha256"], "producer_hash")

    def test_native_direction_flag_drift_fails_closed(self):
        native = native_fixture()
        native["combined_and_median_direction_agree"] = "TRUE"
        with self.assertRaises(MODULE.ContractError):
            MODULE.semantic_effect_row(
                effect_fixture(),
                native,
                producer="producer",
                producer_sha256="hash",
                source_manifest_sha256="source",
            )

    def test_robust_and_terminal_states_cannot_carry_contradictory_flags(self):
        self.assertEqual(MODULE.semantic_evidence_state("robust", True), "robust")
        self.assertEqual(
            MODULE.semantic_evidence_state("untestable", False), "untestable"
        )
        with self.assertRaises(MODULE.ContractError):
            MODULE.semantic_evidence_state("robust", False)
        with self.assertRaises(MODULE.ContractError):
            MODULE.semantic_evidence_state("untestable", True)

    def test_gate_is_indeterminate_for_mixed_nonrobust_family(self):
        status, reason = MODULE._gate_from_states(
            ["indeterminate", "indeterminate"], True
        )
        self.assertEqual(status, "indeterminate")
        self.assertIn("without_prespecified_adequate_negative_rule", reason)


if __name__ == "__main__":
    unittest.main()

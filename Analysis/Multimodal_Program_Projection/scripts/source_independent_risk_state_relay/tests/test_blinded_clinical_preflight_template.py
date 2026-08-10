from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


SCRIPT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = SCRIPT_ROOT / "63_build_blinded_clinical_preflight_template.py"
SPEC = importlib.util.spec_from_file_location("clinical_preflight", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class BlindedClinicalPreflightTemplateTests(unittest.TestCase):
    def test_primary_rule_is_exact_and_not_substitutable(self) -> None:
        primary = MODULE.RESPONDER_ROWS[0]
        self.assertEqual(primary["endpoint_id"], "primary_histologic_improvement")
        self.assertIn("NAS_decrease_ge_2", primary["improver_definition"])
        self.assertIn("no_fibrosis_worsening", primary["improver_definition"])
        self.assertEqual(primary["substitution_allowed"], "false")

    def test_primary_state_effect_has_expected_negative_direction(self) -> None:
        primary = MODULE.ESTIMANDS[0]
        self.assertEqual(primary["analysis_unit"], "participant")
        self.assertEqual(primary["expected_direction"], "negative")
        self.assertIn("permuted_within_frozen_trial_and_intervention_strata", primary["model"])

    def test_meta_analysis_cannot_rescue_a_failed_cohort(self) -> None:
        meta = MODULE.ESTIMANDS[2]
        self.assertIn("both_individual_p_lt_0.05", meta["inference"])
        self.assertIn("cannot_rescue", meta["promotion_role"])

    def test_power_policy_prohibits_state_outcome_variance(self) -> None:
        policies = {parameter: value for parameter, value, _ in MODULE.POWER_POLICIES}
        self.assertEqual(policies["variance_source"], "external_or_response_blinded_only")
        self.assertEqual(policies["underpowered_action"], "fail_source_gate_or_later_manuscript")

    def test_all_orthogonal_assays_prohibit_substitution(self) -> None:
        self.assertEqual(
            {row["assay_id"] for row in MODULE.ORTHOGONAL_ROWS},
            {"tissue_proteomics", "liver_linked_secreted_proteomics", "physical_niche"},
        )
        self.assertTrue(all(row["substitution_allowed"] == "false" for row in MODULE.ORTHOGONAL_ROWS))


if __name__ == "__main__":
    unittest.main()


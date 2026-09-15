"""The seal must carry the honesty constraints, not merely mention them.

Three things are pinned. External support is set-level and the document has to
say so where a reader will meet it. The `samples: 164` note has to be present
and has to name 106, because that field in the frozen substrate artifact is the
exact shape of error this project has hit three times. And S2 has to be power
guarded, because it accepts a null and an unguarded accepted null is the
difference between tested_negative and indeterminate.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "scripts" / "freeze_two_axis_external_prespecification.py"
SUBSTRATE = (
    ROOT / "executions" / "model-data-953-21182514-two-axis-substrate"
    / "substrate" / "substrate.json"
)
ACTIVATION = (
    ROOT / "config" / "evaluation"
    / "two_axis_activity_fibrosis_external_activation.json"
)


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


freezer = _load(MODULE, "two_axis_external_freezer_tested")


def _payload() -> dict:
    return freezer.build(
        json.loads(SUBSTRATE.read_text(encoding="utf-8")),
        json.loads(ACTIVATION.read_text(encoding="utf-8")),
    )


class HonestyConstraintTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = _payload()

    def test_external_support_is_declared_set_level(self) -> None:
        block = self.payload["the_deliverable_and_its_honesty_constraint"]
        self.assertTrue(block["external_support_is_set_level_not_per_gene"])
        self.assertIn("will not reproduce", block["why"])
        self.assertTrue(self.payload["claim_boundary"]["external_support_is_set_level"])

    def test_the_samples_164_note_is_present_and_names_106(self) -> None:
        unit = self.payload["unit_of_inference"]["GSE193066"]
        note = unit["REQUIRED_NOTE_ON_THE_SUBSTRATE_ARTIFACT"]
        self.assertIn("164", note)
        self.assertIn("106", note)
        self.assertIn("NOT the analysis unit", note)
        self.assertEqual(unit["n"], 106)
        self.assertEqual(unit["rows_in_the_deposited_matrix"], 164)
        self.assertTrue(unit["the_analysis_must_assert_exactly_106_distinct_participants"])

    def test_the_suffix_trap_is_named_rather_than_worked_around(self) -> None:
        rule = self.payload["unit_of_inference"]["GSE193066"]["selection_rule"]
        self.assertIn("pairs zero of 58", rule)
        self.assertIn("!Sample_title", rule)

    def test_gse130970_people_are_not_assertable(self) -> None:
        unit = self.payload["unit_of_inference"]["GSE130970"]
        self.assertFalse(unit["distinct_people_assertable"])
        self.assertEqual(unit["unit"], "bulk_rna_sample")

    def test_the_fibrosis_claim_boundary_is_verbatim(self) -> None:
        self.assertEqual(
            self.payload["claim_boundary"][
                "the_fibrosis_axis_has_exactly_one_external_arm"],
            "The fibrosis axis has exactly one external arm, n=78 samples, "
            "with no assertable donor key.",
        )


class CriteriaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = _payload()

    def test_s2_is_power_guarded_against_accepting_a_null(self) -> None:
        s2 = self.payload["criteria"][
            "s2_the_assigned_set_is_not_elevated_on_the_other_axis"]
        self.assertIn("indeterminate rather than met", s2["threshold"])
        self.assertIn("tested_negative", s2["s2_accepts_a_null_so_it_is_power_guarded"])

    def test_both_criteria_exclude_the_both_cell(self) -> None:
        for name in freezer.FROZEN_THRESHOLDS:
            self.assertEqual(
                self.payload["criteria"][name]["evaluated_on"],
                ["activity_only", "fibrosis_only"],
            )

    def test_the_set_level_mde_is_operative_and_the_single_gene_one_is_not(self) -> None:
        power = self.payload["power"]
        self.assertFalse(power["set_level_mde_is_the_operative_reference"]["is_a_gate"])
        self.assertIn(
            "method-to-estimand mismatch",
            power["set_level_mde_is_the_operative_reference"][
                "why_it_is_computed_rather_than_borrowed"],
        )
        self.assertIn(
            "does not make",
            power["single_correlation_mde_is_context_only"][
                "what_they_are_the_reference_for"],
        )

    def test_the_ceiling_caveat_is_recorded(self) -> None:
        self.assertIn("not that an effect exists", self.payload["power"]["ceilings"]["caveat"])


class ArmScopeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = _payload()

    def test_gse193066_serves_activity_only(self) -> None:
        self.assertEqual(freezer.ARM_SCOPE["GSE193066"]["axes_served"], ["activity"])
        self.assertEqual(
            freezer.ARM_SCOPE["GSE130970"]["axes_served"], ["activity", "fibrosis"])

    def test_the_omitted_level_is_the_ground_and_the_others_are_caveats(self) -> None:
        block = self.payload["arm_scope"]["gse193066_is_restricted_to_the_activity_axis"]
        self.assertIn("omits stage 4", block["ground_that_carries_the_decision"])
        self.assertIn("categorical", block["ground_that_carries_the_decision"])
        caveats = " ".join(block["supporting_coverage_caveats_not_independent_grounds"])
        self.assertIn("it does not bias S2", caveats)
        self.assertIn("weakest of the three", caveats)

    def test_the_fibrosis_arm_is_computed_but_uncitable_either_way(self) -> None:
        diagnostic = self.payload["diagnostics_never_gates"]["gse193066_fibrosis_arm"]
        self.assertFalse(diagnostic["is_a_gate"])
        self.assertTrue(diagnostic["enters_no_outcome_cell"])
        self.assertEqual(len(diagnostic["not_citable_as"]), 2)


class SubstrateGuardTests(unittest.TestCase):
    def test_a_substrate_that_does_not_reconstruct_the_marginals_is_refused(self) -> None:
        payload = json.loads(SUBSTRATE.read_text(encoding="utf-8"))
        payload["assignment"]["cell_sizes"]["activity_only"] += 1
        with self.assertRaises(freezer.PrespecificationError):
            freezer.build(payload, json.loads(ACTIVATION.read_text(encoding="utf-8")))

    def test_a_substrate_without_the_164_column_matrix_is_refused(self) -> None:
        payload = json.loads(SUBSTRATE.read_text(encoding="utf-8"))
        payload["external_arms"]["GSE193066"]["samples"] = 106
        with self.assertRaises(freezer.PrespecificationError):
            freezer.build(payload, json.loads(ACTIVATION.read_text(encoding="utf-8")))


if __name__ == "__main__":
    unittest.main()

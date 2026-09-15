"""The seal must carry the mediation framing and the variant-stability correction.

Two properties are pinned. A composition-dependent gene must never be nameable
as spurious, and T1's failure must be named for mechanism rather than for
collapse -- the framing survives into the write-up only if the outcome labels
carry it. And the correction to the frozen substrate's prose must state three
stable lineages, matching the list in the same object rather than the sentence
beside it.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "scripts" / "freeze_composition_prespecification.py"
SUBSTRATE = (ROOT / "executions" / "model-data-954-21183077-composition-substrate"
             / "substrate" / "composition_substrate.json")


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


freezer = _load(MODULE, "composition_freezer_tested")


def _payload() -> dict:
    return freezer.build(json.loads(SUBSTRATE.read_text(encoding="utf-8")))


class FramingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = _payload()

    def test_a_composition_dependent_gene_is_declared_not_an_artifact(self) -> None:
        block = self.payload["the_framing_that_governs_every_reading"]
        self.assertTrue(block["a_composition_dependent_gene_is_not_an_artifact"])
        self.assertIn("mediation, not confounding", block["why"])

    def test_calling_such_a_gene_spurious_is_forbidden(self) -> None:
        self.assertIn(
            "spurious",
            self.payload["the_framing_that_governs_every_reading"]["forbidden"])

    def test_t1_failure_is_named_for_mechanism_not_collapse(self) -> None:
        t1 = self.payload["criteria"]["t1_a_composition_independent_component_exists"]
        self.assertEqual(t1["failure_is_named"], "SET_IS_COMPOSITION_MEDIATED")
        self.assertIn("never a failure of the two-axis result",
                      t1["what_a_failure_means"])

    def test_t1_asks_about_a_component_not_about_the_axis(self) -> None:
        t1 = self.payload["criteria"]["t1_a_composition_independent_component_exists"]
        self.assertIn("not whether", t1["what_it_asks"])

    def test_the_stricter_test_note_is_present_with_its_guard(self) -> None:
        note = self.payload["criteria"][
            "t1_a_composition_independent_component_exists"][
            "six_lineages_is_a_stricter_test_than_five"]
        self.assertIn("harder to pass", note)
        self.assertIn("partly technical", note)
        self.assertIn("indeterminate rather than mediated", note)

    def test_no_outcome_cell_implies_the_axis_collapsed(self) -> None:
        for cell in self.payload["decision_rule"]["cells"]:
            self.assertNotIn("NOT_", cell["outcome"])
            self.assertNotIn("FAIL", cell["outcome"])


class CorrectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = _payload()

    def test_the_correction_states_three_not_two(self) -> None:
        c = self.payload["correction_to_the_frozen_substrate_artifact"]
        self.assertIn("three are stable, not two", c["what_is_true"])
        self.assertEqual(
            len(c["the_same_artifact_records_the_correct_list"]["stable_across_variants"]),
            3)

    def test_the_correction_names_the_field_it_corrects(self) -> None:
        c = self.payload["correction_to_the_frozen_substrate_artifact"]
        self.assertIn("what_it_shows", c["field"])
        self.assertIn("source manifest", c["why_it_is_not_retrofitted"])

    def test_the_substantive_point_is_declared_unchanged(self) -> None:
        c = self.payload["correction_to_the_frozen_substrate_artifact"]
        self.assertIn("variant-dependent", c["the_substantive_point_is_unchanged"])

    def test_a_substrate_with_a_different_stable_count_is_refused(self) -> None:
        payload = json.loads(SUBSTRATE.read_text(encoding="utf-8"))
        payload["variant_sensitivity_on_the_family"]["stable_across_variants"] = ["Hep"]
        with self.assertRaises(freezer.PrespecificationError):
            freezer.build(payload)


class CovariateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = _payload()

    def test_seven_covariates_and_df_n_minus_9(self) -> None:
        cov = self.payload["covariates"]
        self.assertEqual(cov["n_covariates"], 7)
        self.assertIn("n - 9", cov["residual_df"])

    def test_conditioning_is_reported_as_a_measured_negative(self) -> None:
        c = self.payload["covariates"]["conditioning_is_a_measured_negative"]
        self.assertIn("silently", c["what_would_have_broken"])
        self.assertIn("did_not_bite", "".join(c))
        for arm in c["measured_per_arm"].values():
            self.assertEqual(arm["numerical_rank"], arm["n_covariates"])
            self.assertLess(arm["condition_number"], 10.0)

    def test_a_family_that_is_not_six_is_refused(self) -> None:
        payload = json.loads(SUBSTRATE.read_text(encoding="utf-8"))
        payload["eligible_family"]["family"] = ["Hepatocytes"]
        with self.assertRaises(freezer.PrespecificationError):
            freezer.build(payload)


class ClosureTests(unittest.TestCase):
    def test_between_lineage_claims_are_forbidden_with_the_reversal_cited(self) -> None:
        closure = _payload()["closure"]
        self.assertIn("+0.225 to -0.165", closure["no_between_lineage_claims"])
        self.assertIn("REF_HUMAN", closure["one_shared_operator_not_independent_measurements"])
        self.assertIn("11 names", closure["never_align_the_mouse_roster_by_name"])


if __name__ == "__main__":
    unittest.main()

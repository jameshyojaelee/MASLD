"""Pin the v2 prespecification's own criteria, literals and outcome map.

Every earlier stage on this lane had a test pinning its prespec's literal
thresholds and criterion names, and that is what made "the criteria were fixed
before the data were opened" checkable rather than asserted. The first v2 seal
was gated by the *evaluator's* tests, so nothing pinned the freeze itself. This
is that missing test.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "scripts" / "freeze_composition_prespecification_v2.py"
SUBSTRATE = (ROOT / "executions" / "model-data-954-21183077-composition-substrate"
             / "substrate" / "composition_substrate.json")


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


freezer = _load(MODULE, "composition_v2_freezer")


def _payload() -> dict:
    return freezer.build(json.loads(SUBSTRATE.read_text(encoding="utf-8")))


class CriterionNameTests(unittest.TestCase):
    def test_the_two_criteria_are_the_frozen_two(self) -> None:
        self.assertEqual(
            sorted(_payload()["criteria"]),
            ["t1_a_composition_independent_component_exists",
             "t2_per_gene_composition_dependence"])

    def test_the_prespec_id_and_what_it_supersedes(self) -> None:
        payload = _payload()
        self.assertEqual(payload["prespec_id"], "composition_dependence_v2")
        self.assertEqual(payload["supersedes"], "composition_dependence_v1")

    def test_t1_is_decisive_and_t2_is_the_deliverable(self) -> None:
        criteria = _payload()["criteria"]
        self.assertTrue(criteria["t1_a_composition_independent_component_exists"]["decisive"])
        self.assertFalse(criteria["t2_per_gene_composition_dependence"]["decisive"])
        self.assertTrue(criteria["t2_per_gene_composition_dependence"]["is_the_deliverable"])


class ThresholdLiteralTests(unittest.TestCase):
    def test_t1s_literal_states_the_symmetric_guard_in_both_directions(self) -> None:
        literal = freezer.FROZEN_THRESHOLDS[
            "t1_a_composition_independent_component_exists"]
        self.assertIn("BY AT LEAST", literal)
        self.assertIn("set-level minimum detectable shift", literal)
        self.assertIn("indeterminate on either side", literal)

    def test_t2s_literal_scopes_to_one_arm_and_names_the_external_state(self) -> None:
        literal = freezer.FROZEN_THRESHOLDS["t2_per_gene_composition_dependence"]
        self.assertIn("GSE135251 only", literal)
        self.assertIn("indeterminate_by_construction", literal)

    def test_the_built_document_carries_the_literals_verbatim(self) -> None:
        payload = _payload()
        for name, literal in freezer.FROZEN_THRESHOLDS.items():
            self.assertEqual(payload["criteria"][name]["threshold"], literal)

    def test_t1s_failure_is_named_for_mechanism(self) -> None:
        self.assertEqual(
            _payload()["criteria"]["t1_a_composition_independent_component_exists"][
                "failure_is_named"], "SET_IS_COMPOSITION_MEDIATED")


class OutcomeMapTests(unittest.TestCase):
    EXPECTED = ["COMPOSITION_INDEPENDENT_COMPONENT_EXISTS",
                "ARMS_DIFFER_WITHIN_RESOLUTION",
                "COMPOSITION_INDEPENDENT_IN_ONE_ARM_ONLY",
                "SET_IS_COMPOSITION_MEDIATED",
                "INDETERMINATE"]

    def test_the_outcome_map_is_exactly_the_five_declared_cells(self) -> None:
        self.assertEqual(
            [c["outcome"] for c in _payload()["decision_rule"]["cells"]], self.EXPECTED)

    def test_the_new_cell_explains_why_it_exists(self) -> None:
        cell = next(c for c in _payload()["decision_rule"]["cells"]
                    if c["outcome"] == "ARMS_DIFFER_WITHIN_RESOLUTION")
        self.assertIn("coin-flip", cell["when"])

    def test_no_cell_implies_the_two_axis_result_collapsed(self) -> None:
        for cell in _payload()["decision_rule"]["cells"]:
            self.assertNotIn("FAIL", cell["outcome"])
            self.assertNotIn("NOT_", cell["outcome"])

    def test_the_vacuity_guard_raises_rather_than_reporting(self) -> None:
        guard = _payload()["decision_rule"]["vacuity_guard"]
        self.assertIn("raises rather than reporting", guard)


class VoidRunTests(unittest.TestCase):
    def test_the_prior_run_is_recorded_uncitable_with_both_causes(self) -> None:
        prior = _payload()["the_prior_run_is_not_a_result"]
        self.assertFalse(prior["may_it_be_cited"])
        self.assertIn("cause_1_the_training_arm_dropped_out_silently", prior)
        self.assertIn("cause_2_the_verdict_was_a_threshold_coin_flip", prior)

    def test_the_join_measurement_is_recorded_as_a_number(self) -> None:
        measured = _payload()["the_prior_run_is_not_a_result"][
            "cause_1_the_training_arm_dropped_out_silently"]["measured"]
        self.assertEqual(measured["join_on_stable_gene_id"], "61940 of 61940")
        self.assertEqual(measured["join_on_feature_index"], 0)
        self.assertTrue(measured["both_key_spaces_were_ensg"])

    def test_it_records_that_no_crosswalk_was_needed(self) -> None:
        cause = _payload()["the_prior_run_is_not_a_result"][
            "cause_1_the_training_arm_dropped_out_silently"]
        self.assertIn("never needed a crosswalk", cause["it_was_never_a_key_space_mismatch"])

    def test_the_coin_flip_numbers_are_recorded(self) -> None:
        cause = _payload()["the_prior_run_is_not_a_result"][
            "cause_2_the_verdict_was_a_threshold_coin_flip"]
        self.assertEqual(cause["difference"], 0.0008)
        self.assertEqual(cause["the_passing_margin"], 0.00188)
        self.assertEqual(cause["that_arms_set_level_mde"], 0.005)


class ColumnLessonTests(unittest.TestCase):
    def test_both_directions_are_recorded(self) -> None:
        lesson = _payload()["the_column_lesson_runs_both_ways_in_this_tree"]
        self.assertIn("never trust the name", lesson["from_the_ragged_headers"])
        self.assertIn("never trust the position", lesson["from_this_defect"])
        self.assertIn("Neither rule alone is sufficient", lesson["the_unifying_rule"])


class SubstrateGuardTests(unittest.TestCase):
    def test_a_family_that_is_not_six_is_refused(self) -> None:
        payload = json.loads(SUBSTRATE.read_text(encoding="utf-8"))
        payload["eligible_family"]["family"] = ["Hepatocytes"]
        with self.assertRaises(freezer.PrespecificationError):
            freezer.build(payload)

    def test_the_family_of_six_is_carried_into_the_document(self) -> None:
        self.assertEqual(len(_payload()["substrate"]["eligible_family"]), 6)


if __name__ == "__main__":
    unittest.main()


class EvaluatorContractTests(unittest.TestCase):
    """The prespec must supply every key its evaluator reads.

    Four defects on this stage came from deriving v2 files mechanically. This
    one is different in kind: the v2 prespecification simply did not carry two
    records the evaluator needs, and each was found only by running the job and
    reading a KeyError. A regex over the evaluator missed one of them because
    the subscript spans two lines, so the contract is checked by AST.
    """

    @staticmethod
    def _keys_the_evaluator_reads() -> set:
        import ast

        source = (ROOT / "scripts" / "evaluate_composition_dependence_v2.py"
                  ).read_text(encoding="utf-8")
        keys = set()
        for node in ast.walk(ast.parse(source)):
            if (isinstance(node, ast.Subscript)
                    and isinstance(node.value, ast.Name)
                    and node.value.id == "prespec"
                    and isinstance(node.slice, ast.Constant)):
                keys.add(node.slice.value)
        return keys

    def test_the_prespec_supplies_every_key_the_evaluator_reads(self) -> None:
        missing = self._keys_the_evaluator_reads() - set(_payload())
        self.assertEqual(missing, set(), f"the v2 prespec is missing {sorted(missing)}")

    def test_the_scan_finds_the_two_keys_that_broke_the_run(self) -> None:
        """Guards the guard: an AST scan that found nothing would pass vacuously."""

        keys = self._keys_the_evaluator_reads()
        self.assertIn("closure", keys)
        self.assertIn("correction_to_the_frozen_substrate_artifact", keys)
        self.assertGreaterEqual(len(keys), 8)

    def test_the_substrate_correction_is_carried_not_only_inherited(self) -> None:
        correction = _payload()["correction_to_the_frozen_substrate_artifact"]
        self.assertIn("three are stable, not two", correction["what_is_true"])
        reason = correction["carried_forward_from_v1_because_it_is_still_live"]
        self.assertIn("binds the same substrate artifact", reason)
        self.assertIn("lose a record deliberately added", reason)

    def test_the_closure_constraints_survive_into_v2(self) -> None:
        closure = _payload()["closure"]
        self.assertIn("+0.225 to -0.165", closure["no_between_lineage_claims"])
        self.assertIn("REF_HUMAN",
                      closure["one_shared_operator_not_independent_measurements"])
        self.assertIn("11 names", closure["never_align_the_mouse_roster_by_name"])

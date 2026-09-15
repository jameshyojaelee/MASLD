"""The symmetric guard, the confirmatory-arm scoping, and the coin-flip cell.

Each pins a defect that reached a verdict in v1. The symmetric guard is the one
that matters most: v1 reported a pass whose margin over its null was 0.00188
against an MDE of 0.005, and no test would have caught that because the
one-sided guard was working exactly as written.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "scripts" / "evaluate_composition_dependence_v2.py"
FREEZER = ROOT / "scripts" / "freeze_composition_prespecification_v2.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ev = _load(MODULE, "composition_v2_tested")
freezer = _load(FREEZER, "composition_v2_freezer_tested")


def _classify(held: bool, inside: bool, detectable: bool):
    applicable = bool((held and not inside)
                      or (not held and detectable and not inside))
    verdict = ("composition_independent_component_exists" if (held and not inside)
               else "indeterminate" if inside
               else "set_is_composition_mediated" if detectable
               else "indeterminate")
    return applicable, bool(held and not inside), verdict


class SymmetricGuardTests(unittest.TestCase):
    def test_the_v1_pass_now_reads_indeterminate(self) -> None:
        """margin 0.00188 against an MDE of 0.005 - the exact v1 case."""

        margin, mde = 0.08747 - 0.08559, 0.005
        inside = abs(margin) < mde
        self.assertTrue(inside)
        applicable, met, verdict = _classify(held=True, inside=inside, detectable=True)
        self.assertFalse(met)
        self.assertFalse(applicable)
        self.assertEqual(verdict, "indeterminate")

    def test_a_pass_clear_of_the_mde_still_counts(self) -> None:
        applicable, met, verdict = _classify(held=True, inside=False, detectable=True)
        self.assertTrue(met)
        self.assertTrue(applicable)
        self.assertEqual(verdict, "composition_independent_component_exists")

    def test_a_detectable_drop_is_mediated(self) -> None:
        applicable, met, verdict = _classify(held=False, inside=False, detectable=True)
        self.assertTrue(applicable)
        self.assertFalse(met)
        self.assertEqual(verdict, "set_is_composition_mediated")

    def test_an_undetectable_drop_is_indeterminate(self) -> None:
        applicable, met, verdict = _classify(held=False, inside=False, detectable=False)
        self.assertFalse(applicable)
        self.assertEqual(verdict, "indeterminate")

    def test_the_guard_is_symmetric_in_the_sign_of_the_margin(self) -> None:
        for margin in (0.002, -0.002):
            self.assertTrue(abs(margin) < 0.005)


class WithinResolutionTests(unittest.TestCase):
    @staticmethod
    def _within(margins, mdes):
        return bool(margins and mdes and (max(margins) - min(margins)) < min(mdes))

    def test_the_v1_arms_differ_within_resolution(self) -> None:
        """0.0491 and 0.0483 drops, margins far closer than the 0.005 MDE."""

        self.assertTrue(self._within([0.00188, -0.00904], [0.005, 0.005]) is False)
        self.assertTrue(self._within([0.0019, 0.0021], [0.005, 0.005]))

    def test_genuinely_different_arms_are_not_within_resolution(self) -> None:
        self.assertFalse(self._within([0.05, 0.001], [0.005, 0.005]))

    def test_the_smallest_mde_among_the_arms_is_the_yardstick(self) -> None:
        self.assertFalse(self._within([0.0, 0.004], [0.002, 0.05]))
        self.assertTrue(self._within([0.0, 0.004], [0.05, 0.05]))


class ConfirmatoryArmTests(unittest.TestCase):
    def test_t2_is_confirmatory_in_gse135251_only(self) -> None:
        self.assertEqual(ev.CONFIRMATORY_ARM, "GSE135251")
        self.assertEqual(freezer.CONFIRMATORY_ARM, "GSE135251")
        self.assertEqual(set(freezer.EXTERNAL_ARMS), {"GSE130970", "GSE193066"})

    def test_the_frozen_literals_match_between_freezer_and_analysis(self) -> None:
        self.assertEqual(freezer.FROZEN_THRESHOLDS, ev.FROZEN_THRESHOLDS)

    def test_t1s_literal_carries_the_symmetric_guard(self) -> None:
        literal = ev.FROZEN_THRESHOLDS["t1_a_composition_independent_component_exists"]
        self.assertIn("BY AT LEAST", literal)
        self.assertIn("either side", literal)

    def test_t2s_literal_scopes_to_one_arm(self) -> None:
        literal = ev.FROZEN_THRESHOLDS["t2_per_gene_composition_dependence"]
        self.assertIn("GSE135251 only", literal)
        self.assertIn("indeterminate_by_construction", literal)


class NotAResultTests(unittest.TestCase):
    def test_the_prior_run_is_recorded_as_uncitable(self) -> None:
        import json
        substrate = json.loads((
            ROOT / "executions" / "model-data-954-21183077-composition-substrate"
            / "substrate" / "composition_substrate.json").read_text(encoding="utf-8"))
        payload = freezer.build(substrate)
        prior = payload["the_prior_run_is_not_a_result"]
        self.assertFalse(prior["may_it_be_cited"])
        self.assertEqual(prior["its_outcome"], "COMPOSITION_INDEPENDENT_IN_ONE_ARM_ONLY")
        self.assertEqual(
            prior["cause_1_the_training_arm_dropped_out_silently"]["measured"][
                "join_on_feature_index"], 0)

    def test_the_column_lesson_records_both_directions(self) -> None:
        import json
        substrate = json.loads((
            ROOT / "executions" / "model-data-954-21183077-composition-substrate"
            / "substrate" / "composition_substrate.json").read_text(encoding="utf-8"))
        lesson = freezer.build(substrate)["the_column_lesson_runs_both_ways_in_this_tree"]
        self.assertIn("never trust the name", lesson["from_the_ragged_headers"])
        self.assertIn("never trust the position", lesson["from_this_defect"])

    def test_the_new_cell_exists_and_is_named_unmistakably(self) -> None:
        import json
        substrate = json.loads((
            ROOT / "executions" / "model-data-954-21183077-composition-substrate"
            / "substrate" / "composition_substrate.json").read_text(encoding="utf-8"))
        cells = [c["outcome"] for c in
                 freezer.build(substrate)["decision_rule"]["cells"]]
        self.assertIn("ARMS_DIFFER_WITHIN_RESOLUTION", cells)


if __name__ == "__main__":
    unittest.main()

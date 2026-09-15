"""Arm B must run on Arm A's instrument, and E2 must never override E1.

The point of Stage 0c is a comparison across substrates, which is only
meaningful if the statistic is identical. So the first thing pinned is that this
script defines none of the machinery itself and imports all of it. The second is
E2's classification rule, which is the criterion Stage 0b showed decides the
interpretation: below the measured floor an empty count is underpowered, at or
above it the same empty count is evidence. The boundary case is pinned
explicitly because "reaches or clears" is the frozen literal.
"""

from __future__ import annotations

import ast
import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "scripts" / "evaluate_gse135251_axis_count.py"
ARM_A = ROOT / "scripts" / "evaluate_gse267145_axis_count.py"
FREEZER = ROOT / "scripts" / "freeze_arm_b_axis_count_prespecification.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


arm_b = _load(MODULE, "arm_b_analysis_tested")
freezer = _load(FREEZER, "arm_b_freezer_for_analysis_tested")


class SameInstrumentTests(unittest.TestCase):
    def test_thresholds_match_the_freezer_literals(self) -> None:
        self.assertEqual(freezer.FROZEN_THRESHOLDS, arm_b.FROZEN_THRESHOLDS)

    def test_the_statistics_are_imported_not_redefined(self) -> None:
        defined = {
            node.name
            for node in ast.walk(ast.parse(MODULE.read_text(encoding="utf-8")))
            if isinstance(node, ast.FunctionDef)
        }
        for owned_by_arm_a in (
            "partial_family", "unit_ranks", "residualize", "unit_columns",
            "average_ranks", "centred_ranks", "bh_count", "bh_critical_abs_r",
            "evaluate_gate", "direction_condition", "marginal_family",
        ):
            self.assertNotIn(owned_by_arm_a, defined)

    def test_arm_a_exposes_everything_arm_b_imports(self) -> None:
        analysis = arm_b._import(ARM_A, "arm_a_probe")
        for name in (
            "partial_family", "unit_ranks", "residualize", "unit_columns",
            "average_ranks", "centred_ranks", "bh_count_from_abs_r",
            "bh_critical_abs_r", "evaluate_gate", "direction_condition",
            "marginal_family", "tie_profile", "validate_rank_agreement",
            "read_table",
        ):
            self.assertTrue(hasattr(analysis, name), f"Arm A lacks {name}")

    def test_the_two_families_are_the_frozen_two(self) -> None:
        self.assertEqual(
            arm_b.FAMILIES,
            (("nas_score", "fibrosis_stage"), ("fibrosis_stage", "nas_score")),
        )


class EmptyDirectionClassificationTests(unittest.TestCase):
    @staticmethod
    def _result(count: int, observed_max: float, floor: float) -> dict:
        return {
            "exposure": "fibrosis_stage",
            "adjusted_for": "nas_score",
            "genes_bh_below_0_05": count,
            "observed_max_abs_partial_r": observed_max,
            "familywise_floor": {"detectable_abs_partial_r_p95": floor},
        }

    def test_below_the_floor_is_underpowered(self) -> None:
        entry = arm_b.classify_empty_direction(self._result(0, 0.3767, 0.4502))
        self.assertTrue(entry["applicable"])
        self.assertEqual(entry["classification"], "underpowered")
        self.assertIn("indeterminate", entry["reading"])

    def test_above_the_floor_is_not_independent(self) -> None:
        entry = arm_b.classify_empty_direction(self._result(0, 0.52, 0.4502))
        self.assertEqual(entry["classification"], "not_independent")
        self.assertIn("evidence", entry["reading"])

    def test_exactly_at_the_floor_is_not_independent(self) -> None:
        """The frozen literal says 'reaches or clears', so equality is not underpowered."""

        entry = arm_b.classify_empty_direction(self._result(0, 0.4502, 0.4502))
        self.assertEqual(entry["classification"], "not_independent")

    def test_a_non_empty_direction_is_inapplicable(self) -> None:
        entry = arm_b.classify_empty_direction(self._result(805, 0.6516, 0.4592))
        self.assertFalse(entry["applicable"])
        self.assertIsNone(entry["classification"])
        self.assertIn("not empty", entry["why_not_applicable"])

    def test_the_arm_a_numbers_reproduce_their_recorded_reading(self) -> None:
        """Stage 0b's two empty fibrosis directions must both read underpowered."""

        for observed, floor in ((0.3767, 0.4502), (0.4026, 0.4487), (0.3959, 0.4540)):
            self.assertEqual(
                arm_b.classify_empty_direction(self._result(0, observed, floor))[
                    "classification"],
                "underpowered",
            )


class JackknifeTests(unittest.TestCase):
    def test_the_jackknife_runs_one_fit_per_participant(self) -> None:
        analysis = arm_b._import(ARM_A, "arm_a_for_jackknife")
        generator = np.random.default_rng(2)
        n, g = 25, 120
        axes = {
            "nas_score": generator.integers(0, 9, size=n).astype(float),
            "fibrosis_stage": generator.integers(0, 5, size=n).astype(float),
        }
        summary = arm_b.jackknife_counts(
            analysis, generator.normal(size=(n, g)), axes, arm_b.FAMILIES
        )
        self.assertEqual(summary["method"], "leave_one_participant_out")
        self.assertEqual(summary["n_leave_one_out_fits"], n)
        for entry in summary["per_family"].values():
            self.assertEqual(entry["n_fits"], n)

    def test_the_reason_the_bootstrap_was_dropped_is_recorded(self) -> None:
        analysis = arm_b._import(ARM_A, "arm_a_for_jackknife_reason")
        generator = np.random.default_rng(3)
        summary = arm_b.jackknife_counts(
            analysis,
            generator.normal(size=(20, 60)),
            {
                "nas_score": generator.integers(0, 9, size=20).astype(float),
                "fibrosis_stage": generator.integers(0, 5, size=20).astype(float),
            },
            arm_b.FAMILIES,
        )
        self.assertIn("duplicat", summary["why_not_a_bootstrap"])


if __name__ == "__main__":
    unittest.main()

"""Lock the power machinery, above all that it never reads a realised effect."""

from __future__ import annotations

import ast
import inspect
from pathlib import Path
import unittest

from masld_bench.evaluators import detectable_effect as power
from masld_bench.evaluators.detectable_effect import (
    BH_SINGLE_TRUE_ALPHA,
    CONFIRMATORY_FAMILY_SIZE,
    NOMINAL_ALPHA,
    auprc_critical_values,
    auprc_power_curve,
    minimum_detectable_effect,
    required_sample_size_binary,
    resample_outcome,
    spearman_critical_values,
    spearman_power_curve,
)


class NonCircularityTests(unittest.TestCase):
    """The whole deliverable is void if an observed result leaks in here."""

    def test_module_cannot_perform_file_io_at_all(self) -> None:
        """The strongest available guarantee: it cannot read a result.

        A module with no file access, no network access, and no import that
        could reach either, cannot have consulted an observed metric however
        carelessly it were later edited.
        """

        tree = ast.parse(Path(power.__file__).read_text())
        called: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                target = node.func
                if isinstance(target, ast.Name):
                    called.add(target.id)
                elif isinstance(target, ast.Attribute):
                    called.add(target.attr)
        for forbidden in (
            "open", "read_text", "read_bytes", "load", "loads",
            "urlopen", "run", "check_output",
        ):
            self.assertNotIn(forbidden, called)

    def test_module_names_no_result_artifact(self) -> None:
        source = Path(power.__file__).read_text()
        for forbidden in (
            "predictions.tsv",
            "confirmatory_and_comparator_results",
            "standardized_metrics",
            "phenotype-scores",
            "21107640",
        ):
            self.assertNotIn(forbidden, source)

    def test_module_imports_nothing_that_could_carry_a_result(self) -> None:
        tree = ast.parse(Path(power.__file__).read_text())
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module)
            elif isinstance(node, ast.Import):
                imported.update(alias.name for alias in node.names)
        # ast records a relative "from .metrics import" as module "metrics".
        self.assertEqual(
            imported,
            {"__future__", "dataclasses", "math", "metrics", "random", "typing"},
        )

    def test_every_public_entry_point_takes_the_design_not_a_result(self) -> None:
        for name in (
            "spearman_power_curve",
            "auprc_power_curve",
            "required_sample_size",
            "required_sample_size_binary",
        ):
            signature = inspect.signature(getattr(power, name))
            for parameter in signature.parameters:
                self.assertNotIn("observed", parameter)
                self.assertNotIn("realised", parameter)
                self.assertNotIn("estimate", parameter)


class DecisionRuleTests(unittest.TestCase):
    def test_bh_single_true_reduces_to_bonferroni(self) -> None:
        self.assertEqual(CONFIRMATORY_FAMILY_SIZE, 12)
        self.assertAlmostEqual(BH_SINGLE_TRUE_ALPHA, NOMINAL_ALPHA / 12)

    def test_the_stricter_rule_has_the_higher_critical_value(self) -> None:
        critical = auprc_critical_values(15, 22, n_draws=8_000)
        self.assertGreater(critical["bh_single_true"], critical["nominal"])
        rank = spearman_critical_values(list(range(38)), n_draws=8_000)
        self.assertGreater(rank["bh_single_true"], rank["nominal"])

    def test_critical_values_sit_above_prevalence(self) -> None:
        critical = auprc_critical_values(15, 22, n_draws=8_000)
        self.assertGreater(critical["nominal"], 15 / 37)


class CurveShapeTests(unittest.TestCase):
    def test_power_is_monotone_in_the_planted_effect(self) -> None:
        points = auprc_power_curve(
            15, 22, [0.0, 0.5, 1.0, 1.5, 2.0], replicates=400, null_draws=6_000
        )
        nominal = [p.power_nominal for p in points]
        self.assertEqual(nominal, sorted(nominal))
        self.assertLess(nominal[0], 0.15)
        self.assertGreater(nominal[-1], nominal[0])

    def test_a_null_effect_recovers_the_stated_size(self) -> None:
        points = auprc_power_curve(
            15, 22, [0.0], replicates=4_000, null_draws=20_000
        )
        self.assertAlmostEqual(points[0].power_nominal, NOMINAL_ALPHA, delta=0.02)

    def test_rank_power_is_monotone_and_null_calibrated(self) -> None:
        outcome = [0.0] * 10 + [1.0] * 6 + [5.0] * 8 + [40.0] * 14
        points = spearman_power_curve(
            outcome, [0.0, 0.3, 0.6, 0.9], replicates=1_500, null_draws=12_000
        )
        nominal = [p.power_nominal for p in points]
        self.assertEqual(nominal, sorted(nominal))
        self.assertAlmostEqual(nominal[0], NOMINAL_ALPHA, delta=0.03)

    def test_bh_power_never_exceeds_nominal_power(self) -> None:
        points = auprc_power_curve(
            15, 22, [0.5, 1.0, 1.5], replicates=600, null_draws=6_000
        )
        for point in points:
            self.assertLessEqual(point.power_bh_single_true, point.power_nominal)

    def test_attenuation_lowers_power(self) -> None:
        strong = auprc_power_curve(
            15, 22, [1.5], replicates=800, null_draws=6_000, attenuation=1.0
        )[0]
        weak = auprc_power_curve(
            15, 22, [1.5], replicates=800, null_draws=6_000, attenuation=0.5
        )[0]
        self.assertLess(weak.power_nominal, strong.power_nominal)


class InversionTests(unittest.TestCase):
    def test_minimum_detectable_effect_interpolates(self) -> None:
        points = auprc_power_curve(
            15, 22, [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
            replicates=800, null_draws=8_000,
        )
        mde = minimum_detectable_effect(points, rule="nominal")
        self.assertIsNotNone(mde)
        assert mde is not None
        self.assertGreater(mde["latent_parameter"], 0.0)

    def test_power_rises_with_cohort_size(self) -> None:
        rows = required_sample_size_binary(
            15 / 37, [37, 150], 0.8, replicates=400, null_draws=6_000
        )
        self.assertLess(rows[0]["power_nominal"], rows[1]["power_nominal"])

    def test_resampling_preserves_the_outcome_shape(self) -> None:
        outcome = [0.0] * 10 + [85.0] * 2 + [5.0] * 6
        extended = resample_outcome(outcome, 200)
        self.assertEqual(len(extended), 200)
        self.assertTrue(set(extended).issubset(set(outcome)))


if __name__ == "__main__":
    unittest.main()


class VerdictTests(unittest.TestCase):
    """A verdict that can only say 'underpowered' is not a verdict."""

    def test_the_adequate_branch_is_reachable(self) -> None:
        self.assertEqual(
            power.power_verdict([39, 45], cohort_n=39),
            "adequately_powered_at_the_realised_cohort_size",
        )

    def test_a_moderate_shortfall_is_distinguished(self) -> None:
        self.assertEqual(
            power.power_verdict([80, 150], cohort_n=39),
            "underpowered_by_roughly_a_factor_of_two_to_five",
        )

    def test_an_order_of_magnitude_shortfall_is_named(self) -> None:
        self.assertEqual(
            power.power_verdict([100, 400], cohort_n=39),
            "underpowered_by_an_order_of_magnitude",
        )

    def test_unreached_effects_are_not_silently_treated_as_adequate(self) -> None:
        self.assertEqual(
            power.power_verdict([None, None], cohort_n=39),
            "no_reference_effect_reaches_target_power_within_the_tested_range",
        )

    def test_the_worst_reference_effect_governs(self) -> None:
        self.assertEqual(
            power.power_verdict([39, None, 600], cohort_n=39),
            "underpowered_by_an_order_of_magnitude",
        )

"""The axis-count analysis must agree with the package and with its own freezer.

Three things are pinned here, chosen because each is a place where this
analysis could be silently wrong and still produce a plausible number.

*The fast BH path.*  The count null recomputes Benjamini-Hochberg tens of
thousands of times, so the analysis reads the count straight off the
descending-sorted ``|partial r|`` against a precomputed critical vector instead
of evaluating an inverse-t inside the loop. That shortcut is only valid because
the two-sided t p-value is strictly decreasing in ``|r|`` at fixed df. It is
pinned against ``masld_bench.evaluators.stats.benjamini_hochberg`` on random
vectors, on the all-null case and on the everything-significant case.

*Partial Spearman with ties.*  The vectorized residual construction is pinned
against the textbook three-correlation formula computed with the package's own
``spearman_correlation``, on data with heavy ties, which is the regime these
histology scores actually live in.

*Vacuity.*  ``all([])`` is ``True``. A gate with no applicable condition must
return ``NO_APPLICABLE_CONDITIONS`` and must not report itself as passed.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np

from masld_bench.evaluators.metrics import spearman_correlation
from masld_bench.evaluators.stats import benjamini_hochberg


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "scripts" / "evaluate_gse267145_axis_count.py"
FREEZER = ROOT / "scripts" / "freeze_axis_count_prespecification.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


analysis = _load(MODULE, "axis_count_analysis_tested")
freezer = _load(FREEZER, "axis_count_freezer_tested")


class ThresholdsMatchTheFreezerTests(unittest.TestCase):
    def test_the_analysis_expects_exactly_the_literals_the_freezer_writes(self) -> None:
        self.assertEqual(freezer.FROZEN_THRESHOLDS, analysis.FROZEN_THRESHOLDS)

    def test_d1_is_the_decisive_criterion_in_both_files(self) -> None:
        self.assertEqual(
            sorted(analysis.FROZEN_THRESHOLDS),
            [
                "d1_activity_and_fibrosis_are_independent",
                "d2_steatosis_and_saf_activity_are_independent",
            ],
        )

    def test_the_two_composites_are_different_vectors_by_definition(self) -> None:
        self.assertEqual(analysis.NAS_ACTIVITY, freezer.NAS_ACTIVITY)
        self.assertEqual(analysis.SAF_ACTIVITY, freezer.SAF_ACTIVITY)
        self.assertIn("steatosis", analysis.NAS_ACTIVITY)
        self.assertNotIn("steatosis", analysis.SAF_ACTIVITY)

    def test_necrosis_enters_no_composite(self) -> None:
        for composite in (analysis.NAS_ACTIVITY, analysis.SAF_ACTIVITY):
            self.assertNotIn("lobular_necrosis", composite)
            self.assertNotIn("fibrosis", composite)

    def test_both_directions_of_each_criterion_are_declared(self) -> None:
        self.assertEqual(
            analysis.D1_FAMILIES,
            (("nas_activity_sum", "fibrosis"), ("fibrosis", "nas_activity_sum")),
        )
        self.assertEqual(
            analysis.D2_FAMILIES,
            (("steatosis", "saf_activity_sum"), ("saf_activity_sum", "steatosis")),
        )


class FastBenjaminiHochbergTests(unittest.TestCase):
    """The shortcut in the inner loop must equal the package's adjustment."""

    def _package_count(self, p_values: np.ndarray) -> int:
        adjusted = np.asarray(benjamini_hochberg(p_values.tolist()), dtype=float)
        return int(np.sum(adjusted < analysis.BH_ALPHA))

    def test_step_up_matches_the_package_on_random_p_values(self) -> None:
        generator = np.random.default_rng(11)
        for _ in range(20):
            p_values = generator.random(500)
            self.assertEqual(analysis.bh_count(p_values), self._package_count(p_values))

    def test_step_up_matches_the_package_when_signal_is_present(self) -> None:
        generator = np.random.default_rng(12)
        p_values = np.concatenate(
            [generator.random(400), generator.random(100) * 1e-6]
        )
        self.assertEqual(analysis.bh_count(p_values), self._package_count(p_values))

    def test_nothing_significant_gives_zero_not_a_silent_one(self) -> None:
        self.assertEqual(analysis.bh_count(np.full(300, 0.9)), 0)

    def test_everything_significant_gives_the_whole_family(self) -> None:
        self.assertEqual(analysis.bh_count(np.full(300, 1e-12)), 300)

    def test_an_empty_family_counts_zero(self) -> None:
        self.assertEqual(analysis.bh_count(np.asarray([], dtype=float)), 0)

    def test_the_abs_r_shortcut_equals_the_p_value_path(self) -> None:
        generator = np.random.default_rng(13)
        for df in (96, 97):
            for m in (50, 500, 5000):
                critical = analysis.bh_critical_abs_r(m, df)
                for scale in (0.05, 0.2, 0.5):
                    values = np.clip(generator.normal(0.0, scale, m), -0.99, 0.99)
                    self.assertEqual(
                        analysis.bh_count_from_abs_r(np.abs(values), critical),
                        analysis.bh_count(analysis.t_p_values(values, df)),
                    )

    def test_the_vectorized_block_equals_the_row_by_row_shortcut(self) -> None:
        generator = np.random.default_rng(14)
        critical = analysis.bh_critical_abs_r(400, 96)
        block = np.clip(generator.normal(0.0, 0.25, (37, 400)), -0.99, 0.99)
        vectorized = analysis.bh_counts_from_abs_r(block, critical)
        for row, count in zip(block, vectorized):
            self.assertEqual(
                int(count), analysis.bh_count_from_abs_r(np.abs(row), critical)
            )

    def test_a_block_with_no_rejection_anywhere_returns_zeros(self) -> None:
        critical = analysis.bh_critical_abs_r(400, 96)
        block = np.full((5, 400), 1e-6)
        self.assertTrue(np.all(analysis.bh_counts_from_abs_r(block, critical) == 0))


class PartialSpearmanTests(unittest.TestCase):
    """The residual construction must be the textbook partial rank correlation."""

    @staticmethod
    def _textbook(gene, exposure, covariate) -> float:
        r_ge = spearman_correlation(list(gene), list(exposure))
        r_gc = spearman_correlation(list(gene), list(covariate))
        r_ec = spearman_correlation(list(exposure), list(covariate))
        return (r_ge - r_gc * r_ec) / np.sqrt((1 - r_gc**2) * (1 - r_ec**2))

    def _vectorized(self, gene, exposure, covariate) -> float:
        unit_covariate = analysis.unit_ranks(np.asarray(covariate, dtype=float))
        gene_residual = analysis.residualize(
            analysis.centred_ranks(np.asarray(gene, dtype=float)), unit_covariate
        )
        exposure_residual = analysis.residualize(
            analysis.centred_ranks(np.asarray(exposure, dtype=float)), unit_covariate
        )
        return float(
            gene_residual
            @ exposure_residual
            / np.sqrt((gene_residual**2).sum() * (exposure_residual**2).sum())
        )

    def test_matches_the_textbook_formula_without_ties(self) -> None:
        generator = np.random.default_rng(21)
        for _ in range(20):
            covariate = generator.normal(size=60)
            exposure = covariate + generator.normal(size=60)
            gene = 0.5 * exposure + generator.normal(size=60)
            self.assertAlmostEqual(
                self._vectorized(gene, exposure, covariate),
                self._textbook(gene, exposure, covariate),
                places=10,
            )

    def test_matches_the_textbook_formula_with_heavy_ties(self) -> None:
        """The regime these histology scores actually live in."""

        generator = np.random.default_rng(22)
        for _ in range(20):
            covariate = generator.integers(0, 4, size=99).astype(float)
            covariate[covariate > 0] *= generator.integers(0, 2, size=int((covariate > 0).sum()))
            exposure = np.clip(covariate + generator.integers(0, 5, size=99), 0, 7).astype(float)
            gene = generator.normal(size=99) + 0.3 * exposure
            self.assertAlmostEqual(
                self._vectorized(gene, exposure, covariate),
                self._textbook(gene, exposure, covariate),
                places=10,
            )

    def test_a_gene_collinear_with_the_covariate_has_a_zero_residual(self) -> None:
        covariate = np.asarray([0, 0, 1, 1, 2, 2, 3, 3], dtype=float)
        unit_covariate = analysis.unit_ranks(covariate)
        residual = analysis.residualize(analysis.centred_ranks(covariate), unit_covariate)
        self.assertAlmostEqual(float(np.sqrt((residual**2).sum())), 0.0, places=12)

    def test_unit_columns_reports_a_zero_norm_rather_than_dividing_by_it(self) -> None:
        matrix = np.column_stack([np.asarray([1.0, -1.0, 0.0]), np.zeros(3)])
        scaled, norms = analysis.unit_columns(matrix)
        self.assertGreater(norms[0], 0.0)
        self.assertEqual(norms[1], 0.0)
        self.assertTrue(np.all(scaled[:, 1] == 0.0))

    def test_a_permuted_exposure_stays_orthogonal_to_the_covariate(self) -> None:
        generator = np.random.default_rng(23)
        covariate = generator.integers(0, 4, size=50).astype(float)
        exposure = generator.integers(0, 8, size=50).astype(float)
        unit_covariate = analysis.unit_ranks(covariate)
        residual = analysis.residualize(
            analysis.centred_ranks(exposure), unit_covariate
        )
        for _ in range(20):
            permuted = analysis._permuted_unit_exposure(
                residual, unit_covariate, generator.permutation(50)
            )
            self.assertIsNotNone(permuted)
            self.assertAlmostEqual(float(unit_covariate @ permuted), 0.0, places=12)
            self.assertAlmostEqual(float(permuted @ permuted), 1.0, places=12)


class StratifiedPermutationTests(unittest.TestCase):
    def test_a_within_stratum_order_never_leaves_its_stratum(self) -> None:
        strata = ["F"] * 85 + ["M"] * 14
        groups = list(analysis._stratum_index_groups(strata).values())
        generator = np.random.default_rng(31)
        for _ in range(50):
            order = analysis._stratum_orders(groups, generator, len(strata))
            self.assertEqual(sorted(order.tolist()), list(range(len(strata))))
            for index, source in enumerate(order):
                self.assertEqual(strata[index], strata[int(source)])

    def test_a_singleton_stratum_is_held_fixed(self) -> None:
        strata = ["a", "a", "b"]
        groups = list(analysis._stratum_index_groups(strata).values())
        generator = np.random.default_rng(32)
        for _ in range(20):
            order = analysis._stratum_orders(groups, generator, 3)
            self.assertEqual(int(order[2]), 2)


class VacuityTests(unittest.TestCase):
    """all([]) is True. A gate must never inherit that."""

    def test_no_applicable_condition_is_not_a_pass(self) -> None:
        gate = analysis.evaluate_gate(
            "d1",
            [
                {"condition": "a", "applicable": False, "met": False},
                {"condition": "b", "applicable": False, "met": False},
            ],
        )
        self.assertEqual(gate["verdict"], "NO_APPLICABLE_CONDITIONS")
        self.assertFalse(gate["passed"])
        self.assertEqual(gate["n_conditions_applicable"], 0)

    def test_an_empty_condition_list_is_not_a_pass(self) -> None:
        gate = analysis.evaluate_gate("d1", [])
        self.assertEqual(gate["verdict"], "NO_APPLICABLE_CONDITIONS")
        self.assertFalse(gate["passed"])

    def test_the_verdict_is_met_over_applicable_not_over_total(self) -> None:
        gate = analysis.evaluate_gate(
            "d1",
            [
                {"condition": "a", "applicable": True, "met": True},
                {"condition": "b", "applicable": False, "met": False},
            ],
        )
        self.assertTrue(gate["passed"])
        self.assertEqual(gate["n_conditions_total"], 2)
        self.assertEqual(gate["n_conditions_applicable"], 1)
        self.assertEqual(gate["n_conditions_met"], 1)

    def test_one_unmet_applicable_condition_fails_the_gate(self) -> None:
        gate = analysis.evaluate_gate(
            "d1",
            [
                {"condition": "a", "applicable": True, "met": True},
                {"condition": "b", "applicable": True, "met": False},
            ],
        )
        self.assertEqual(gate["verdict"], "NOT_MET")
        self.assertFalse(gate["passed"])

    def test_a_count_at_or_below_the_floor_does_not_meet_a_direction(self) -> None:
        for count, floor, expected in ((0, 0.0, False), (5, 5.0, False), (6, 5.0, True)):
            condition = analysis.direction_condition({
                "exposure": "x",
                "adjusted_for": "z",
                "family_size": 100,
                "genes_bh_below_0_05": count,
                "count_null": {
                    "null_percentile_95_count": floor,
                    "count_exceedance_p_value": 0.001,
                },
            })
            self.assertEqual(condition["met"], expected)

    def test_an_empty_family_is_inapplicable_rather_than_failed(self) -> None:
        condition = analysis.direction_condition({
            "exposure": "x",
            "adjusted_for": "z",
            "family_size": 0,
            "genes_bh_below_0_05": 0,
            "count_null": {
                "null_percentile_95_count": 0.0,
                "count_exceedance_p_value": 1.0,
            },
        })
        self.assertFalse(condition["applicable"])
        self.assertFalse(condition["met"])
        self.assertIsNotNone(condition["why_not_applicable"])


class FrozenDocumentTests(unittest.TestCase):
    """The freezer must record the provenance and the complete outcome map."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = freezer.build(_endpoint_rows())

    def test_the_provenance_disclosure_is_present_and_affirmative(self) -> None:
        provenance = self.payload["honest_provenance"]
        self.assertTrue(
            provenance["this_question_was_derived_after_seeing_stage_0s_answer"]
        )
        self.assertTrue(
            provenance["criteria_fixed_before_any_expression_value_was_read"]
        )
        self.assertEqual(provenance["stage_0_decision"], "STOP")
        self.assertEqual(
            provenance["stage_0_prespecification_artifacts_sha256"],
            freezer.STAGE0_PRESPEC_ARTIFACTS_SHA256,
        )

    def test_no_gene_expression_was_opened_before_the_freeze(self) -> None:
        self.assertTrue(
            self.payload["what_was_inspected_before_freezing"][
                "no_gene_expression_was_opened"
            ]
        )

    def test_the_outcome_map_covers_all_four_cells(self) -> None:
        outcomes = [
            cell["outcome"] for cell in self.payload["decision_rule"]["cells"]
        ]
        self.assertEqual(
            outcomes,
            [
                "THREE_AXIS_SAF",
                "TWO_AXIS_ACTIVITY_FIBROSIS",
                "TWO_AXIS_STEATOSIS_ACTIVITY",
                "ONE_AXIS",
            ],
        )

    def test_the_determination_audit_is_recorded_as_numbers(self) -> None:
        audit = self.payload["determination_audit_run_before_the_freeze"]
        fibrosis = audit["rank_fibrosis_in_the_span_of_the_three_aspects"]
        composite = audit[
            "rank_saf_activity_sum_in_the_span_of_its_two_components"
        ]
        self.assertLess(fibrosis, 0.6)
        self.assertGreater(composite, 0.99)

    def test_d1_is_marked_decisive_and_d2_is_not(self) -> None:
        criteria = self.payload["criteria"]
        self.assertTrue(criteria["d1_activity_and_fibrosis_are_independent"]["decisive"])
        self.assertFalse(
            criteria["d2_steatosis_and_saf_activity_are_independent"]["decisive"]
        )

    def test_d3_is_declared_a_diagnostic_and_never_a_gate(self) -> None:
        self.assertFalse(
            self.payload["diagnostics_never_gates"][
                "d3_does_a_composite_earn_its_place"
            ]["is_a_gate"]
        )


def _endpoint_rows() -> list[dict[str, str]]:
    path = (
        ROOT
        / "executions"
        / "model-data-061-21079623"
        / "activation"
        / "participant_endpoints.tsv"
    )
    return freezer.read_table(path)


if __name__ == "__main__":
    unittest.main()

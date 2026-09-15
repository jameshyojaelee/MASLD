from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.analyze_bulk_lineage_composition_fibrosis import (
    COLLAPSE_CONFOUNDING_ALPHA,
    DETECTION_FLOOR,
    MIN_FRACTION_ABOVE_FLOOR,
    collapse_confounding,
    benjamini_hochberg,
    fraction_above_floor,
    permutation_null,
    spearman,
)
from scripts.freeze_lineage_composition_substrate import (
    LINEAGES,
    STAGE_MAP,
    SubstrateFreezeError,
    read_proportions,
    read_rowname_offset_tsv,
)


class OffsetReaderTests(unittest.TestCase):
    def test_metadata_offset_required(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "m.tsv"
            path.write_text("a\tb\nrow\t1\t2\n", encoding="utf-8")
            header, rows = read_rowname_offset_tsv(path)
            self.assertEqual(header, ["a", "b"])
            self.assertEqual(rows[0]["a"], "1")
            path.write_text("a\tb\n1\t2\n", encoding="utf-8")
            with self.assertRaises(SubstrateFreezeError):
                read_rowname_offset_tsv(path)

    def test_proportions_header_is_sixteen_lineages_with_no_key_name(self) -> None:
        """The trap: 16 header names, 17 data fields, no key column name."""

        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "p.tsv"
            vals = [1.0 / 16] * 16
            path.write_text(
                "\t".join(LINEAGES) + "\n"
                + "SRR1\t" + "\t".join(f"{v:.17g}" for v in vals) + "\n",
                encoding="utf-8",
            )
            keys, matrix = read_proportions(path)
            self.assertEqual(keys, ["SRR1"])
            self.assertEqual(matrix.shape, (1, 16))

    def test_a_fifteen_lineage_parse_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "p.tsv"
            path.write_text(
                "key\t" + "\t".join(LINEAGES[:15]) + "\n"
                + "SRR1\t" + "\t".join(["0.0625"] * 15) + "\n",
                encoding="utf-8",
            )
            with self.assertRaises(SubstrateFreezeError):
                read_proportions(path)

    def test_rows_that_do_not_sum_to_one_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            path = Path(value) / "p.tsv"
            path.write_text(
                "\t".join(LINEAGES) + "\n" + "SRR1\t" + "\t".join(["0.5"] * 16) + "\n",
                encoding="utf-8",
            )
            with self.assertRaises(SubstrateFreezeError):
                read_proportions(path)


class StageMapTests(unittest.TestCase):
    def test_control_arms_map_to_none_never_to_stage_zero(self) -> None:
        self.assertIsNone(STAGE_MAP["GSE162694"]["Control"])
        self.assertIsNone(STAGE_MAP["GSE213621"]["Control"])
        self.assertEqual(STAGE_MAP["GSE162694"]["NASH_F0"], 0)

    def test_graded_cohorts_carry_all_five_levels(self) -> None:
        for cohort in ("GSE135251", "GSE162694", "GSE240729"):
            graded = {v for v in STAGE_MAP[cohort].values() if v is not None}
            self.assertEqual(graded, {0, 1, 2, 3, 4}, cohort)

    def test_the_pooled_cohort_is_three_levels_not_five(self) -> None:
        graded = {v for v in STAGE_MAP["GSE213621"].values() if v is not None}
        self.assertEqual(graded, {0, 1, 2})


class DetectionFloorTests(unittest.TestCase):
    def test_fraction_above_floor_counts_correctly(self) -> None:
        values = np.asarray([1e-2, 1e-3, 1e-5, 1e-22])
        self.assertAlmostEqual(fraction_above_floor(values), 0.5)

    def test_the_criterion_is_the_median_clearing_the_floor(self) -> None:
        """0.50 is exactly the point at which the median sample clears the floor."""

        just_under = np.asarray([1e-2, 1e-9, 1e-9, 1e-9])
        just_over = np.asarray([1e-2, 1e-2, 1e-9, 1e-9])
        self.assertLess(fraction_above_floor(just_under), MIN_FRACTION_ABOVE_FLOOR)
        self.assertGreaterEqual(fraction_above_floor(just_over), MIN_FRACTION_ABOVE_FLOOR)

    def test_a_collapsed_lineage_still_produces_a_finite_rho(self) -> None:
        """This is why the floor rule exists: rank makes noise finite, not informative."""

        rng = np.random.default_rng(0)
        collapsed = rng.uniform(1e-22, 1e-15, size=60)
        stage = np.repeat(np.arange(5.0), 12)
        value = spearman(collapsed, stage)
        self.assertTrue(np.isfinite(value))
        self.assertLess(fraction_above_floor(collapsed), MIN_FRACTION_ABOVE_FLOOR)

    def test_the_floor_is_one_ten_thousandth(self) -> None:
        self.assertEqual(DETECTION_FLOOR, 1e-4)


class NullTests(unittest.TestCase):
    def test_measured_null_agrees_with_the_analytic_reference(self) -> None:
        """Outcome ties do not move a rank-correlation null."""

        rng = np.random.default_rng(3)
        n = 120
        x = rng.normal(size=n)
        stage = np.repeat(np.arange(5.0), n // 5)
        null = permutation_null(x, stage, replicates=1500, seed=11)
        self.assertAlmostEqual(null["null_sd"], null["analytic_null_se"], delta=0.02)
        self.assertAlmostEqual(null["null_mean"], 0.0, delta=0.02)

    def test_a_heavily_tied_outcome_does_not_shrink_the_null(self) -> None:
        rng = np.random.default_rng(5)
        n = 100
        x = rng.normal(size=n)
        few_ties = rng.permutation(np.arange(float(n)))
        many_ties = np.repeat(np.arange(2.0), n // 2)
        a = permutation_null(x, few_ties, replicates=1200, seed=2)
        b = permutation_null(x, many_ties, replicates=1200, seed=2)
        self.assertAlmostEqual(a["null_sd"], b["null_sd"], delta=0.03)


class BenjaminiHochbergTests(unittest.TestCase):
    def test_monotone_and_bounded(self) -> None:
        q = benjamini_hochberg([0.001, 0.01, 0.04, 0.5])
        self.assertTrue(all(0.0 <= v <= 1.0 for v in q))
        self.assertEqual(q, sorted(q))

    def test_a_family_of_three_is_less_stringent_than_a_family_of_sixteen(self) -> None:
        small = benjamini_hochberg([0.02, 0.5, 0.6])
        large = benjamini_hochberg([0.02] + [0.5] * 15)
        self.assertLess(small[0], large[0])

    def test_an_empty_family_returns_empty(self) -> None:
        self.assertEqual(benjamini_hochberg([]), [])


class CollapseConfoundingTests(unittest.TestCase):
    """Distinctness settles the tie mechanism; this settles the confounding one."""

    STAGE = np.repeat(np.arange(5.0), 20)  # n=100, five stages

    def test_collapse_that_tracks_stage_is_flagged(self) -> None:
        """The failure mode: degraded high-stage tissue collapses more often."""

        rng = np.random.default_rng(1)
        values = np.where(
            self.STAGE >= 3.0,
            rng.uniform(1e-20, 1e-12, size=100),   # collapsed in advanced stages
            rng.uniform(1e-3, 1e-2, size=100),     # measured in early stages
        )
        result = collapse_confounding(
            values, self.STAGE, floor=DETECTION_FLOOR, replicates=800, seed=3
        )
        self.assertTrue(result["applicable"])
        self.assertTrue(result["confounded"])
        self.assertLess(result["permutation_p"], COLLAPSE_CONFOUNDING_ALPHA)
        self.assertGreater(abs(result["rho_collapse_vs_stage"]), 0.5)

    def test_collapse_unrelated_to_stage_is_not_flagged(self) -> None:
        rng = np.random.default_rng(2)
        collapsed = rng.permutation(np.arange(100) < 40)
        values = np.where(
            collapsed, rng.uniform(1e-20, 1e-12, size=100), rng.uniform(1e-3, 1e-2, size=100)
        )
        result = collapse_confounding(
            values, self.STAGE, floor=DETECTION_FLOOR, replicates=800, seed=5
        )
        self.assertTrue(result["applicable"])
        self.assertFalse(result["confounded"])

    def test_constant_collapse_status_is_not_testable(self) -> None:
        rng = np.random.default_rng(4)
        never = rng.uniform(1e-3, 1e-2, size=100)
        always = rng.uniform(1e-20, 1e-12, size=100)
        for values, fraction in ((never, 0.0), (always, 1.0)):
            result = collapse_confounding(
                values, self.STAGE, floor=DETECTION_FLOOR, replicates=200, seed=6
            )
            self.assertFalse(result["applicable"])
            self.assertFalse(result["confounded"])
            self.assertAlmostEqual(result["collapsed_fraction"], fraction)

    def test_a_lineage_can_pass_the_floor_rule_and_fail_the_confounding_test(self) -> None:
        """This is why both guards exist rather than one."""

        rng = np.random.default_rng(7)
        values = np.where(
            self.STAGE >= 4.0,
            rng.uniform(1e-20, 1e-12, size=100),
            rng.uniform(1e-3, 1e-2, size=100),
        )
        self.assertGreaterEqual(fraction_above_floor(values), MIN_FRACTION_ABOVE_FLOOR)
        result = collapse_confounding(
            values, self.STAGE, floor=DETECTION_FLOOR, replicates=800, seed=8
        )
        self.assertTrue(result["confounded"])

    def test_the_binary_indicator_null_is_narrower_than_a_continuous_one(self) -> None:
        """Why the null is measured against the indicator, not assumed."""

        from scripts.analyze_bulk_lineage_composition_fibrosis import permutation_null

        rng = np.random.default_rng(9)
        collapsed = (rng.permutation(np.arange(100)) < 10).astype(float)
        continuous = rng.normal(size=100)
        tied = permutation_null(collapsed, self.STAGE, replicates=800, seed=10)
        cont = permutation_null(continuous, self.STAGE, replicates=800, seed=10)
        self.assertLess(tied["null_p95_abs"], cont["null_p95_abs"])


class LineageAxisTests(unittest.TestCase):
    def test_there_are_sixteen_lineages(self) -> None:
        self.assertEqual(len(LINEAGES), 16)
        self.assertEqual(len(set(LINEAGES)), 16)
        self.assertIn("Neutrophils", LINEAGES)
        self.assertIn("Hepatocytes", LINEAGES)


class FrozenGateAlignmentTests(unittest.TestCase):
    """Regression test for the stale-key defect: the script must read the gate as frozen."""

    GATE = Path(
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/"
        "MASLD_Model_Benchmark/config/evaluation/"
        "bulk_lineage_composition_fibrosis_promotion_gate.json"
    )
    EVALUATOR_CONDITION_IDS = {
        "qualifying_lineage_count_minimum",
        "qualifying_set_dominance_guard",
        "held_out_magnitude_retention_minimum",
        "co_primary_set_change_maximum",
        "seed_direction_consistency",
    }

    def setUp(self) -> None:
        import json

        if not self.GATE.is_file():
            self.skipTest("frozen gate not present")
        self.gate = json.loads(self.GATE.read_text(encoding="utf-8"))

    def test_every_registered_condition_is_computed_and_no_others(self) -> None:
        registered = set(self.gate["development_advance_conditions"])
        self.assertEqual(
            registered,
            self.EVALUATOR_CONDITION_IDS,
            "a gate condition the evaluator does not compute would report as absent "
            "rather than fail, which is the defect this test exists to catch",
        )

    def test_the_effect_floor_is_a_filter_not_a_condition(self) -> None:
        self.assertNotIn(
            "qualifying_lineage_effect_floor", self.gate["development_advance_conditions"]
        )
        self.assertIn("effect_floor", self.gate["qualifying_set_filters"])
        self.assertEqual(
            self.gate["qualifying_set_filters"]["effect_floor"]["role"],
            "filter, not a scored condition",
        )

    def test_the_filters_the_evaluator_hardcodes_match_the_gate(self) -> None:
        filters = self.gate["qualifying_set_filters"]
        self.assertEqual(float(filters["detection_floor"]["floor"]), DETECTION_FLOOR)
        self.assertEqual(
            float(filters["detection_floor"]["minimum_fraction_above_floor"]),
            MIN_FRACTION_ABOVE_FLOOR,
        )
        self.assertEqual(
            float(filters["collapse_confounding"]["alpha"]), COLLAPSE_CONFOUNDING_ALPHA
        )

    def test_the_dominance_guard_names_the_lineage_the_evaluator_tests(self) -> None:
        from scripts.analyze_bulk_lineage_composition_fibrosis import DOMINANT_LINEAGE

        requirement = self.gate["development_advance_conditions"][
            "qualifying_set_dominance_guard"
        ]["requirement"]
        self.assertIn(DOMINANT_LINEAGE, requirement)


if __name__ == "__main__":
    unittest.main()

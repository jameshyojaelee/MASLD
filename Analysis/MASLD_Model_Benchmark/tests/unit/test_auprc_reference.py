"""Lock the two facts the ranking-metric null reference exists to state."""

from __future__ import annotations

import unittest

from masld_bench.evaluators.auprc_reference import (
    GUIDANCE,
    bias_table,
    permutation_reference,
    random_score_reference,
    spearman_reference,
    stratified_permutation_reference,
    stratified_spearman_reference,
    tie_structure_table,
)
from masld_bench.evaluators.metrics import MetricError


def _structured_case() -> tuple[list[int], list[float], list[str]]:
    """Four strata whose membership predicts both the label and the score.

    Within each stratum the score carries no information about the label, so an
    honest null must report no signal.  A global permutation destroys the
    stratum structure too and therefore credits the observed statistic for it.
    """

    labels: list[int] = []
    scores: list[float] = []
    strata: list[str] = []
    for stratum, (positive_rate, offset) in enumerate(
        ((0.9, 10.0), (0.7, 6.0), (0.3, 3.0), (0.1, 0.0))
    ):
        for index in range(20):
            label = 1 if index < int(round(positive_rate * 20)) else 0
            labels.append(label)
            # Score depends only on the stratum plus a within-stratum ramp that
            # is deliberately unrelated to the label ordering.
            scores.append(offset + (index % 5) * 0.01)
            strata.append(f"s{stratum}")
    return labels, scores, strata


class PrevalenceIsNotTheNullTests(unittest.TestCase):
    def test_a_random_scorer_beats_prevalence_at_every_small_size(self) -> None:
        for n_positive, n_negative in ((5, 15), (15, 22), (20, 20), (30, 70)):
            reference = random_score_reference(
                n_positive, n_negative, n_draws=2_000
            )
            self.assertGreater(
                reference.mean,
                reference.prevalence,
                f"{n_positive}/{n_negative} should show upward bias",
            )
            self.assertGreater(reference.bias_over_prevalence, 0.0)

    def test_the_bias_shrinks_as_the_cohort_grows(self) -> None:
        small = random_score_reference(5, 15, n_draws=4_000)
        medium = random_score_reference(25, 75, n_draws=4_000)
        large = random_score_reference(100, 300, n_draws=4_000)
        self.assertAlmostEqual(small.prevalence, medium.prevalence, places=6)
        self.assertAlmostEqual(medium.prevalence, large.prevalence, places=6)
        self.assertGreater(small.bias_over_prevalence, medium.bias_over_prevalence)
        self.assertGreater(medium.bias_over_prevalence, large.bias_over_prevalence)

    def test_the_gap_that_matters_is_the_upper_tail_not_the_mean(self) -> None:
        """A one-in-twenty random result is far above prevalence at n=37."""

        reference = random_score_reference(15, 22, n_draws=10_000)
        self.assertGreater(reference.percentile_95 - reference.prevalence, 0.15)

    def test_the_campaign_configuration_reproduces(self) -> None:
        reference = random_score_reference(15, 22, n_draws=10_000)
        self.assertAlmostEqual(reference.prevalence, 15 / 37, places=6)
        self.assertAlmostEqual(reference.mean, 0.458, places=2)
        self.assertAlmostEqual(reference.percentile_95, 0.61, places=2)


class TieStructureTests(unittest.TestCase):
    def test_fewer_distinct_scores_give_a_tighter_null(self) -> None:
        rows = tie_structure_table(15, 22, [2, 3, 5, 10, None], n_draws=4_000)
        spreads = [row["null_percentile_95"] for row in rows]
        self.assertEqual(spreads, sorted(spreads), "null must widen with k")
        self.assertLess(rows[0]["null_percentile_95"], rows[-1]["null_percentile_95"])

    def test_a_binary_scorer_cannot_reach_the_continuous_floor(self) -> None:
        binary = random_score_reference(
            15, 22, n_draws=4_000, distinct_score_values=2
        )
        continuous = random_score_reference(15, 22, n_draws=4_000)
        self.assertLess(binary.percentile_95, continuous.percentile_95)
        self.assertLess(binary.sd, continuous.sd)

    def test_permutation_reference_records_the_observed_tie_count(self) -> None:
        labels = [1] * 5 + [0] * 7
        scores = [1.0] * 5 + [0.0] * 7  # a perfectly separating binary scorer
        reference, observed, p_value = permutation_reference(
            labels, scores, n_permutations=2_000
        )
        self.assertEqual(reference.distinct_score_values, 2)
        self.assertAlmostEqual(observed, 1.0, places=6)
        self.assertLess(p_value, 0.05)

    def test_a_degenerate_scorer_is_rejected(self) -> None:
        with self.assertRaises(MetricError):
            random_score_reference(5, 5, n_draws=10, distinct_score_values=1)


class GuardTests(unittest.TestCase):
    def test_both_classes_are_required(self) -> None:
        with self.assertRaises(MetricError):
            random_score_reference(0, 10, n_draws=10)
        with self.assertRaises(MetricError):
            random_score_reference(10, 0, n_draws=10)

    def test_spearman_reference_respects_outcome_ties(self) -> None:
        spread = spearman_reference(list(range(38)), n_draws=4_000)
        tied = spearman_reference([0.0] * 20 + list(range(18)), n_draws=4_000)
        self.assertGreater(spread.percentile_95, 0.2)
        self.assertNotAlmostEqual(
            spread.percentile_95, tied.percentile_95, places=3
        )

    def test_guidance_names_both_references(self) -> None:
        self.assertIn("random-score reference", GUIDANCE)
        self.assertIn("never against class prevalence", GUIDANCE)
        self.assertIn("permutation", GUIDANCE)

    def test_bias_table_is_lane_agnostic(self) -> None:
        rows = bias_table([(15, 22), (3, 30)], n_draws=1_000)
        self.assertEqual(len(rows), 2)
        self.assertEqual(set(rows[0]) & {"dataset_id", "donor_id"}, set())


class StratifiedNullTests(unittest.TestCase):
    """A global permutation destroys nuisance structure along with the signal."""

    def test_the_stratified_null_is_harder_when_structure_exists(self) -> None:
        labels, scores, strata = _structured_case()
        globally, observed, global_p = permutation_reference(
            labels, scores, n_permutations=2000, seed=11
        )
        within, observed_within, stratified_p = stratified_permutation_reference(
            labels, scores, strata, n_permutations=2000, seed=11
        )
        self.assertAlmostEqual(observed, observed_within, places=12)
        # The stratified null keeps the between-stratum association, so it sits
        # far higher and the observed value is much less impressive against it.
        self.assertGreater(within.mean, globally.mean)
        self.assertGreater(stratified_p, global_p)

    def test_the_two_nulls_converge_when_strata_carry_no_structure(self) -> None:
        """The control: the gap must come from structure, not the code path."""

        labels, scores, _ = _structured_case()
        shuffled = [f"r{index % 4}" for index in range(len(labels))]
        globally, _, _ = permutation_reference(
            labels, scores, n_permutations=2000, seed=23
        )
        within, _, _ = stratified_permutation_reference(
            labels, scores, shuffled, n_permutations=2000, seed=23
        )
        self.assertAlmostEqual(within.mean, globally.mean, delta=0.05)

    def test_the_null_names_itself_distinctly(self) -> None:
        labels, scores, strata = _structured_case()
        within, _, _ = stratified_permutation_reference(
            labels, scores, strata, n_permutations=200, seed=5
        )
        self.assertEqual(within.null, "label_permutation_within_strata")
        self.assertEqual(within.n_strata, 4)
        self.assertEqual(within.non_contributing_strata, 0)

    def test_non_contributing_strata_are_counted_not_dropped(self) -> None:
        labels = [1, 0, 1, 0, 1, 1]
        scores = [0.9, 0.1, 0.8, 0.2, 0.7, 0.6]
        strata = ["a", "a", "b", "b", "c", "c"]
        within, _, _ = stratified_permutation_reference(
            labels, scores, strata, n_permutations=200, seed=7
        )
        self.assertEqual(within.n_strata, 3)
        self.assertEqual(within.non_contributing_strata, 1)

    def test_all_degenerate_strata_fail_closed(self) -> None:
        with self.assertRaises(MetricError):
            stratified_permutation_reference(
                [1, 1, 0, 0], [0.9, 0.8, 0.2, 0.1], ["a", "a", "b", "b"],
                n_permutations=50,
            )

    def test_misaligned_strata_are_rejected(self) -> None:
        with self.assertRaises(MetricError):
            stratified_permutation_reference(
                [1, 0, 1], [0.9, 0.1, 0.5], ["a", "b"], n_permutations=50
            )

    def test_stratified_spearman_holds_between_stratum_structure_fixed(self) -> None:
        outcome = [float(index % 5) + 10.0 * (index // 20) for index in range(60)]
        scores = [float(index // 20) for index in range(60)]
        strata = [f"s{index // 20}" for index in range(60)]
        within, observed, p_value = stratified_spearman_reference(
            outcome, scores, strata, n_draws=500, seed=3
        )
        self.assertEqual(within.null, "outcome_permutation_within_strata")
        self.assertEqual(within.n_strata, 3)
        # The score is constant inside every stratum, so once the outcome is
        # permuted within strata nothing distinguishes the observed value.
        self.assertGreater(p_value, 0.5)
        self.assertEqual(within.percentile_scale, "absolute_two_sided")
        self.assertTrue(-1.0 <= observed <= 1.0)

    def test_constant_outcome_strata_fail_closed(self) -> None:
        with self.assertRaises(MetricError):
            stratified_spearman_reference(
                [1.0, 1.0, 2.0, 2.0], [0.1, 0.2, 0.3, 0.4], ["a", "a", "b", "b"],
                n_draws=50,
            )


class BackwardCompatibilityTests(unittest.TestCase):
    def test_to_dict_keeps_every_existing_key(self) -> None:
        reference = random_score_reference(15, 22, n_draws=200, seed=1)
        payload = reference.to_dict()
        for key in (
            "metric", "null", "percentile_scale", "null_mean", "null_sd",
            "null_percentile_95", "null_percentile_99", "null_maximum",
            "n_draws", "seed", "n_positive", "n_negative", "prevalence",
            "distinct_score_values", "bias_over_prevalence",
        ):
            self.assertIn(key, payload)

    def test_a_global_reference_reports_no_strata(self) -> None:
        reference = random_score_reference(15, 22, n_draws=200, seed=1)
        self.assertIsNone(reference.n_strata)
        self.assertEqual(reference.non_contributing_strata, 0)

    def test_guidance_names_the_stratified_rule(self) -> None:
        self.assertIn("stratum", GUIDANCE)


if __name__ == "__main__":
    unittest.main()

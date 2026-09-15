from __future__ import annotations

import unittest
import unittest.mock

import numpy as np

from scripts.evaluate_gse83452_baseline_nash_transfer import (
    CLASSES,
    MANDATORY_BASELINES,
    METRIC_IDS,
    PRIMARY_METRIC,
    brier,
    evaluate_gate,
    macro_f1,
    metrics,
    permutation_null,
    seed_direction_consistency,
)

GATE = {
    "gate_id": "gse83452_baseline_nash_external_development_v1",
    "primary_endpoint": "participant_macro_f1_nash_status",
    "champion_eligible": False,
    "external_development_only": True,
    "development_advance_conditions": {
        "brier_score_degradation_maximum": 0.01,
        "gain_direction_required_bootstrap_lower_bound_above_zero": True,
        "minimum_absolute_macro_f1_gain_over_strongest_baseline": 0.02,
        "minimum_seed_direction_consistency": "4_of_5",
    },
}


class MetricTests(unittest.TestCase):
    def test_macro_f1_of_a_constant_all_positive_scorer(self) -> None:
        """The training-prevalence baseline calls every participant NASH."""

        observed = np.asarray([1] * 104 + [0] * 44)
        scores = np.full(148, 131 / 172)
        # F1 for NASH = 2*104/(2*104+44); F1 for no-NASH = 0.
        expected = float(np.mean([0.0, 2 * 104 / (2 * 104 + 44)]))
        self.assertAlmostEqual(macro_f1(observed, scores), expected)

    def test_macro_f1_of_a_perfect_scorer_is_one(self) -> None:
        observed = np.asarray([1, 1, 0, 0])
        self.assertAlmostEqual(
            macro_f1(observed, np.asarray([0.9, 0.8, 0.1, 0.2])), 1.0
        )

    def test_brier_is_on_the_nash_indicator(self) -> None:
        observed = np.asarray([1, 0])
        self.assertAlmostEqual(brier(observed, np.asarray([0.75, 0.25])), 0.0625)

    def test_the_metric_roster_matches_the_taskspec(self) -> None:
        observed = np.asarray([1, 1, 0, 0])
        value = metrics(observed, np.asarray([0.9, 0.4, 0.6, 0.1]))
        self.assertEqual(sorted(value), sorted(METRIC_IDS))
        self.assertEqual(PRIMARY_METRIC, "participant_macro_f1_nash_status")

    def test_the_mandatory_roster_is_the_taskspec_baseline_roster(self) -> None:
        self.assertEqual(
            sorted(MANDATORY_BASELINES),
            [
                "age_sex_logistic",
                "gene_median_elastic_net",
                "gene_median_linear_svm",
                "gene_median_pca_elastic_net",
                "per_array_rank_elastic_net",
                "training_prevalence",
            ],
        )


class NullTests(unittest.TestCase):
    def test_the_auprc_null_mean_exceeds_prevalence(self) -> None:
        """Scoring average precision against prevalence manufactures significance."""

        rng = np.random.default_rng(19)
        observed = np.asarray([1] * 104 + [0] * 44)
        scores = rng.normal(size=148)
        null = permutation_null(observed, scores, replicates=2000, seed=5)
        prevalence = float(np.mean(observed))
        self.assertGreater(null["auprc_null_mean"], prevalence)
        self.assertGreater(null["auprc_null_p95"], null["auprc_null_mean"])

    def test_a_constant_scorer_has_a_degenerate_macro_f1_null(self) -> None:
        observed = np.asarray([1] * 104 + [0] * 44)
        constant = np.full(148, 0.76)
        continuous = np.random.default_rng(3).normal(size=148)
        flat = permutation_null(observed, constant, replicates=500, seed=5)
        wide = permutation_null(observed, continuous, replicates=500, seed=5)
        self.assertEqual(flat["distinct_score_values"], 1)
        self.assertLessEqual(flat["macro_f1_null_p95"], wide["macro_f1_null_p95"])

    def test_both_the_null_and_the_prevalence_reference_are_emitted(self) -> None:
        observed = np.asarray([1] * 10 + [0] * 10)
        null = permutation_null(observed, np.arange(20.0), replicates=200, seed=1)
        for key in (
            "auprc_null_mean",
            "auprc_null_p95",
            "auprc_permutation_p_value",
            "macro_f1_null_mean",
            "macro_f1_null_p95",
            "macro_f1_permutation_p_value",
            "auroc_null_mean",
        ):
            self.assertIn(key, null)


class SeedConsistencyTests(unittest.TestCase):
    def test_five_identical_seed_vectors_carry_no_stability_evidence(self) -> None:
        observed = np.asarray([1, 1, 0, 0])
        seeds = np.tile(np.asarray([0.9, 0.8, 0.1, 0.2]), (5, 1))
        summary = seed_direction_consistency(observed, seeds, 0.0)
        self.assertEqual(summary["distinct_seed_prediction_vectors"], 1)
        self.assertTrue(summary["deterministic_single_fit"])
        self.assertFalse(summary["consistency_is_evidence_of_stability"])
        self.assertEqual(summary["seeds_above_comparator"], 5)

    def test_distinct_seed_vectors_are_a_real_replicate(self) -> None:
        observed = np.asarray([1, 1, 0, 0])
        seeds = np.asarray(
            [
                [0.9, 0.8, 0.1, 0.2],
                [0.7, 0.6, 0.3, 0.4],
                [0.9, 0.7, 0.2, 0.1],
                [0.8, 0.9, 0.1, 0.3],
                [0.6, 0.7, 0.4, 0.2],
            ]
        )
        summary = seed_direction_consistency(observed, seeds, 0.0)
        self.assertEqual(summary["distinct_seed_prediction_vectors"], 5)
        self.assertTrue(summary["consistency_is_evidence_of_stability"])


class GateTests(unittest.TestCase):
    def _verdict(self, *, gain, ci_low, brier_change, distinct, above):
        candidate = {PRIMARY_METRIC: 0.60, "participant_brier_nash_status": 0.20}
        comparator = {
            PRIMARY_METRIC: 0.60 - gain,
            "participant_brier_nash_status": 0.20 - brier_change,
        }
        summary = {
            "seeds": 5,
            "seeds_above_comparator": above,
            "distinct_seed_prediction_vectors": distinct,
            "consistency_is_evidence_of_stability": distinct > 1,
            "interpretation": "n/a",
        }
        return evaluate_gate(
            gate=GATE,
            candidate_id="per_array_rank_elastic_net",
            comparator_id="training_prevalence",
            candidate_metrics=candidate,
            comparator_metrics=comparator,
            paired_gain={"ci_low": ci_low, "ci_high": ci_low + 0.1},
            seed_summary=summary,
        )

    def test_all_four_conditions_must_hold_to_pass(self) -> None:
        verdict = self._verdict(
            gain=0.05, ci_low=0.01, brier_change=0.005, distinct=5, above=5
        )
        self.assertEqual(verdict["verdict"], "PASS")
        self.assertEqual(verdict["conditions_met"], 4)

    def test_a_gain_below_the_threshold_fails(self) -> None:
        verdict = self._verdict(
            gain=0.01, ci_low=0.005, brier_change=0.0, distinct=5, above=5
        )
        self.assertEqual(verdict["verdict"], "FAIL")

    def test_a_bootstrap_bound_at_or_below_zero_fails(self) -> None:
        verdict = self._verdict(
            gain=0.05, ci_low=0.0, brier_change=0.0, distinct=5, above=5
        )
        self.assertEqual(verdict["verdict"], "FAIL")

    def test_brier_degradation_above_the_cap_fails(self) -> None:
        verdict = self._verdict(
            gain=0.05, ci_low=0.01, brier_change=0.02, distinct=5, above=5
        )
        self.assertEqual(verdict["verdict"], "FAIL")

    def test_a_deterministic_single_fit_cannot_satisfy_seed_consistency(self) -> None:
        """5-of-5 by construction is not 4-of-5 stability evidence."""

        verdict = self._verdict(
            gain=0.05, ci_low=0.01, brier_change=0.0, distinct=1, above=5
        )
        self.assertEqual(verdict["verdict"], "FAIL")
        condition = next(
            item
            for item in verdict["conditions"]
            if item["condition_id"] == "minimum_seed_direction_consistency"
        )
        self.assertFalse(condition["met"])

    def test_the_gate_never_grants_a_champion_claim(self) -> None:
        verdict = self._verdict(
            gain=0.5, ci_low=0.4, brier_change=-0.1, distinct=5, above=5
        )
        self.assertFalse(verdict["champion_eligible"])
        self.assertTrue(verdict["external_development_only"])


class ClassOrderTests(unittest.TestCase):
    def test_the_positive_class_index_is_nash(self) -> None:
        self.assertEqual(CLASSES, ("no_nash", "nash"))
        self.assertEqual(CLASSES.index("nash"), 1)


if __name__ == "__main__":
    unittest.main()

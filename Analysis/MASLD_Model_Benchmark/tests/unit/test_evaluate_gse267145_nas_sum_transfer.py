from __future__ import annotations

import unittest

import numpy as np

from scripts.evaluate_gse267145_nas_sum_transfer import (
    PRIMARY_METRIC,
    calibration,
    evaluate_gate,
    mean_absolute_error,
    min_pairwise_correlation,
    permutation_null,
    spearman,
)

GATE = {
    "gate_id": "gse267145_nas_sum_external_development_v1",
    "primary_endpoint": PRIMARY_METRIC,
    "external_development_only": True,
    "champion_eligible": False,
    "development_advance_conditions": {
        "candidate_must_exceed_its_own_permutation_null_p95": True,
        "gain_direction_required_bootstrap_lower_bound_above_zero": True,
        "minimum_absolute_spearman_gain_over_strongest_baseline": 0.05,
        "minimum_seed_direction_consistency": "4_of_5",
    },
    "seed_direction_applicability": {"seed_determinism_correlation_threshold": 0.999},
}


def seed_summary(*, corr: float, above: int = 5):
    return {"seeds": 5, "min_pairwise_correlation": corr,
            "seed_spearman": [0.5] * 5, "seeds_above_comparator": above}


def run_gate(*, corr, gain=0.10, ci_low=0.02, cand=0.60, p95=0.17):
    return evaluate_gate(
        gate=GATE, candidate_id="c", comparator_id="b",
        candidate_metrics={PRIMARY_METRIC: cand},
        comparator_metrics={PRIMARY_METRIC: cand - gain},
        paired_gain={"ci_low": ci_low, "ci_high": 0.3},
        seed_summary=seed_summary(corr=corr),
        candidate_null={"null_p95": p95, "permutation_p_value": 1e-4},
    )


class SeedApplicabilityTests(unittest.TestCase):
    """The required fix: a deterministic fit must never earn the seed condition."""

    def test_deterministic_fit_records_not_applicable(self) -> None:
        v = run_gate(corr=1.0)
        seed = next(c for c in v["conditions"]
                    if c["condition_id"] == "minimum_seed_direction_consistency")
        self.assertEqual(seed["state"], "not_applicable")
        self.assertIn("deterministic", seed["reason"])
        self.assertEqual(v["conditions_applicable"], 3)
        self.assertEqual(v["conditions_total"], 4)
        self.assertEqual(v["verdict_form"], "3 of 3 applicable (4 total)")

    def test_near_deterministic_fit_also_not_applicable(self) -> None:
        v = run_gate(corr=0.99996)
        seed = next(c for c in v["conditions"]
                    if c["condition_id"] == "minimum_seed_direction_consistency")
        self.assertEqual(seed["state"], "not_applicable")

    def test_a_deterministic_fit_can_still_pass_on_the_other_three(self) -> None:
        self.assertEqual(run_gate(corr=1.0)["verdict"], "PASS")

    def test_genuinely_varying_seeds_make_the_condition_live(self) -> None:
        v = run_gate(corr=0.90)
        seed = next(c for c in v["conditions"]
                    if c["condition_id"] == "minimum_seed_direction_consistency")
        self.assertEqual(seed["state"], "met")
        self.assertEqual(v["conditions_applicable"], 4)
        self.assertEqual(v["verdict_form"], "4 of 4 applicable (4 total)")

    def test_a_live_seed_condition_can_fail(self) -> None:
        v = evaluate_gate(
            gate=GATE, candidate_id="c", comparator_id="b",
            candidate_metrics={PRIMARY_METRIC: 0.60},
            comparator_metrics={PRIMARY_METRIC: 0.50},
            paired_gain={"ci_low": 0.02, "ci_high": 0.3},
            seed_summary=seed_summary(corr=0.90, above=2),
            candidate_null={"null_p95": 0.17, "permutation_p_value": 1e-4},
        )
        self.assertEqual(v["verdict"], "FAIL")
        self.assertEqual(v["conditions_met"], 3)
        self.assertEqual(v["conditions_applicable"], 4)


class GateConditionTests(unittest.TestCase):
    def test_gain_below_threshold_fails(self) -> None:
        self.assertEqual(run_gate(corr=1.0, gain=0.02)["verdict"], "FAIL")

    def test_bootstrap_bound_at_or_below_zero_fails(self) -> None:
        self.assertEqual(run_gate(corr=1.0, ci_low=0.0)["verdict"], "FAIL")

    def test_candidate_inside_its_own_null_fails(self) -> None:
        v = run_gate(corr=1.0, cand=0.15, p95=0.17)
        cond = next(c for c in v["conditions"]
                    if c["condition_id"] == "candidate_must_exceed_its_own_permutation_null_p95")
        self.assertEqual(cond["state"], "not_met")
        self.assertEqual(v["verdict"], "FAIL")


class MetricTests(unittest.TestCase):
    def test_spearman_is_invariant_to_the_severity_shift(self) -> None:
        rng = np.random.default_rng(3)
        truth = rng.integers(0, 8, size=99).astype(float)
        pred = truth + rng.normal(size=99)
        self.assertAlmostEqual(spearman(truth, pred), spearman(truth, pred + 1.688), places=12)

    def test_mae_is_not_invariant_to_the_severity_shift(self) -> None:
        truth = np.zeros(10)
        self.assertLess(mean_absolute_error(truth, np.zeros(10)),
                        mean_absolute_error(truth, np.zeros(10) + 1.688))

    def test_perfect_calibration_gives_slope_one_intercept_zero(self) -> None:
        rng = np.random.default_rng(5)
        truth = rng.normal(size=200)
        slope, intercept = calibration(truth, truth)
        self.assertAlmostEqual(slope, 1.0, places=10)
        self.assertAlmostEqual(intercept, 0.0, places=10)

    def test_constant_prediction_gives_undefined_calibration(self) -> None:
        """np.var of a constant vector is 7.9e-31, not 0; the guard must be exact.

        Dividing a zero covariance by that residue reports a calibration line
        for the training-mean baseline, which has none.
        """

        self.assertGreater(float(np.var(np.full(10, 4.34))), 0.0)
        slope, intercept = calibration(np.arange(10.0), np.full(10, 4.34))
        self.assertTrue(np.isnan(slope) and np.isnan(intercept))

    def test_the_real_training_mean_constant_is_also_undefined(self) -> None:
        slope, intercept = calibration(
            np.arange(99.0) % 8, np.full(99, 4.344444444444444)
        )
        self.assertTrue(np.isnan(slope) and np.isnan(intercept))

    def test_null_p95_matches_the_analytic_reference(self) -> None:
        rng = np.random.default_rng(7)
        truth = rng.integers(0, 8, size=99).astype(float)
        null = permutation_null(truth, rng.normal(size=99), replicates=4000, seed=1)
        self.assertAlmostEqual(null["analytic_null_sd"], 1 / np.sqrt(98), places=12)
        self.assertAlmostEqual(null["null_p95"], 1.645 / np.sqrt(98), delta=0.03)

    def test_p_value_floor_is_reported(self) -> None:
        rng = np.random.default_rng(9)
        truth = rng.integers(0, 8, size=99).astype(float)
        null = permutation_null(truth, truth, replicates=2000, seed=1)
        self.assertAlmostEqual(null["p_value_resolution_floor"], 1 / 2001, places=12)
        self.assertGreaterEqual(null["permutation_p_value"], null["p_value_resolution_floor"])

    def test_min_pairwise_correlation_reports_the_minimum(self) -> None:
        rng = np.random.default_rng(11)
        base = rng.normal(size=99)
        stack = np.stack([base, base, base, base, rng.normal(size=99)])
        self.assertLess(min_pairwise_correlation(stack), 0.999)


if __name__ == "__main__":
    unittest.main()

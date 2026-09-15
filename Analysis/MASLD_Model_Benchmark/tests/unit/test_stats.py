"""Tests for inferential-unit-aware statistical helpers."""

from __future__ import annotations

from inspect import signature
import unittest

from masld_bench.evaluators.stats import (
    DEFAULT_BOOTSTRAP_SEED,
    DEFAULT_TRAINING_SEEDS,
    empirical_bootstrap_power_gate,
    five_seed_direction_stability,
    gain_concentration,
    ld_block_bootstrap,
    ld_block_paired_bootstrap,
    two_way_donor_block_bootstrap,
)


class EmpiricalPowerGateTests(unittest.TestCase):
    def test_endpoint_bootstrap_power_is_deterministic_and_formula_agnostic(self) -> None:
        # A frozen endpoint-recomputed error distribution, not row-level
        # observations. Its asymmetric values exercise the empirical tail.
        errors = [(-0.02 + (index % 11) * 0.004) for index in range(1_000)]
        observed_effect = 0.03
        distribution = [observed_effect + error for error in errors]
        first = empirical_bootstrap_power_gate(
            paired_difference_bootstrap=distribution,
            observed_effect=observed_effect,
            n_units=48,
            minimum_effect=0.05,
            n_primary_claims=5,
        )
        second = empirical_bootstrap_power_gate(
            paired_difference_bootstrap=distribution,
            observed_effect=observed_effect,
            n_units=48,
            minimum_effect=0.05,
            n_primary_claims=5,
        )
        self.assertEqual(first, second)
        self.assertEqual(first.method_id, "empirical_paired_endpoint_bootstrap_v1")
        self.assertAlmostEqual(first.effective_alpha, 0.01)
        self.assertEqual(first.n_resamples, 1_000)
        self.assertGreaterEqual(first.estimated_power, 0.80)
        self.assertTrue(first.passed)

    def test_power_gate_fails_closed_on_small_n_or_invalid_distribution(self) -> None:
        result = empirical_bootstrap_power_gate(
            paired_difference_bootstrap=[0.03] * 100,
            observed_effect=0.03,
            n_units=5,
            min_units=10,
            minimum_effect=0.05,
        )
        self.assertFalse(result.passed)
        self.assertIn("below the locked floor", result.reason)
        with self.assertRaisesRegex(ValueError, "at least 100"):
            empirical_bootstrap_power_gate(
                paired_difference_bootstrap=[0.03] * 99,
                observed_effect=0.03,
                n_units=48,
                minimum_effect=0.05,
            )


class TwoWayBootstrapTests(unittest.TestCase):
    def test_two_way_bootstrap_is_deterministic_and_cell_balanced(self) -> None:
        candidate = [1.0, 1.0, 2.0, -1.0, 0.0, 3.0, 1.0]
        baseline = [0.0] * len(candidate)
        donors = ["d1", "d1", "d1", "d2", "d2", "d3", "d3"]
        blocks = ["b1", "b1", "b2", "b1", "b2", "b1", "b2"]
        first = two_way_donor_block_bootstrap(
            candidate,
            baseline,
            donors,
            blocks,
            n_resamples=1_000,
            seed=4317,
            expected_n_donors=3,
            expected_n_genomic_blocks=2,
        )
        second = two_way_donor_block_bootstrap(
            candidate,
            baseline,
            donors,
            blocks,
            n_resamples=1_000,
            seed=4317,
            expected_n_donors=3,
            expected_n_genomic_blocks=2,
        )

        self.assertEqual(first, second)
        self.assertAlmostEqual(first.estimate, 1.0)
        self.assertEqual(first.n_donors, 3)
        self.assertEqual(first.n_genomic_blocks, 2)
        self.assertEqual(first.n_cells, 6)
        self.assertEqual(first.n_observations, 7)

    def test_two_way_bootstrap_validates_alignment_counts_and_crossing(self) -> None:
        with self.assertRaisesRegex(ValueError, "equal lengths"):
            two_way_donor_block_bootstrap(
                [1.0, 2.0], [0.0], ["d1", "d2"], ["b1", "b1"]
            )
        with self.assertRaisesRegex(ValueError, "independent donors"):
            two_way_donor_block_bootstrap(
                [1.0, 2.0],
                [0.0, 0.0],
                ["d1", "d1"],
                ["b1", "b2"],
            )
        with self.assertRaisesRegex(ValueError, "complete donor-by-genomic-block"):
            two_way_donor_block_bootstrap(
                [1.0, 2.0, 3.0],
                [0.0, 0.0, 0.0],
                ["d1", "d1", "d2"],
                ["b1", "b2", "b1"],
            )
        with self.assertRaisesRegex(ValueError, "does not match expected 3"):
            two_way_donor_block_bootstrap(
                [1.0, 2.0, 3.0, 4.0],
                [0.0, 0.0, 0.0, 0.0],
                ["d1", "d1", "d2", "d2"],
                ["b1", "b2", "b1", "b2"],
                expected_n_donors=3,
            )


class LDBlockBootstrapTests(unittest.TestCase):
    def test_ld_block_specialization_weights_blocks_not_rows(self) -> None:
        first = ld_block_paired_bootstrap(
            [1.0, 1.0, 4.0, 0.0, 0.0, 0.0],
            [0.0] * 6,
            ["b1", "b1", "b2", "b3", "b3", "b3"],
            n_resamples=1_000,
            seed=97,
            expected_n_ld_blocks=3,
        )
        second = ld_block_bootstrap(
            [1.0, 1.0, 4.0, 0.0, 0.0, 0.0],
            [0.0] * 6,
            ["b1", "b1", "b2", "b3", "b3", "b3"],
            n_resamples=1_000,
            seed=97,
            expected_n_ld_blocks=3,
        )

        self.assertEqual(first, second)
        self.assertAlmostEqual(first.estimate, 5.0 / 3.0)
        self.assertEqual(first.n_clusters, 3)
        self.assertEqual(first.n_observations, 6)

    def test_ld_block_specialization_rejects_invalid_unit_counts(self) -> None:
        with self.assertRaisesRegex(ValueError, "independent LD blocks"):
            ld_block_paired_bootstrap([1.0, 2.0], [0.0, 0.0], ["b1", "b1"])
        with self.assertRaisesRegex(ValueError, "does not match expected 3"):
            ld_block_paired_bootstrap(
                [1.0, 2.0],
                [0.0, 0.0],
                ["b1", "b2"],
                expected_n_ld_blocks=3,
            )
        with self.assertRaisesRegex(ValueError, "equal lengths"):
            ld_block_paired_bootstrap(
                [1.0], [0.0, 0.0], ["b1", "b2"]
            )


class SeedDirectionStabilityTests(unittest.TestCase):
    def test_four_of_five_passes_without_creating_replicates(self) -> None:
        result = five_seed_direction_stability(
            [0.2, 0.1, -0.1, 0.3, 0.2],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        )

        self.assertTrue(result.passed)
        self.assertEqual(result.qualifying_seed_count, 4)
        self.assertEqual(result.evaluated_seed_count, 5)
        self.assertEqual(result.required_seed_count, 4)
        self.assertEqual(result.seed_ids, DEFAULT_TRAINING_SEEDS)
        self.assertFalse(result.seeds_are_biological_replicates)
        self.assertFalse(result.supports_inferential_test)
        self.assertIn("not biological replicates", result.reason)

    def test_three_of_five_and_ties_fail(self) -> None:
        result = five_seed_direction_stability(
            [1.0, 1.0, 1.0, 0.0, -1.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        )
        self.assertFalse(result.passed)
        self.assertEqual(result.qualifying_seed_count, 3)

    def test_direction_stability_validates_exact_alignment_and_unique_seeds(self) -> None:
        with self.assertRaisesRegex(ValueError, "equal lengths"):
            five_seed_direction_stability([1.0] * 5, [0.0] * 4)
        with self.assertRaisesRegex(ValueError, "exactly five aligned seeds"):
            five_seed_direction_stability([1.0] * 4, [0.0] * 4)
        with self.assertRaisesRegex(ValueError, "five unique seeds"):
            five_seed_direction_stability(
                [1.0] * 5,
                [0.0] * 5,
                [1, 1, 2, 3, 4],
            )

    def test_resampling_defaults_are_locked(self) -> None:
        for helper in (two_way_donor_block_bootstrap, ld_block_paired_bootstrap):
            parameters = signature(helper).parameters
            self.assertEqual(parameters["n_resamples"].default, 10_000)
            self.assertEqual(parameters["seed"].default, DEFAULT_BOOTSTRAP_SEED)


class GainConcentrationTests(unittest.TestCase):
    """Where a pooled gain comes from, reported and never gated."""

    def test_below_min_strata_it_returns_not_applicable_and_no_numbers(self) -> None:
        result = gain_concentration([0.9, 0.8, 0.7], [0.8, 0.7, 0.6], ["a", "b", "c"])
        self.assertFalse(result.applicable)
        self.assertIn("not estimable", result.reason)
        self.assertEqual(result.leave_one_out_pooled_gain, ())
        self.assertIsNone(result.most_influential_stratum)
        self.assertIsNone(result.largest_absolute_shift)
        self.assertIsNone(result.sign_flips_when_dropped)

    def test_it_never_supports_a_pass_fail_verdict(self) -> None:
        result = gain_concentration(
            [0.9, 0.8, 0.7, 0.6, 0.5],
            [0.8, 0.7, 0.6, 0.5, 0.4],
            ["a", "b", "c", "d", "e"],
        )
        self.assertTrue(result.applicable)
        self.assertFalse(result.supports_pass_fail_verdict)

    def test_an_evenly_spread_gain_shows_little_leverage(self) -> None:
        result = gain_concentration(
            [0.60, 0.61, 0.59, 0.60, 0.61],
            [0.50, 0.51, 0.49, 0.50, 0.51],
            ["a", "b", "c", "d", "e"],
        )
        self.assertTrue(result.applicable)
        self.assertLess(result.largest_absolute_shift, 0.005)
        self.assertFalse(result.sign_flips_when_dropped)

    def test_it_reproduces_the_gse189600_collapse(self) -> None:
        """One stratum carrying a pooled margin is the signature it exists for.

        Established independently by the 2026-08-25 audit: dropping GSE189600
        moved elastic_net from a 0.066 margin to 0.0004.
        """

        strata = ["s1", "s2", "s3", "s4", "s5", "s6", "GSE189600"]
        # Six strata at parity, one carrying the entire pooled margin.
        candidate = [0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.9620]
        baseline = [0.4999, 0.4999, 0.4999, 0.4999, 0.4999, 0.4999, 0.5000]
        result = gain_concentration(candidate, baseline, strata)
        self.assertTrue(result.applicable)
        self.assertAlmostEqual(result.pooled_gain, 0.066, places=3)
        self.assertEqual(result.most_influential_stratum, "GSE189600")
        self.assertAlmostEqual(result.pooled_gain_without_it, 0.0001, places=4)
        self.assertGreater(result.largest_absolute_shift, 0.06)

    def test_a_sign_flip_is_reported(self) -> None:
        result = gain_concentration(
            [0.50, 0.50, 0.50, 0.50, 0.90],
            [0.52, 0.52, 0.52, 0.52, 0.50],
            ["a", "b", "c", "d", "carrier"],
        )
        self.assertTrue(result.applicable)
        self.assertGreater(result.pooled_gain, 0.0)
        self.assertLess(result.pooled_gain_without_it, 0.0)
        self.assertTrue(result.sign_flips_when_dropped)

    def test_weights_change_the_pooled_gain(self) -> None:
        candidate = [0.9, 0.5, 0.5, 0.5]
        baseline = [0.5, 0.5, 0.5, 0.5]
        strata = ["big", "b", "c", "d"]
        unweighted = gain_concentration(candidate, baseline, strata)
        weighted = gain_concentration(
            candidate, baseline, strata, weights=[100.0, 1.0, 1.0, 1.0]
        )
        self.assertGreater(weighted.pooled_gain, unweighted.pooled_gain)

    def test_repeated_strata_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            gain_concentration(
                [0.9, 0.8, 0.7, 0.6], [0.8, 0.7, 0.6, 0.5], ["a", "a", "b", "c"]
            )

    def test_misaligned_inputs_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            gain_concentration([0.9, 0.8, 0.7, 0.6], [0.8, 0.7], ["a", "b", "c", "d"])


if __name__ == "__main__":
    unittest.main()

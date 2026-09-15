from __future__ import annotations

import unittest

import numpy as np

from masld_bench.observed_multiome_specialist_gate import (
    GENOMIC_BLOCKS,
    LINEAGES,
    N_DONORS,
    N_UNIT_ROWS,
    SEEDS,
    SpecialistGateError,
    UnitAxis,
    compare_specialist,
    validate_candidate_metadata,
)


def _axis() -> UnitAxis:
    donors = tuple(
        f"donor-{donor:02d}"
        for block in GENOMIC_BLOCKS
        for donor in range(N_DONORS)
        for _lineage in LINEAGES
    )
    lineages = tuple(
        lineage
        for _block in GENOMIC_BLOCKS
        for _donor in range(N_DONORS)
        for lineage in LINEAGES
    )
    blocks = tuple(
        block
        for block in GENOMIC_BLOCKS
        for _donor in range(N_DONORS)
        for _lineage in LINEAGES
    )
    return UnitAxis(donors, lineages, blocks)


def _metadata(native: str = "peak_accessibility_rate") -> dict[str, object]:
    return {
        "schema_version": "masld-bench-observed-multiome-specialist-candidate-v1",
        "dataset_id": "gse296875",
        "stage": "development",
        "model_id": "peakvi",
        "native_output_family": native,
        "comparison_output_family": "donor_lineage_target_profile",
        "common_profile_projection": {
            "projection_id": "outer-fold-frozen-peak-to-target-map-v1",
            "source_sha256": "a" * 64,
            "fit_scope": "parameter_free",
            "frozen_before_scoring": True,
            "held_outcome_access": False,
            "sealed_outcome_access": False,
            "native_lane_preserved": True,
        },
        "fixed_seeds": list(SEEDS),
        "seed_ensemble": "arithmetic_mean_prediction_before_scoring",
        "prediction_committed_before_outcome": True,
        "partial_rectangle_used": False,
        "supports_external_claim": False,
        "supports_champion_claim": False,
        "sealed_outcomes_read": False,
    }


def _skills(gain: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    baseline = np.full(N_UNIT_ROWS, 0.2)
    candidate = baseline + gain * (1.0 - baseline)
    baseline_seeds = np.tile(baseline, (len(SEEDS), 1))
    candidate_seeds = np.tile(candidate, (len(SEEDS), 1))
    return candidate, candidate_seeds, baseline, baseline_seeds


class ObservedMultiomeSpecialistGateTests(unittest.TestCase):
    def test_qualifying_candidate_passes_all_prespecified_rules(self) -> None:
        candidate, candidate_seeds, baseline, baseline_seeds = _skills(0.06)
        receipt, arrays = compare_specialist(
            metadata=_metadata(),
            candidate_axis=_axis(),
            baseline_axis=_axis(),
            candidate_ensemble_skill=candidate,
            candidate_seed_skill=candidate_seeds,
            baseline_ensemble_skill=baseline,
            baseline_seed_skill=baseline_seeds,
            n_resamples=100,
            bootstrap_seed=7,
        )
        self.assertTrue(receipt["development_gate_passed"])
        self.assertEqual(receipt["improved_lineages"], 5)
        self.assertEqual(receipt["positive_gain_seeds"], 5)
        self.assertEqual(arrays["ensemble_relative_residual_deviance_gain"].shape, (975,))

    def test_direct_cross_family_entry_is_rejected(self) -> None:
        metadata = _metadata()
        metadata["comparison_output_family"] = "peak_accessibility_rate"
        with self.assertRaises(SpecialistGateError):
            validate_candidate_metadata(metadata)

    def test_projection_that_reads_held_outcomes_is_rejected(self) -> None:
        metadata = _metadata()
        metadata["common_profile_projection"]["held_outcome_access"] = True
        with self.assertRaises(SpecialistGateError):
            validate_candidate_metadata(metadata)

    def test_missing_seed_is_rejected(self) -> None:
        candidate, candidate_seeds, baseline, baseline_seeds = _skills(0.06)
        with self.assertRaises(SpecialistGateError):
            compare_specialist(
                metadata=_metadata(),
                candidate_axis=_axis(),
                baseline_axis=_axis(),
                candidate_ensemble_skill=candidate,
                candidate_seed_skill=candidate_seeds[:-1],
                baseline_ensemble_skill=baseline,
                baseline_seed_skill=baseline_seeds,
                n_resamples=10,
            )

    def test_axis_reordering_is_rejected(self) -> None:
        candidate, candidate_seeds, baseline, baseline_seeds = _skills(0.06)
        axis = _axis()
        changed = UnitAxis(axis.donor_hash[::-1], axis.lineage[::-1], axis.genomic_block[::-1])
        with self.assertRaises(SpecialistGateError):
            compare_specialist(
                metadata=_metadata(),
                candidate_axis=changed,
                baseline_axis=axis,
                candidate_ensemble_skill=candidate,
                candidate_seed_skill=candidate_seeds,
                baseline_ensemble_skill=baseline,
                baseline_seed_skill=baseline_seeds,
                n_resamples=10,
            )

    def test_four_lineage_rule_and_degradation_cap_are_both_enforced(self) -> None:
        candidate, candidate_seeds, baseline, baseline_seeds = _skills(0.07)
        axis = _axis()
        bad = np.asarray(axis.lineage) == "t_cell"
        candidate[bad] = baseline[bad] - 0.03 * (1.0 - baseline[bad])
        receipt, _ = compare_specialist(
            metadata=_metadata(),
            candidate_axis=axis,
            baseline_axis=axis,
            candidate_ensemble_skill=candidate,
            candidate_seed_skill=candidate_seeds,
            baseline_ensemble_skill=baseline,
            baseline_seed_skill=baseline_seeds,
            n_resamples=100,
            bootstrap_seed=7,
        )
        self.assertEqual(receipt["improved_lineages"], 4)
        self.assertFalse(receipt["development_gate_checks"]["no_lineage_worse_by_more_than_0_02"])
        self.assertFalse(receipt["development_gate_passed"])

    def test_gain_below_five_percent_fails(self) -> None:
        candidate, candidate_seeds, baseline, baseline_seeds = _skills(0.049)
        receipt, _ = compare_specialist(
            metadata=_metadata(),
            candidate_axis=_axis(),
            baseline_axis=_axis(),
            candidate_ensemble_skill=candidate,
            candidate_seed_skill=candidate_seeds,
            baseline_ensemble_skill=baseline,
            baseline_seed_skill=baseline_seeds,
            n_resamples=100,
            bootstrap_seed=7,
        )
        self.assertFalse(receipt["development_gate_checks"]["minimum_relative_deviance_gain"])
        self.assertFalse(receipt["development_gate_passed"])


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import unittest

import numpy as np

from masld_bench.observed_multiome_specialist_gate_v2 import (
    GENOMIC_BLOCKS,
    LINEAGES,
    MODEL_NATIVE_OUTPUT_FAMILY,
    N_DONORS,
    N_UNIT_ROWS,
    SEEDS,
    SpecialistGateError,
    UnitAxis,
    compare_specialist,
    validate_candidate_metadata,
)


def _axis() -> UnitAxis:
    return UnitAxis(
        tuple(f"donor-{donor:02d}" for block in GENOMIC_BLOCKS for donor in range(N_DONORS) for lineage in LINEAGES),
        tuple(lineage for block in GENOMIC_BLOCKS for donor in range(N_DONORS) for lineage in LINEAGES),
        tuple(block for block in GENOMIC_BLOCKS for donor in range(N_DONORS) for lineage in LINEAGES),
    )


def _metadata(model_id: str = "peakvi") -> dict[str, object]:
    return {
        "schema_version": "masld-bench-observed-multiome-specialist-candidate-v2",
        "dataset_id": "gse296875",
        "stage": "development",
        "model_id": model_id,
        "native_output_family": MODEL_NATIVE_OUTPUT_FAMILY[model_id],
        "comparison_output_family": "donor_lineage_target_profile",
        "common_profile_projection": {
            "projection_id": "outer-fold-frozen-common-profile-v2",
            "source_sha256": "a" * 64,
            "fit_scope": "outer_training_only",
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
    return candidate, np.tile(candidate, (5, 1)), baseline, np.tile(baseline, (5, 1))


class ObservedMultiomeSpecialistGateV2Tests(unittest.TestCase):
    def test_authoritative_native_families_are_distinct(self) -> None:
        self.assertEqual(MODEL_NATIVE_OUTPUT_FAMILY["multivi"], "joint_representation")
        self.assertEqual(MODEL_NATIVE_OUTPUT_FAMILY["peakvi"], "joint_representation")
        self.assertEqual(MODEL_NATIVE_OUTPUT_FAMILY["scooby_epicardioids"], "functional_track")
        self.assertEqual(MODEL_NATIVE_OUTPUT_FAMILY["scooby_neurips"], "functional_track")
        for model_id in MODEL_NATIVE_OUTPUT_FAMILY:
            validate_candidate_metadata(_metadata(model_id))

    def test_mislabeled_peakvi_and_generic_scooby_are_rejected(self) -> None:
        peakvi = _metadata("peakvi")
        peakvi["native_output_family"] = "peak_accessibility_rate"
        with self.assertRaises(SpecialistGateError):
            validate_candidate_metadata(peakvi)
        generic = _metadata("scooby_neurips")
        generic["model_id"] = "scooby"
        with self.assertRaises(SpecialistGateError):
            validate_candidate_metadata(generic)

    def test_qualifying_candidate_passes_exact_gate(self) -> None:
        candidate, candidate_seeds, baseline, baseline_seeds = _skills(0.06)
        receipt, arrays = compare_specialist(
            metadata=_metadata(), candidate_axis=_axis(), baseline_axis=_axis(),
            candidate_ensemble_skill=candidate, candidate_seed_skill=candidate_seeds,
            baseline_ensemble_skill=baseline, baseline_seed_skill=baseline_seeds,
            n_resamples=100, bootstrap_seed=7,
        )
        self.assertTrue(receipt["development_gate_passed"])
        self.assertEqual(receipt["improved_lineages"], 5)
        self.assertEqual(receipt["positive_gain_seeds"], 5)
        self.assertEqual(arrays["ensemble_relative_residual_deviance_gain"].shape, (975,))

    def test_three_improved_lineages_fails_four_of_five_rule(self) -> None:
        candidate, candidate_seeds, baseline, baseline_seeds = _skills(0.06)
        axis = _axis()
        lineages = np.asarray(axis.lineage)
        for lineage in ("macrophage", "t_cell"):
            query = lineages == lineage
            candidate[query] = baseline[query]
        receipt, _ = compare_specialist(
            metadata=_metadata(), candidate_axis=axis, baseline_axis=axis,
            candidate_ensemble_skill=candidate, candidate_seed_skill=candidate_seeds,
            baseline_ensemble_skill=baseline, baseline_seed_skill=baseline_seeds,
            n_resamples=100, bootstrap_seed=7,
        )
        self.assertEqual(receipt["improved_lineages"], 3)
        self.assertFalse(receipt["development_gate_checks"]["minimum_four_improved_lineages"])
        self.assertFalse(receipt["development_gate_passed"])

    def test_lineage_degradation_and_missing_seed_fail_closed(self) -> None:
        candidate, candidate_seeds, baseline, baseline_seeds = _skills(0.07)
        axis = _axis()
        query = np.asarray(axis.lineage) == "t_cell"
        candidate[query] = baseline[query] - 0.03 * (1.0 - baseline[query])
        receipt, _ = compare_specialist(
            metadata=_metadata(), candidate_axis=axis, baseline_axis=axis,
            candidate_ensemble_skill=candidate, candidate_seed_skill=candidate_seeds,
            baseline_ensemble_skill=baseline, baseline_seed_skill=baseline_seeds,
            n_resamples=100, bootstrap_seed=7,
        )
        self.assertFalse(receipt["development_gate_checks"]["no_lineage_worse_by_more_than_0_02"])
        with self.assertRaises(SpecialistGateError):
            compare_specialist(
                metadata=_metadata(), candidate_axis=axis, baseline_axis=axis,
                candidate_ensemble_skill=candidate, candidate_seed_skill=candidate_seeds[:-1],
                baseline_ensemble_skill=baseline, baseline_seed_skill=baseline_seeds,
                n_resamples=10,
            )

    def test_axis_order_and_direct_cross_family_entry_are_rejected(self) -> None:
        candidate, candidate_seeds, baseline, baseline_seeds = _skills(0.06)
        axis = _axis()
        changed = UnitAxis(axis.donor_hash[::-1], axis.lineage[::-1], axis.genomic_block[::-1])
        with self.assertRaises(SpecialistGateError):
            compare_specialist(
                metadata=_metadata(), candidate_axis=changed, baseline_axis=axis,
                candidate_ensemble_skill=candidate, candidate_seed_skill=candidate_seeds,
                baseline_ensemble_skill=baseline, baseline_seed_skill=baseline_seeds,
                n_resamples=10,
            )
        metadata = _metadata()
        metadata["comparison_output_family"] = "joint_representation"
        with self.assertRaises(SpecialistGateError):
            validate_candidate_metadata(metadata)


if __name__ == "__main__":
    unittest.main()

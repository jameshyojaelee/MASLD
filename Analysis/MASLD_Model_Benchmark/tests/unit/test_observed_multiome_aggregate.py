from __future__ import annotations

import unittest

import numpy as np

from masld_bench.adapters.observed_multiome_factorized import MODEL_IDS
from masld_bench.observed_multiome_aggregate import DONOR_FOLDS, GENOMIC_FOLDS, LINEAGES, SEEDS, ObservedMultiomeAggregateError, SurfaceRecord, aggregate_records, expected_surface_keys


def _records() -> list[SurfaceRecord]:
    records = []
    donors_by_fold = {0: range(0, 8), 1: range(8, 16), 2: range(16, 24), 3: range(24, 32), 4: range(32, 39)}
    for seed in SEEDS:
        for genomic in GENOMIC_FOLDS:
            for donor_fold in DONOR_FOLDS:
                units = tuple(f"donor-{donor:02d}" for donor in donors_by_fold[donor_fold] for _ in LINEAGES)
                lineages = tuple(lineage for _ in donors_by_fold[donor_fold] for lineage in LINEAGES)
                observed = np.tile(np.asarray([8.0, 2.0, 1.0]), (len(units), 1))
                predictions = {}
                for model_index, model_id in enumerate(MODEL_IDS):
                    base = np.asarray([8.0 - 0.6 * model_index, 2.0 + 0.3 * model_index, 1.0 + 0.3 * model_index])
                    predictions[model_id] = np.tile(base + (seed - SEEDS[0]) * 0.001, (len(units), 1))
                records.append(SurfaceRecord(seed, genomic, donor_fold, units, lineages, observed, predictions))
    return records


class ObservedMultiomeAggregateTests(unittest.TestCase):
    def test_lineage_roster_matches_the_frozen_production_identifier_authority(self) -> None:
        self.assertEqual(LINEAGES, ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell"))

    def test_expected_rectangle_contains_125_unique_surfaces(self) -> None:
        self.assertEqual(len(expected_surface_keys()), 125)

    def test_complete_rectangle_uses_ensemble_and_two_way_units(self) -> None:
        receipt, arrays = aggregate_records(_records(), n_resamples=100, bootstrap_seed=19)
        self.assertEqual(receipt["surface_seed_runs"], 125)
        self.assertEqual(receipt["biological_donors"], 39)
        self.assertEqual(receipt["genomic_blocks"], 5)
        self.assertEqual(receipt["unit_rows"], 975)
        self.assertEqual(receipt["strongest_development_control"], MODEL_IDS[0])
        self.assertFalse(receipt["seeds_are_biological_replicates"])
        self.assertEqual(arrays["ensemble_profile_deviance_skill"].shape, (975, 5))
        self.assertEqual(arrays["seed_profile_deviance_skill"].shape, (5, 975, 5))
        self.assertEqual(receipt["models"][MODEL_IDS[0]]["two_way_bootstrap_vs_uniform"]["n_donors"], 39)

    def test_partial_rectangle_is_rejected_before_values_are_aggregated(self) -> None:
        records = _records()
        with self.assertRaises(ObservedMultiomeAggregateError):
            aggregate_records(records[:-1], n_resamples=100)

    def test_seed_join_drift_is_rejected(self) -> None:
        records = _records()
        changed = list(records)
        target = changed[25]
        changed[25] = SurfaceRecord(target.seed, target.held_genomic_fold, target.held_donor_fold, ("different",) + target.unit_hash[1:], target.lineage, target.observed, target.predictions)
        with self.assertRaises(ObservedMultiomeAggregateError):
            aggregate_records(changed, n_resamples=100)


if __name__ == "__main__":
    unittest.main()

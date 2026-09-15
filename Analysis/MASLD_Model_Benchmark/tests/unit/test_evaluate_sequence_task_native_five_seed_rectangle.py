from __future__ import annotations

import unittest

from scripts.evaluate_sequence_task_native_five_seed_rectangle import (
    FiveSeedEvaluationError,
    FOLDS,
    LINEAGES,
    SEEDS,
    VIEWS,
    _macro_fisher_z,
    stratified_two_way_bootstrap,
    validate_candidate_roster,
)


class SequenceTaskNativeFiveSeedRectangleEvaluationTests(unittest.TestCase):
    def test_candidate_roster_requires_every_view_and_seed(self) -> None:
        rows = [
            {"candidate_id": f"{view}__s{seed}", "seed": str(seed)}
            for view in VIEWS
            for seed in SEEDS
        ]
        validate_candidate_roster(rows)
        with self.assertRaises(FiveSeedEvaluationError):
            validate_candidate_roster(rows[:-1])

    def test_bootstrap_reuses_fold_draws_across_lineages(self) -> None:
        rows = []
        for fold in FOLDS:
            for lineage_index, lineage in enumerate(LINEAGES):
                for donor in ("d0", "d1"):
                    for block in ("b0", "b1"):
                        rows.append(
                            {
                                "outer_fold": fold,
                                "lineage_id": lineage,
                                "donor_hash": f"f{fold}-{donor}",
                                "block_hash": f"f{fold}-{block}",
                                "skill": 0.1 + lineage_index * 0.01,
                            }
                        )
        result = stratified_two_way_bootstrap(rows, n_resamples=100, seed=17)
        self.assertAlmostEqual(result["estimate"], 0.12)
        self.assertEqual(result["n_resamples"], 100)
        self.assertTrue(result["lineage_macro_weighted"])

    def test_bootstrap_rejects_missing_donor_block_cell(self) -> None:
        rows = []
        for fold in FOLDS:
            for lineage in LINEAGES:
                for donor in ("d0", "d1"):
                    for block in ("b0", "b1"):
                        if fold == 0 and lineage == LINEAGES[0] and donor == "d1" and block == "b1":
                            continue
                        rows.append(
                            {
                                "outer_fold": fold,
                                "lineage_id": lineage,
                                "donor_hash": f"f{fold}-{donor}",
                                "block_hash": f"f{fold}-{block}",
                                "skill": 0.1,
                            }
                        )
        with self.assertRaises(FiveSeedEvaluationError):
            stratified_two_way_bootstrap(rows, n_resamples=100, seed=17)

    def test_correlation_macro_uses_fisher_z_with_equal_stratum_weight(self) -> None:
        rows = [
            {
                "outer_fold": fold,
                "lineage_id": lineage,
                "correlation": 0.5,
            }
            for fold in FOLDS
            for lineage in LINEAGES
        ]
        self.assertAlmostEqual(_macro_fisher_z(rows, "correlation"), 0.5)


if __name__ == "__main__":
    unittest.main()

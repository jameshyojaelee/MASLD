from __future__ import annotations

from pathlib import Path
import unittest

from scripts.prepare_chrombpnet_bundle_partial_evaluation import (
    PartialEvaluationError,
    candidate_row,
    validate_expected_units,
)


class ChromBPNetBundlePartialEvaluationTests(unittest.TestCase):
    def test_unit_roster_preserves_actual_seed_and_diagonal_split(self) -> None:
        units = [
            {
                "logical_work_id": f"seq__chrombpnet__lineage_{fold}__d{fold}g{fold}__s{100 + fold}",
                "lineage_id": f"lineage_{fold}",
                "outer_fold": fold,
                "split_id": f"donor{fold}_genomic{fold}",
                "seed": 100 + fold,
            }
            for fold in range(5)
        ]
        units.extend(
            {
                "logical_work_id": f"seq__chrombpnet__extra_{index}__d0g0__s{200 + index}",
                "lineage_id": f"extra_{index}",
                "outer_fold": 0,
                "split_id": "donor0_genomic0",
                "seed": 200 + index,
            }
            for index in range(3)
        )
        validated = validate_expected_units(units)
        self.assertEqual([row["seed"] for row in validated], [100, 101, 102, 103, 104, 200, 201, 202])

    def test_non_diagonal_split_is_rejected(self) -> None:
        unit = {
            "logical_work_id": "seq__chrombpnet__x__d0g1__s1",
            "lineage_id": "x",
            "outer_fold": 0,
            "split_id": "donor0_genomic1",
            "seed": 1,
        }
        with self.assertRaises(PartialEvaluationError):
            validate_expected_units([unit] * 8)

    def test_native_views_are_singleton_candidates_with_actual_seed(self) -> None:
        root = Path("/tmp/frozen-prediction")
        full = candidate_row(
            view_id="chrombpnet_full",
            prediction_subdir="predictions/full_model",
            root=root,
            artifacts_sha256="a" * 64,
            seed=20260825,
        )
        nobias = candidate_row(
            view_id="chrombpnet_nobias",
            prediction_subdir="predictions/nobias_model",
            root=root,
            artifacts_sha256="a" * 64,
            seed=20260825,
        )
        self.assertEqual(full["candidate_id"], "chrombpnet_full__s20260825")
        self.assertEqual(nobias["candidate_id"], "chrombpnet_nobias__s20260825")
        self.assertEqual(full["seed"], "20260825")
        self.assertNotEqual(full["prediction_subdir"], nobias["prediction_subdir"])


if __name__ == "__main__":
    unittest.main()

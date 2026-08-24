from __future__ import annotations

import unittest

from scripts.build_sequence_training_mean_baseline import (
    TrainingMeanBaselineError,
    role_windows,
    select_training_row,
)


class SequenceTrainingMeanBaselineTests(unittest.TestCase):
    def test_selects_exact_training_folds(self) -> None:
        split = {
            "donor_test_fold": "0",
            "donor_valid_fold": "1",
            "donor_train_folds": "2,3,4",
        }
        rows = [
            {
                "donor_test_fold": "0",
                "donor_valid_fold": "1",
                "donor_train_folds": "2,3,4",
                "lineage_id": "hepatocyte",
                "donors": "19",
            }
        ]
        self.assertEqual(
            select_training_row(rows, split=split, lineage="hepatocyte")["donors"],
            "19",
        )
        with self.assertRaises(TrainingMeanBaselineError):
            select_training_row(rows, split=split, lineage="macrophage")

    def test_role_windows_rejects_nonfrozen_census(self) -> None:
        split = {"genomic_test_fold": "0", "genomic_valid_fold": "1"}
        row = {
            "genomic_fold": "0",
            "output_start": "100",
            "output_end": "1100",
        }
        with self.assertRaises(TrainingMeanBaselineError):
            role_windows([row], split=split, role="test")


if __name__ == "__main__":
    unittest.main()

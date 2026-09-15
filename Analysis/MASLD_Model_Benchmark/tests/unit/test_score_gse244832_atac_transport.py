from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.score_gse244832_atac_transport import (
    ATACTransportScoringError,
    deviance_per_insertion,
    score_arrays,
    score_committed,
)


class ScoreGSE244832ATACTransportTests(unittest.TestCase):
    def test_exact_prediction_has_lower_deviance_than_uniform(self) -> None:
        truth = np.asarray([1, 3, 8, 2], dtype=float)
        self.assertLess(
            deviance_per_insertion(truth, truth),
            deviance_per_insertion(truth, np.ones(4)),
        )

    def test_count_family_scores_eighteen_donor_lineage_blocks(self) -> None:
        truth = np.ones((18, 4, 16000), dtype=np.uint32)
        score = truth.astype(np.float32)
        missing = np.zeros((18, 4, 16000), dtype=np.uint8)
        rows = score_arrays(
            truth,
            score,
            missing,
            block_ids=["chr1"] * 8000 + ["chr2"] * 8000,
            output_family="masked_accessibility_count",
        )
        self.assertEqual(len(rows), 18 * 4 * 2)
        self.assertTrue(all(abs(row["deviance_per_insertion"]) < 1e-12 for row in rows))

    def test_profile_family_preserves_twenty_bin_axis(self) -> None:
        truth = np.ones((18, 4, 16000, 20), dtype=np.uint32)
        score = truth.astype(np.float32)
        missing = np.zeros((18, 4, 16000), dtype=np.uint8)
        rows = score_arrays(
            truth,
            score,
            missing,
            block_ids=["chr1"] * 16000,
            output_family="functional_track_profile",
        )
        self.assertEqual(len(rows), 72)

    def test_invalid_commit_stops_before_outcome_loader(self) -> None:
        called = False

        def loader(*_args: object) -> tuple[object, object, list[str]]:
            nonlocal called
            called = True
            raise AssertionError("outcome loader must not run")

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaises((FileNotFoundError, ATACTransportScoringError)):
                score_committed(
                    commit_root=root / "absent-commit",
                    commit_artifacts_sha256="0" * 64,
                    prediction_root=root / "prediction",
                    outcome_root=root / "outcome",
                    outcome_artifacts_sha256="1" * 64,
                    axis_root=root / "axis",
                    axis_artifacts_sha256="2" * 64,
                    registration_root=root / "registration",
                    registration_artifacts_sha256="3" * 64,
                    output=root / "score",
                    outcome_loader=loader,
                )
        self.assertFalse(called)


if __name__ == "__main__":
    unittest.main()

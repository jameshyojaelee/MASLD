from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from scripts.corgi_predict_outcome_aligned_tiles import (
    ATAC_CHANNEL,
    CorgiPredictionError,
    CHANNEL_RC,
    overlap_weights,
    read_valid_lane_contract,
    regional_mean,
    regional_orientation_mean_softplus_sum,
    regional_softplus_sum,
)


class CorgiOutcomeAlignedPredictionTests(unittest.TestCase):
    def test_atac_channel_is_reverse_complement_invariant(self) -> None:
        self.assertEqual(int(CHANNEL_RC[ATAC_CHANNEL]), ATAC_CHANNEL)

    def test_mapper_receipt_defines_synchronized_valid_lane(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            mapper = Path(value)
            receipt = {
                "schema_version": "masld-bench-corgi-masked-context-mapper-v1",
                "status": "pass_outcome_free_mapper",
                "outer_fold": 4,
                "valid_fold": 0,
                "seed": 20260824,
                "held_ATAC_or_other_outcomes_used": False,
                "test_or_sealed_outcomes_read": False,
            }
            (mapper / "receipt.json").write_text(json.dumps(receipt), encoding="utf-8")
            observed, donor_fold, genomic_fold = read_valid_lane_contract(mapper, 4)
            self.assertEqual(observed, receipt)
            self.assertEqual(donor_fold, 0)
            self.assertEqual(genomic_fold, 0)

    def test_mapper_outer_fold_mismatch_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as value:
            mapper = Path(value)
            receipt = {
                "schema_version": "masld-bench-corgi-masked-context-mapper-v1",
                "status": "pass_outcome_free_mapper",
                "outer_fold": 0,
                "valid_fold": 1,
                "seed": 20260824,
                "held_ATAC_or_other_outcomes_used": False,
                "test_or_sealed_outcomes_read": False,
            }
            (mapper / "receipt.json").write_text(json.dumps(receipt), encoding="utf-8")
            with self.assertRaises(CorgiPredictionError):
                read_valid_lane_contract(mapper, 1)

    def test_partial_bin_overlap_is_exact(self) -> None:
        bins, overlap = overlap_weights(10, 110)
        self.assertEqual(bins.tolist(), [0, 1])
        self.assertEqual(overlap.tolist(), [54.0, 46.0])

    def test_regional_reductions(self) -> None:
        track = np.zeros((2, 6_144), dtype=np.float32)
        track[0, 0:2] = 2.0
        track[1, 0:2] = 4.0
        mean = regional_mean(track, 10, 110)
        self.assertTrue(np.allclose(mean, [2.0, 4.0]))
        expected = np.logaddexp(0.0, np.asarray([2.0, 4.0])) * (100.0 / 64.0)
        self.assertTrue(np.allclose(regional_softplus_sum(track, 10, 110), expected))

    def test_softplus_precedes_orientation_average(self) -> None:
        forward = np.full((1, 6_144), -2.0, dtype=np.float32)
        reverse = np.full((1, 6_144), 2.0, dtype=np.float32)
        observed = regional_orientation_mean_softplus_sum(forward, reverse, 10, 110)
        expected = (
            np.logaddexp(0.0, -2.0) + np.logaddexp(0.0, 2.0)
        ) / 2.0 * (100.0 / 64.0)
        self.assertTrue(np.allclose(observed, [expected]))
        jensen_underestimate = np.logaddexp(0.0, 0.0) * (100.0 / 64.0)
        self.assertGreater(float(observed[0]), float(jensen_underestimate))


if __name__ == "__main__":
    unittest.main()

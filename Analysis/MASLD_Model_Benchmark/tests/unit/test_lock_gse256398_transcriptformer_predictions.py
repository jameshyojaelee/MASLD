from __future__ import annotations

from hashlib import sha256
import unittest

import numpy as np

from scripts.lock_gse256398_transcriptformer_predictions import (
    arm_mask,
    leave_s35_masks_identical,
    update_prediction_digest,
)


class GSE256398TranscriptFormerPredictionLockTests(unittest.TestCase):
    def test_leave_s35_out_masks_both_arms(self) -> None:
        source = np.asarray([True, False], dtype=bool)
        self.assertFalse(arm_mask(arm="prospective_capped_leave_s35_out", donor="S35", source_mask=source, rows=2).any())
        self.assertFalse(arm_mask(arm="source_exact_leave_s35_out", donor="S35", source_mask=source, rows=2).any())

    def test_prediction_digest_depends_on_row_identity(self) -> None:
        embedding = np.ones((1, 2), dtype=np.float32)
        first = sha256()
        second = sha256()
        update_prediction_digest(first, ["row-a"], embedding)
        update_prediction_digest(second, ["row-b"], embedding)
        self.assertNotEqual(first.hexdigest(), second.hexdigest())

    def test_leave_s35_stream_equality_is_checked_before_namespaced_hashing(self) -> None:
        source = np.asarray([True, False], dtype=bool)
        self.assertTrue(
            leave_s35_masks_identical(donor="S35", source_mask=source, rows=2)
        )
        self.assertFalse(
            leave_s35_masks_identical(donor="S1", source_mask=source, rows=2)
        )
        self.assertTrue(
            leave_s35_masks_identical(
                donor="S1", source_mask=np.ones(2, dtype=bool), rows=2
            )
        )


if __name__ == "__main__":
    unittest.main()

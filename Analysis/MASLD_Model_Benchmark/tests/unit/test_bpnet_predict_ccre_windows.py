from __future__ import annotations

import unittest

import numpy as np

from scripts import bpnet_predict_ccre_windows as prediction


class BPNetPredictionTests(unittest.TestCase):
    def test_log1p_count_inverse_is_nonnegative(self) -> None:
        observed = prediction.count_mass(np.asarray([0.0, np.log(2.0), -2.0]))
        np.testing.assert_allclose(observed, [0.0, 1.0, 0.0])

    def test_nonfinite_or_extreme_count_fails(self) -> None:
        for values in (np.asarray([np.nan]), np.asarray([51.0])):
            with self.assertRaises(prediction.BPNetContractError):
                prediction.count_mass(values)

    def test_batch_partition_is_complete(self) -> None:
        self.assertEqual(
            list(prediction.batches(list(range(7)), 3)),
            [[0, 1, 2], [3, 4, 5], [6]],
        )


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import numpy as np
import unittest

from scripts.sei_predict_gse281364 import SeiPredictionError, scalarize


class SeiGSE281364PredictionTests(unittest.TestCase):
    def test_scalarization_is_outcome_free_and_tie_stable(self) -> None:
        scores = np.zeros((2, 40), dtype=np.float32)
        scores[0, 3] = -2.0
        scores[0, 7] = 2.0
        scores[1, 5] = 1.5
        index, signed, maximum, mean = scalarize(scores)
        np.testing.assert_array_equal(index, [3, 5])
        np.testing.assert_allclose(signed, [-2.0, 1.5])
        np.testing.assert_allclose(maximum, [2.0, 1.5])
        np.testing.assert_allclose(mean, [0.0, 1.5 / 40])

    def test_rejects_wrong_native_axis(self) -> None:
        with self.assertRaises(SeiPredictionError):
            scalarize(np.zeros((2, 39), dtype=np.float32))


if __name__ == "__main__":
    unittest.main()

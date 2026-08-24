from __future__ import annotations

import unittest

import numpy as np

from scripts.corgi_predict_valid_contexts import (
    CHANNEL_RC,
    one_hot,
    restore_reverse_complement,
    reverse_complement_one_hot,
)


class CorgiValidPredictionTests(unittest.TestCase):
    def test_reverse_complement_one_hot(self) -> None:
        sequence = ("ACGT" * 131_072)
        encoded = one_hot(sequence)
        self.assertTrue(np.array_equal(reverse_complement_one_hot(encoded), encoded))

    def test_stranded_channel_and_bin_restoration(self) -> None:
        values = np.arange(1 * 1 * 22 * 6_144, dtype=np.float32).reshape(1, 1, 22, 6_144)
        restored = restore_reverse_complement(values)
        self.assertEqual(restored[0, 0, 12, 0], values[0, 0, 13, -1])
        self.assertEqual(CHANNEL_RC[16], 17)


if __name__ == "__main__":
    unittest.main()

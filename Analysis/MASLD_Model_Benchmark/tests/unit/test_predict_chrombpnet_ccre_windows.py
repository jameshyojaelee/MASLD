from __future__ import annotations

import unittest

import numpy as np

from scripts import predict_chrombpnet_ccre_windows as prediction


class ChromBPNetCCREPredictionTests(unittest.TestCase):
    def test_reverse_complement_one_hot_is_an_involution(self) -> None:
        values = prediction.one_hot_dna(["ACGT", "TTAG"])
        restored = prediction.reverse_complement_one_hot(
            prediction.reverse_complement_one_hot(values)
        )
        np.testing.assert_array_equal(restored, values)

    def test_profile_and_count_transforms_are_finite(self) -> None:
        probabilities = prediction.softmax(
            np.asarray([[0.0, 1.0, -1.0], [1000.0, 999.0, 998.0]])
        )
        np.testing.assert_allclose(
            probabilities.sum(axis=1), np.ones(2), rtol=0, atol=1e-7
        )
        np.testing.assert_allclose(
            prediction.regional_mass_from_logscore(
                np.asarray([0.0, np.log(3.0)]), "log1p_absolute"
            ),
            np.asarray([0.0, 2.0]),
            rtol=0,
            atol=1e-12,
        )
        np.testing.assert_allclose(
            prediction.regional_mass_from_logscore(
                np.asarray([0.0, np.log(3.0)]), "log_component"
            ),
            np.asarray([1.0, 3.0]),
            rtol=0,
            atol=1e-12,
        )

    def test_ambiguous_dna_is_rejected(self) -> None:
        with self.assertRaises(prediction.ChromBPNetPredictionError):
            prediction.one_hot_dna(["ACNT"])


if __name__ == "__main__":
    unittest.main()

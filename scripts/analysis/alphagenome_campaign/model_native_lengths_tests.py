#!/usr/bin/env python3
"""Scientific invariants for context-length native calibration."""
import unittest
import numpy as np
import pandas as pd
from model_native_lengths import calibrate_scores


class CalibrationTests(unittest.TestCase):
    def setUp(self):
        self.train = pd.DataFrame(dict(key=["a", "b", "c"], heldout_fold=[2, 3, 4],
                                      beta_alt=[2., -4., 6.], score=[1., -2., 3.]))
        self.valid = pd.DataFrame(dict(key=["d", "e", "f"], heldout_fold=[1, 1, 1],
                                      beta_alt=[100., 50., -20.], score=[0., 2., -2.]))

    def test_training_only_and_allele_symmetry(self):
        first, receipt = calibrate_scores(self.train, self.valid, ["score"])
        changed = self.valid.assign(beta_alt=[-800., 900., 1234.])
        second, _ = calibrate_scores(self.train, changed, ["score"])
        np.testing.assert_array_equal(first.score_calibrated, second.score_calibrated)
        np.testing.assert_array_equal(first.score_calibrated, [0., 4., -4.])
        self.assertEqual(receipt[0]["slope"], 2.)

    def test_forbidden_folds_and_duplicate_identity(self):
        with self.assertRaises(ValueError):
            calibrate_scores(self.train.assign(heldout_fold=[0, 3, 4]), self.valid, ["score"])
        with self.assertRaises(ValueError):
            calibrate_scores(self.train, self.valid.assign(key=["a", "e", "f"]), ["score"])


if __name__ == "__main__":
    unittest.main()

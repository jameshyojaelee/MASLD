import itertools
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from h48_stats import empirical_two_sided, fraction_below, read_null_array


class EmpiricalPTests(unittest.TestCase):
    def test_complete_rank_permutation_oracle(self):
        scores = np.array([1.0, 3.0, 2.0, 4.0])
        weights = np.arange(4.0)
        observed = spearmanr(scores, weights).statistic
        null = np.array([spearmanr(p, weights).statistic for p in itertools.permutations(scores)])
        expected = (1 + sum(abs(v) >= abs(observed) for v in null)) / 25
        self.assertEqual(empirical_two_sided(observed, null), expected)
        self.assertEqual(empirical_two_sided(-observed, null), expected)

    def test_ties_and_actual_number_of_draws(self):
        self.assertEqual(empirical_two_sided(2, [-2, 0, 2]), 3 / 4)
        self.assertEqual(empirical_two_sided(3, [-2, 0, 2]), 1 / 4)
        self.assertEqual(empirical_two_sided(0, [-2, 0, 2]), 1)

    def test_undefined_observed_or_draw_cannot_be_significant(self):
        for value in (np.nan, np.inf, -np.inf):
            self.assertTrue(np.isnan(empirical_two_sided(value, [0, 1, 2])))
            self.assertTrue(np.isnan(empirical_two_sided(3, [0, value, 2])))

    def test_invalid_null_shape(self):
        for values in ([], [[1, 2]]):
            with self.assertRaises(ValueError):
                empirical_two_sided(1, values)

    def test_missing_tests_do_not_count_as_nonsignificant_calibration(self):
        self.assertEqual(fraction_below([0.01, 0.2], 0.05), 0.5)
        for values in ([], [np.nan], [0.01, np.nan]):
            self.assertTrue(np.isnan(fraction_below(values, 0.05)))


class NullArrayTests(unittest.TestCase):
    def test_values_order_dtype_and_readonly(self):
        values = np.arange(30, dtype=np.float32).reshape(2, 3, 5) / 7
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'null.bin'
            values.tofile(path)
            mapped = read_null_array(path, values.shape)
            self.assertIsInstance(mapped, np.memmap)
            self.assertFalse(mapped.flags.writeable)
            np.testing.assert_array_equal(mapped, values)
            np.testing.assert_array_equal(mapped[1, :, [4, 0]], values[1, :, [4, 0]])
            for index in range(2):
                block = read_null_array(path, values.shape, index=index)
                np.testing.assert_array_equal(block, values[index])
                self.assertEqual(block.nbytes, values[index].nbytes)
            for index in (-1, 2):
                with self.assertRaises(IndexError):
                    read_null_array(path, values.shape, index=index)
            with self.assertRaises(ValueError):
                mapped[0, 0, 0] = 99

    def test_truncated_or_extra_input_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'null.bin'
            for size in (29, 31):
                np.arange(size, dtype=np.float32).tofile(path)
                with self.assertRaises(ValueError):
                    read_null_array(path, (2, 3, 5))


if __name__ == '__main__':
    unittest.main()

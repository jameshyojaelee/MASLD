from __future__ import annotations

import unittest

import numpy as np

from scripts.fit_corgi_masked_context_mapper import derangement, masked_quantile_map


class CorgiMaskedMapperTests(unittest.TestCase):
    def test_masked_quantile_uses_neutral_nonzero_imputation(self) -> None:
        values = np.asarray([0.0, 1.0, 2.0, 0.0])
        observed = np.asarray([True, True, True, False])
        reference = np.asarray([-1.0, 0.0, 1.0, 2.0])
        mapped = masked_quantile_map(values, observed, reference)
        self.assertEqual(mapped[-1], np.median(reference))
        self.assertTrue(np.all(np.diff(mapped[:3]) > 0))

    def test_derangement_is_complete(self) -> None:
        values = np.arange(8)
        observed = derangement(values, np.random.default_rng(17))
        self.assertEqual(set(observed), set(values))
        self.assertTrue(np.all(observed != values))


if __name__ == "__main__":
    unittest.main()

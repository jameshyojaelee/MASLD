from __future__ import annotations

import unittest

import numpy as np

from masld_cl.raw_pca import log_cpm_pca


class TestRawPCA(unittest.TestCase):
    def test_pca_is_donor_level_and_deterministic(self):
        counts = np.array([[10, 1, 0], [8, 2, 1], [0, 2, 12], [1, 1, 9]])
        left = log_cpm_pca(counts, 3)
        right = log_cpm_pca(counts, 3)
        self.assertEqual(left.shape, (4, 3))
        np.testing.assert_allclose(np.abs(left), np.abs(right), rtol=0, atol=1e-12)

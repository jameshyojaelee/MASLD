from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from masld_cl.metrics import donor_distance_spearman, standardized_case_control_separation
from masld_cl.orthogonal_bridge import (
    apply_dataset_control_offsets,
    apply_scaled_orthogonal_bridge,
    fit_scaled_orthogonal_bridge,
)


class TestOrthogonalBridge(unittest.TestCase):
    def test_bridge_preserves_distance_ranks_and_standardized_separation(self):
        rng = np.random.default_rng(17)
        source = rng.normal(size=(40, 4))
        q, _ = np.linalg.qr(rng.normal(size=(4, 4)))
        target = source @ q * 1.7 + 2.0
        sm, tm, rotation, scale = fit_scaled_orthogonal_bridge(source, target)
        mapped = apply_scaled_orthogonal_bridge(source, sm, tm, rotation, scale)
        self.assertAlmostEqual(donor_distance_spearman(source, mapped), 1.0)
        control = np.arange(len(source)) < 10
        left = standardized_case_control_separation(source[~control], source[control])
        right = standardized_case_control_separation(mapped[~control], mapped[control])
        self.assertAlmostEqual(left / right, 1.0, places=6)

    def test_dataset_offset_is_uniform_and_reference_is_exact(self):
        latent = np.asarray([[0.0, 0.0], [2.0, 0.0], [5.0, 1.0], [7.0, 1.0], [9.0, 1.0]])
        cells = pd.DataFrame({
            "donor_id": ["r1", "r2", "q1", "q2", "q3"],
            "dataset": ["R", "R", "Q", "Q", "Q"],
            "preparation": ["u"] * 5,
            "strict_reference": [True, True, False, False, False],
            "primary_query": [False, False, True, True, True],
            "query_control": [False, False, True, True, False],
        })
        out, rows = apply_dataset_control_offsets(latent, cells, ["Q"], 1.0, 2)
        np.testing.assert_array_equal(out[:2], latent[:2])
        np.testing.assert_allclose(out[3] - out[2], latent[3] - latent[2])
        np.testing.assert_allclose(out[4] - out[2], latent[4] - latent[2])
        self.assertLessEqual(rows[0]["maximum_centered_coordinate_error"], 1e-12)

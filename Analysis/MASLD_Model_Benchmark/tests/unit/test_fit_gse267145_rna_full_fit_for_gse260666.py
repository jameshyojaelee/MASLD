from __future__ import annotations

import numpy as np
import unittest

from scripts.fit_gse267145_rna_full_fit_for_gse260666 import (
    log1p_cpm,
    select_feature_indices,
)


class GSE267145RNAFullFitTests(unittest.TestCase):
    def test_log1p_cpm_is_per_participant(self) -> None:
        matrix = np.asarray([[1.0, 1.0], [10.0, 0.0]])
        transformed = log1p_cpm(matrix)
        self.assertAlmostEqual(transformed[0, 0], np.log1p(500_000.0))
        self.assertAlmostEqual(transformed[1, 0], np.log1p(1_000_000.0))
        self.assertEqual(transformed[1, 1], 0.0)

    def test_source_variance_then_stable_id_selects_features(self) -> None:
        matrix = np.asarray(
            [
                [0.0, 0.0, 1.0, 4.0],
                [1.0, 1.0, 1.0, 1.0],
                [2.0, 2.0, 1.0, 4.0],
                [3.0, 3.0, 1.0, 1.0],
            ]
        )
        selected = select_feature_indices(
            matrix,
            ["ENSG_B", "ENSG_A", "ENSG_C", "ENSG_D"],
            2,
        )
        self.assertEqual(selected.tolist(), [3, 1])


if __name__ == "__main__":
    unittest.main()

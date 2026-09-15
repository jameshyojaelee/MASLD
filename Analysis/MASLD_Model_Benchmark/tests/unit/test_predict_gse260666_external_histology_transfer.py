from __future__ import annotations

import numpy as np
import unittest

from scripts.predict_gse260666_external_histology_transfer import infer_seed, softmax


class GSE260666ExternalPredictionTests(unittest.TestCase):
    def test_softmax_is_stable_and_row_normalized(self) -> None:
        values = softmax(np.asarray([[1000.0, 999.0, 998.0], [-1000.0, -999.0, -998.0]]))
        np.testing.assert_allclose(values.sum(axis=1), 1.0)
        self.assertTrue(np.all(values > 0))

    def test_infer_seed_emits_three_valid_probability_models(self) -> None:
        rng = np.random.default_rng(7)
        state = {
            "center": np.zeros(4),
            "scale": np.ones(4),
            "pca_mean": np.zeros(4),
            "pca_components": np.eye(2, 4),
            "svm_coef": rng.normal(size=(3, 2)),
            "svm_intercept": rng.normal(size=3),
            "svm_calibrator_coef": rng.normal(size=(3, 3)),
            "svm_calibrator_intercept": rng.normal(size=3),
            "elastic_coef": rng.normal(size=(3, 2)),
            "elastic_intercept": rng.normal(size=3),
            "centroid_centroids": rng.normal(size=(3, 2)),
        }
        result = infer_seed(rng.normal(size=(5, 4)), state)
        self.assertEqual(
            set(result),
            {
                "rna_hvg_pca_elastic_net",
                "rna_hvg_pca_linear_svm",
                "rna_hvg_pca_nearest_centroid",
            },
        )
        for probabilities in result.values():
            self.assertEqual(probabilities.shape, (5, 3))
            np.testing.assert_allclose(probabilities.sum(axis=1), 1.0)


if __name__ == "__main__":
    unittest.main()

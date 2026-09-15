from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy import sparse

from masld_bench.observed_multiome_surface import (
    SurfacePreprocessState,
    apply_nonnegative_reduction,
    apply_sparse_reduction,
    apply_target_reduction,
    export_preprocess,
    fit_nonnegative_reduction,
    fit_sparse_reduction,
    fit_target_reduction,
    load_preprocess,
    log_cpm,
    profile_deviance_skill,
)


class ObservedMultiomeSurfaceTests(unittest.TestCase):
    def test_fold_local_reductions_round_trip_without_pickle(self) -> None:
        generator = np.random.default_rng(5)
        counts = sparse.csr_matrix(generator.poisson(2, size=(12, 9)))
        normalized, zero = log_cpm(counts)
        train, basis = fit_sparse_reduction(normalized, components=3, seed=7, iterations=5)
        np.testing.assert_allclose(train, apply_sparse_reduction(normalized, basis), rtol=1e-10, atol=1e-10)
        nonnegative, nonnegative_basis = fit_nonnegative_reduction(normalized, components=3, seed=7, max_iter=200)
        np.testing.assert_allclose(nonnegative, apply_nonnegative_reduction(normalized, nonnegative_basis), rtol=1e-10, atol=1e-10)
        self.assertTrue(np.all(nonnegative >= 0))
        target = generator.normal(size=(20, 6))
        embedded, mean, scale, target_basis = fit_target_reduction(target, components=3)
        np.testing.assert_allclose(embedded, apply_target_reduction(target, mean, scale, target_basis), rtol=1e-10, atol=1e-10)
        self.assertEqual(zero, 0)
        state = SurfacePreprocessState(basis, nonnegative_basis, mean, scale, target_basis, 3, 7)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "state"
            export_preprocess(state, path, {"training_rows": 12})
            restored, metadata = load_preprocess(path)
            np.testing.assert_array_equal(restored.rna_components, basis)
            self.assertEqual(metadata["training_rows"], 12)

    def test_zero_library_is_observed_zero_not_missing(self) -> None:
        matrix = sparse.csr_matrix([[0, 0, 0], [1, 0, 2]], dtype=float)
        transformed, zero = log_cpm(matrix)
        self.assertEqual(zero, 1)
        np.testing.assert_array_equal(transformed.toarray()[0], np.zeros(3))

    def test_profile_deviance_rewards_exact_profile(self) -> None:
        observed = np.asarray([[5, 3, 2], [0, 4, 6]], dtype=float)
        exact, _, null = profile_deviance_skill(observed, observed, pseudocount=1e-8)
        uniform, _, uniform_null = profile_deviance_skill(observed, np.ones_like(observed), pseudocount=1e-8)
        self.assertTrue(np.all(exact > 0.999999))
        np.testing.assert_allclose(uniform, np.zeros(2), atol=1e-12)
        np.testing.assert_array_equal(null, uniform_null)


if __name__ == "__main__":
    unittest.main()

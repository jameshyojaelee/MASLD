from __future__ import annotations

import unittest

import numpy as np

from scripts.fit_gse135251_nas_sum_transfer_source_models import (
    NasSourceFitError,
    apply_pipeline,
    fit_pipeline,
    log2_cpm_complete_axis,
    select_feature_indices,
    select_hyperparameters,
    spearman,
    spearman_null,
)


class ScoringTests(unittest.TestCase):
    def test_spearman_matches_scipy(self) -> None:
        from scipy.stats import spearmanr

        rng = np.random.default_rng(3)
        a = rng.normal(size=200)
        b = a * 0.6 + rng.normal(size=200)
        self.assertAlmostEqual(spearman(a, b), float(spearmanr(a, b).statistic), places=10)

    def test_spearman_is_invariant_to_a_monotone_shift(self) -> None:
        """The severity gap between the two cohorts must not move the primary metric."""

        rng = np.random.default_rng(11)
        truth = rng.integers(0, 9, size=99).astype(float)
        pred = truth + rng.normal(scale=1.0, size=99)
        shifted = pred + 1.688
        self.assertAlmostEqual(spearman(truth, pred), spearman(truth, shifted), places=12)

    def test_spearman_is_invariant_to_positive_rescaling(self) -> None:
        rng = np.random.default_rng(13)
        truth = rng.integers(0, 9, size=99).astype(float)
        pred = rng.normal(size=99)
        self.assertAlmostEqual(spearman(truth, pred), spearman(truth, pred * 3.7 + 2.0), places=12)

    def test_null_mean_is_near_zero_and_p95_near_the_analytic_value(self) -> None:
        rng = np.random.default_rng(17)
        truth = rng.integers(0, 9, size=99).astype(float)
        scores = rng.normal(size=99)
        null = spearman_null(truth, scores, replicates=4000, seed=5)
        self.assertAlmostEqual(null["null_mean"], 0.0, delta=0.02)
        # analytic reference 1.645/sqrt(98) = 0.1662
        self.assertAlmostEqual(null["null_p95"], 0.1662, delta=0.03)

    def test_outcome_ties_do_not_narrow_the_null(self) -> None:
        """24 tied zeros in the target must not be used to argue for a tighter null."""

        rng = np.random.default_rng(23)
        scores = rng.normal(size=99)
        untied = rng.normal(size=99)
        tied = np.concatenate([np.zeros(24), rng.integers(1, 9, size=75).astype(float)])
        a = spearman_null(untied, scores, replicates=4000, seed=7)
        b = spearman_null(tied, scores, replicates=4000, seed=7)
        self.assertAlmostEqual(a["null_p95"], b["null_p95"], delta=0.03)


class TransformTests(unittest.TestCase):
    def test_log2_cpm_removes_depth(self) -> None:
        counts = np.asarray([[1.0, 2.0, 3.0], [10.0, 20.0, 30.0]])
        out = log2_cpm_complete_axis(counts)
        self.assertTrue(np.allclose(out[0], out[1]))

    def test_negative_counts_raise(self) -> None:
        with self.assertRaises(NasSourceFitError):
            log2_cpm_complete_axis(np.asarray([[-1.0, 1.0]]))


class PipelineTests(unittest.TestCase):
    def _fixture(self):
        rng = np.random.default_rng(29)
        targets = rng.integers(0, 9, size=90).astype(float)
        rep = rng.normal(size=(90, 60))
        rep[:, :5] += targets[:, None] * 0.8
        return targets, rep, np.ones_like(rep), [f"ENSG{i:05d}" for i in range(60)]

    def test_fit_then_apply_is_row_independent(self) -> None:
        targets, rep, elig, ids = self._fixture()
        config = {"representation": "gene_median", "reducer": "none"}
        state = fit_pipeline(
            config=config, representation=rep, eligibility=elig, gene_ids=ids,
            fitting=np.arange(90, dtype=np.int64), targets=targets,
            feature_count=20, hyperparameters={"alpha": 10.0}, seed=1701,
        )
        state.pop("_design_columns")
        whole = apply_pipeline(state=state, representation=rep)
        for index in (0, 17, 89):
            alone = apply_pipeline(state=state, representation=rep[index:index + 1])
            self.assertAlmostEqual(float(alone[0]), float(whole[index]), places=12)

    def test_ridge_always_keeps_coefficients(self) -> None:
        """Ridge cannot produce the all-zero-coefficient collapse seen before."""

        targets, rep, elig, ids = self._fixture()
        for alpha in (1.0, 100.0, 10000.0, 1e8):
            state = fit_pipeline(
                config={"representation": "gene_median", "reducer": "none"},
                representation=rep, eligibility=elig, gene_ids=ids,
                fitting=np.arange(90, dtype=np.int64), targets=targets,
                feature_count=20, hyperparameters={"alpha": alpha}, seed=1701,
            )
            self.assertGreater(int(state["nonzero_coefficients"][0]), 0)

    def test_pca_reducer_shapes(self) -> None:
        targets, rep, elig, ids = self._fixture()
        state = fit_pipeline(
            config={"representation": "gene_median", "reducer": "pca"},
            representation=rep, eligibility=elig, gene_ids=ids,
            fitting=np.arange(90, dtype=np.int64), targets=targets,
            feature_count=20, hyperparameters={"alpha": 10.0, "pca_components": 5},
            seed=1701,
        )
        self.assertEqual(state["pca_components"].shape, (5, 20))
        self.assertEqual(state["coef"].shape, (1, 5))

    def test_unregistered_reducer_raises(self) -> None:
        targets, rep, elig, ids = self._fixture()
        with self.assertRaises(NasSourceFitError):
            fit_pipeline(
                config={"representation": "gene_median", "reducer": "umap"},
                representation=rep, eligibility=elig, gene_ids=ids,
                fitting=np.arange(90, dtype=np.int64), targets=targets,
                feature_count=20, hyperparameters={"alpha": 1.0}, seed=1,
            )

    def test_selection_returns_a_choice_and_grid_evidence(self) -> None:
        targets, rep, elig, ids = self._fixture()
        folds = np.asarray([i % 5 for i in range(90)], dtype=np.int64)
        chosen, rows = select_hyperparameters(
            model_id="gene_median_ridge",
            config={"representation": "gene_median", "reducer": "none"},
            grid=[{"alpha": 10.0}, {"alpha": 1000.0}],
            representation=rep, eligibility=elig, gene_ids=ids, outer_folds=folds,
            targets=targets, feature_count=20, seed=1701,
            null_replicates=200, full_indices=np.arange(90, dtype=np.int64),
        )
        self.assertIsNotNone(chosen)
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(row["eligible"] for row in rows))
        self.assertIn("source_pooled_oof_excess_over_null_mean", rows[0])


class FeatureSelectionTests(unittest.TestCase):
    def test_ineligible_genes_are_dropped(self) -> None:
        elig = np.zeros((20, 5)); elig[:, [1, 3]] = 1.0
        rep = np.tile(np.arange(20)[:, None], (1, 5)) * np.asarray([9.0, 1.0, 8.0, 2.0, 7.0])
        got = select_feature_indices(
            eligibility=elig, representation=rep,
            gene_ids=[f"G{i}" for i in range(5)], feature_count=2,
        )
        self.assertEqual(sorted(got.tolist()), [1, 3])


if __name__ == "__main__":
    unittest.main()

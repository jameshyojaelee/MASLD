from __future__ import annotations

import unittest

import numpy as np

from scripts.fit_gse267145_fibrosis_transfer_source_models import (
    SourceFibrosisFitError,
    apply_pipeline,
    assert_row_independent,
    average_precision_null,
    build_labels,
    fit_pipeline,
    gene_median_center,
    log2_cpm_complete_axis,
    out_of_fold_scores,
    per_array_rank,
    select_feature_indices,
    select_hyperparameters,
)


MAPPING = {
    "source_field": "fibrosis",
    "source_scale": [0, 1, 2, 3],
    "mild_f0_f1": [0, 1],
    "advanced_f3_f4": [3],
    "excluded": [2],
}


def endpoint_rows(values):
    return [{"participant_id": f"P{index}", "fibrosis": str(value)} for index, value in enumerate(values)]


class RepresentationTests(unittest.TestCase):
    def test_rank_is_scale_free_within_an_array(self) -> None:
        values = np.asarray([[1.0, 5.0, 2.0, 9.0], [10.0, 50.0, 20.0, 90.0]])
        ranks = per_array_rank(values)
        # A per-array monotone rescale must not move any rank.
        self.assertTrue(np.allclose(ranks[0], ranks[1]))
        self.assertTrue(np.all(ranks > 0.0) and np.all(ranks < 1.0))

    def test_rank_handles_ties_by_average(self) -> None:
        values = np.asarray([[3.0, 3.0, 1.0, 4.0]])
        ranks = per_array_rank(values)
        self.assertAlmostEqual(float(ranks[0, 0]), float(ranks[0, 1]))
        self.assertLess(float(ranks[0, 2]), float(ranks[0, 0]))

    def test_gene_median_removes_only_the_array_offset(self) -> None:
        values = np.asarray([[1.0, 2.0, 3.0, 4.0]])
        shifted = values + 7.5
        self.assertTrue(
            np.allclose(gene_median_center(values), gene_median_center(shifted))
        )

    def test_representations_are_row_independent(self) -> None:
        rng = np.random.default_rng(11)
        values = rng.normal(size=(6, 25))
        for name in ("per_array_rank", "gene_median"):
            assert_row_independent(name, values)

    def test_row_independence_check_rejects_a_pooled_transform(self) -> None:
        import scripts.fit_gse267145_fibrosis_transfer_source_models as module

        rng = np.random.default_rng(3)
        values = rng.normal(size=(6, 25))
        original = module.REPRESENTATION_FUNCTIONS["gene_median"]
        try:
            module.REPRESENTATION_FUNCTIONS["gene_median"] = (
                lambda matrix: np.asarray(matrix) - np.median(np.asarray(matrix))
            )
            with self.assertRaises(SourceFibrosisFitError):
                assert_row_independent("gene_median", values)
        finally:
            module.REPRESENTATION_FUNCTIONS["gene_median"] = original

    def test_log2_cpm_uses_the_complete_axis_total(self) -> None:
        counts = np.asarray([[1.0, 1.0, 2.0], [10.0, 10.0, 20.0]])
        transformed = log2_cpm_complete_axis(counts)
        # Proportional libraries must land on the same CPM.
        self.assertTrue(np.allclose(transformed[0], transformed[1]))

    def test_negative_source_values_raise(self) -> None:
        with self.assertRaises(SourceFibrosisFitError):
            log2_cpm_complete_axis(np.asarray([[-1.0, 2.0]]))


class EndpointMappingTests(unittest.TestCase):
    def test_f2_is_excluded_from_both_arms(self) -> None:
        labels, keep, audit = build_labels(
            endpoint_rows=endpoint_rows([0, 1, 2, 3, 3, 0]), mapping=MAPPING
        )
        self.assertEqual(list(labels), [0, 0, 1, 1, 0])
        self.assertEqual(list(keep), [0, 1, 3, 4, 5])
        self.assertEqual(audit["excluded_participants"], 1)
        self.assertEqual(
            audit["class_counts"], {"mild_f0_f1": 3, "advanced_f3_f4": 2}
        )
        self.assertIs(audit["source_scale_has_stage_4"], False)

    def test_prevalence_comes_from_the_fitting_partition(self) -> None:
        _, _, audit = build_labels(
            endpoint_rows=endpoint_rows([0, 0, 0, 2, 3]), mapping=MAPPING
        )
        # F2 is dropped before the denominator, so 1 of 4, not 1 of 5.
        self.assertAlmostEqual(audit["training_prevalence_advanced_f3_f4"], 0.25)

    def test_off_scale_value_raises(self) -> None:
        with self.assertRaises(SourceFibrosisFitError):
            build_labels(endpoint_rows=endpoint_rows([0, 4]), mapping=MAPPING)

    def test_overlapping_mapping_raises(self) -> None:
        mapping = dict(MAPPING, advanced_f3_f4=[1, 3])
        with self.assertRaises(SourceFibrosisFitError):
            build_labels(endpoint_rows=endpoint_rows([0, 1, 3]), mapping=mapping)

    def test_single_class_mapping_raises(self) -> None:
        with self.assertRaises(SourceFibrosisFitError):
            build_labels(endpoint_rows=endpoint_rows([0, 0, 1]), mapping=MAPPING)


class SelectionAndPipelineTests(unittest.TestCase):
    def test_selection_ranks_by_variance_and_breaks_ties_by_gene_id(self) -> None:
        eligibility = np.ones((10, 4))
        representation = np.zeros((10, 4))
        representation[:, 0] = np.arange(10)
        representation[:, 1] = np.arange(10) * 0.5
        representation[:, 2] = np.arange(10) * 0.5
        selected = select_feature_indices(
            eligibility_matrix=eligibility,
            representation=representation,
            gene_ids=["ENSG4", "ENSG1", "ENSG2", "ENSG3"],
            feature_count=2,
        )
        # Highest variance first, then the alphabetically smaller of the tie.
        self.assertEqual(sorted(selected.tolist()), [0, 1])

    def test_ineligible_genes_are_dropped_before_ranking(self) -> None:
        eligibility = np.zeros((10, 4))
        eligibility[:, [1, 2]] = 1.0
        representation = np.tile(np.arange(10)[:, None], (1, 4)) * np.asarray([5.0, 1.0, 2.0, 4.0])
        selected = select_feature_indices(
            eligibility_matrix=eligibility,
            representation=representation,
            gene_ids=["ENSG1", "ENSG2", "ENSG3", "ENSG4"],
            feature_count=2,
        )
        self.assertEqual(sorted(selected.tolist()), [1, 2])

    def test_too_few_eligible_genes_raises(self) -> None:
        with self.assertRaises(SourceFibrosisFitError):
            select_feature_indices(
                eligibility_matrix=np.zeros((10, 4)),
                representation=np.ones((10, 4)),
                gene_ids=["A", "B", "C", "D"],
                feature_count=2,
            )

    def test_fit_then_apply_reproduces_the_fitted_probabilities(self) -> None:
        rng = np.random.default_rng(7)
        labels = np.asarray([0] * 20 + [1] * 6, dtype=np.int64)
        representation = rng.normal(size=(26, 30))
        representation[labels == 1] += 1.5
        eligibility = np.ones_like(representation)
        gene_ids = [f"ENSG{index:05d}" for index in range(30)]
        config = {"representation": "gene_median", "reducer": "pca", "classifier": "elastic_net"}
        state = fit_pipeline(
            model_id="gene_median_pca_elastic_net",
            config=config,
            representation=representation,
            eligibility=eligibility,
            gene_ids=gene_ids,
            fitting=np.arange(26, dtype=np.int64),
            labels=labels,
            feature_count=10,
            pca_components=3,
            hyperparameters={"c": 1.0, "l1_ratio": 0.5},
            seed=1701,
        )
        self.assertTrue(state.pop("_converged"))
        state.pop("_design_columns")
        probabilities = apply_pipeline(
            config=config, state=state, representation=representation
        )
        self.assertEqual(probabilities.shape, (26,))
        self.assertTrue(np.all(probabilities >= 0.0) and np.all(probabilities <= 1.0))
        # Applying to a subset must equal the matching slice of the full apply.
        subset = apply_pipeline(
            config=config, state=state, representation=representation[[3, 1, 25]]
        )
        self.assertTrue(np.allclose(subset, probabilities[[3, 1, 25]]))

    def test_unregistered_reducer_raises(self) -> None:
        with self.assertRaises(SourceFibrosisFitError):
            fit_pipeline(
                model_id="x",
                config={"representation": "gene_median", "reducer": "umap", "classifier": "elastic_net"},
                representation=np.ones((10, 5)),
                eligibility=np.ones((10, 5)),
                gene_ids=[f"G{i}" for i in range(5)],
                fitting=np.arange(10, dtype=np.int64),
                labels=np.asarray([0] * 5 + [1] * 5, dtype=np.int64),
                feature_count=3,
                pca_components=2,
                hyperparameters={"c": 1.0, "l1_ratio": 0.5},
                seed=1,
            )


if __name__ == "__main__":
    unittest.main()


class NullReferencedSelectionTests(unittest.TestCase):
    """The v1 selection let an all-zero L1 fit win; these pin the v2 guards."""

    def _fixture(self):
        rng = np.random.default_rng(19)
        # 90 participants, 4 positives: the real source composition.
        labels = np.asarray([0] * 86 + [1] * 4, dtype=np.int64)
        representation = rng.normal(size=(90, 40))
        representation[labels == 1] += 2.0
        eligibility = np.ones_like(representation)
        gene_ids = [f"ENSG{index:05d}" for index in range(40)]
        outer_folds = np.asarray([index % 5 for index in range(90)], dtype=np.int64)
        # Make sure every fold's training partition keeps both classes.
        outer_folds[labels == 1] = np.asarray([0, 1, 2, 3])
        return labels, representation, eligibility, gene_ids, outer_folds

    def test_null_mean_for_a_continuous_scorer_exceeds_prevalence(self) -> None:
        rng = np.random.default_rng(23)
        labels = np.asarray([0] * 86 + [1] * 4, dtype=np.int64)
        scores = rng.uniform(size=90)
        null = average_precision_null(labels, scores, replicates=2000, seed=5)
        self.assertGreater(null["null_mean"], float(np.mean(labels)))

    def test_all_zero_coefficient_configuration_is_ineligible(self) -> None:
        labels, representation, eligibility, gene_ids, outer_folds = self._fixture()
        config = {
            "representation": "gene_median",
            "reducer": "none",
            "classifier": "elastic_net",
        }
        # C = 1e-6 with pure L1 shrinks every coefficient to exactly zero.
        chosen, rows = select_hyperparameters(
            model_id="gene_median_elastic_net",
            config=config,
            grid=[{"c": 1e-6, "l1_ratio": 1.0}, {"c": 1.0, "l1_ratio": 0.5}],
            representation=representation,
            eligibility=eligibility,
            gene_ids=gene_ids,
            outer_folds=outer_folds,
            labels=labels,
            feature_count=20,
            pca_components=3,
            seed=1701,
            null_replicates=200,
            full_indices=np.arange(90, dtype=np.int64),
        )
        collapsed = next(row for row in rows if row["c"] == 1e-6)
        working = next(row for row in rows if row["c"] == 1.0)
        self.assertEqual(collapsed["full_fit_nonzero_coefficients"], 0)
        self.assertIs(collapsed["eligible"], False)
        self.assertIs(working["eligible"], True)
        # The collapsed fit must not be chosen even though c ascends first.
        self.assertEqual(chosen["c"], 1.0)

    def test_every_configuration_degenerate_returns_no_choice(self) -> None:
        """An unfittable model is a recorded result, not an exception.

        Raising here would push the caller toward widening the grid until
        something survives.  Returning None lets the fit record the model as
        unfittable with its evidence and carry on.
        """

        labels, representation, eligibility, gene_ids, outer_folds = self._fixture()
        chosen, rows = select_hyperparameters(
            model_id="gene_median_elastic_net",
            config={
                "representation": "gene_median",
                "reducer": "none",
                "classifier": "elastic_net",
            },
            grid=[{"c": 1e-8, "l1_ratio": 1.0}],
            representation=representation,
            eligibility=eligibility,
            gene_ids=gene_ids,
            outer_folds=outer_folds,
            labels=labels,
            feature_count=20,
            pca_components=3,
            seed=1701,
            null_replicates=100,
            full_indices=np.arange(90, dtype=np.int64),
        )
        self.assertIsNone(chosen)
        self.assertEqual(rows[0]["eligible_configurations_for_this_model"], 0)
        self.assertIs(rows[0]["eligible"], False)

    def test_out_of_fold_reports_fold_level_coefficient_use(self) -> None:
        labels, representation, eligibility, gene_ids, outer_folds = self._fixture()
        config = {
            "representation": "gene_median",
            "reducer": "none",
            "classifier": "elastic_net",
        }
        _, converged, uses_data = out_of_fold_scores(
            model_id="gene_median_elastic_net",
            config=config,
            representation=representation,
            eligibility=eligibility,
            gene_ids=gene_ids,
            outer_folds=outer_folds,
            labels=labels,
            feature_count=20,
            pca_components=3,
            hyperparameters={"c": 1e-6, "l1_ratio": 1.0},
            seed=1701,
        )
        self.assertTrue(converged)
        self.assertIs(uses_data, False)


class SeparabilityTests(unittest.TestCase):
    def test_separable_classes_are_detected(self) -> None:
        from scripts.fit_gse267145_fibrosis_transfer_source_models import (
            classes_are_linearly_separable,
        )

        design = np.asarray([[0.0], [0.1], [5.0], [5.1]])
        labels = np.asarray([0, 0, 1, 1], dtype=np.int64)
        self.assertIs(classes_are_linearly_separable(design, labels), True)

    def test_interleaved_classes_are_not_separable(self) -> None:
        from scripts.fit_gse267145_fibrosis_transfer_source_models import (
            classes_are_linearly_separable,
        )

        design = np.asarray([[0.0], [1.0], [2.0], [3.0]])
        labels = np.asarray([0, 1, 0, 1], dtype=np.int64)
        self.assertIs(classes_are_linearly_separable(design, labels), False)

    def test_four_positives_in_five_dimensions_are_not_automatically_separable(
        self,
    ) -> None:
        """Separability at 4-of-90 in five dimensions is a fact about the data.

        Four positives sitting inside the negatives' convex hull are not
        separable, so this has to be measured on the real source rather than
        assumed from the class counts.
        """

        from scripts.fit_gse267145_fibrosis_transfer_source_models import (
            classes_are_linearly_separable,
        )

        rng = np.random.default_rng(41)
        design = rng.normal(size=(90, 5))
        labels = np.asarray([0] * 86 + [1] * 4, dtype=np.int64)
        self.assertIs(classes_are_linearly_separable(design, labels), False)

    def test_ninety_points_in_a_thousand_dimensions_are_separable(self) -> None:
        """Where the fitted design has more columns than participants, though,
        a separating hyperplane always exists and the likelihood is unbounded."""

        from scripts.fit_gse267145_fibrosis_transfer_source_models import (
            classes_are_linearly_separable,
        )

        rng = np.random.default_rng(43)
        design = rng.normal(size=(90, 1000))
        labels = np.asarray([0] * 86 + [1] * 4, dtype=np.int64)
        self.assertIs(classes_are_linearly_separable(design, labels), True)

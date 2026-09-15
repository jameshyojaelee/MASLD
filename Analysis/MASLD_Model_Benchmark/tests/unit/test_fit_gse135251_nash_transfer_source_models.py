from __future__ import annotations

import unittest

import numpy as np

from scripts.fit_gse135251_nash_transfer_source_models import (
    CLASSES,
    LEARNED_MODEL_IDS,
    POSITIVE_CLASS,
    SourceNashFitError,
    apply_pipeline,
    assert_row_independent,
    average_precision_null,
    build_labels,
    build_representation,
    class_weights,
    classes_are_linearly_separable,
    fit_pipeline,
    gene_median_center,
    log2_cpm_complete_axis,
    out_of_fold_scores,
    per_array_rank,
    select_feature_indices,
    select_hyperparameters,
)

MAPPING = {
    "source_field": "nash_status__nash_vs_nafl",
    "positive_value": "nash",
    "negative_value": "no_nash",
    "excluded_value": "excluded",
    "positive_groups": ["NASH_F0-F1", "NASH_F2", "NASH_F3", "NASH_F4"],
    "negative_groups": ["NAFL"],
    "excluded_groups": ["control"],
    "mapping_rule": "controls excluded",
}


def toy_source(*, participants: int = 60, genes: int = 40, seed: int = 3):
    rng = np.random.default_rng(seed)
    counts = rng.integers(1, 500, size=(participants, genes)).astype(np.float64)
    labels = np.asarray([1] * (participants // 2) + [0] * (participants // 2))
    counts[labels == 1, :5] += 4000.0
    return counts, labels


class RepresentationTests(unittest.TestCase):
    def test_rank_and_gene_median_are_per_array_functions(self) -> None:
        values = np.asarray(
            [[1.0, 2.0, 3.0, 4.0], [8.0, 6.0, 7.0, 5.0], [0.5, 0.25, 2.0, 1.0]]
        )
        for name in ("per_array_rank", "gene_median"):
            assert_row_independent(name, values)
            built = build_representation(name=name, log2_shared=values)
            self.assertEqual(built.shape, values.shape)

    def test_one_row_alone_reproduces_its_row_in_the_full_matrix(self) -> None:
        values = np.asarray([[1.0, 5.0, 3.0], [9.0, 2.0, 7.0], [4.0, 4.0, 8.0]])
        for function in (per_array_rank, gene_median_center):
            full = function(values)
            for index in range(values.shape[0]):
                alone = function(values[index : index + 1])
                self.assertTrue(np.array_equal(alone[0], full[index]))

    def test_a_pooled_transform_is_caught(self) -> None:
        """A column-centred transform is not per-array and must be rejected."""

        import scripts.fit_gse135251_nash_transfer_source_models as module

        values = np.asarray([[1.0, 5.0, 3.0], [9.0, 2.0, 7.0], [4.0, 4.0, 8.0]])
        original = module.REPRESENTATION_FUNCTIONS["gene_median"]
        module.REPRESENTATION_FUNCTIONS["gene_median"] = lambda x: np.asarray(
            x, dtype=np.float64
        ) - np.mean(np.asarray(x, dtype=np.float64), axis=0, keepdims=True)
        try:
            with self.assertRaises(SourceNashFitError):
                assert_row_independent("gene_median", values)
        finally:
            module.REPRESENTATION_FUNCTIONS["gene_median"] = original

    def test_log2_cpm_uses_a_base_two_logarithm(self) -> None:
        counts = np.asarray([[1_000_000.0, 0.0]])
        transformed = log2_cpm_complete_axis(counts)
        self.assertAlmostEqual(float(transformed[0, 0]), np.log2(1_000_001.0))
        self.assertAlmostEqual(float(transformed[0, 1]), 0.0)

    def test_a_negative_or_empty_library_is_rejected(self) -> None:
        with self.assertRaises(SourceNashFitError):
            log2_cpm_complete_axis(np.asarray([[-1.0, 2.0]]))
        with self.assertRaises(SourceNashFitError):
            log2_cpm_complete_axis(np.asarray([[0.0, 0.0]]))


class LabelTests(unittest.TestCase):
    def _rows(self):
        rows = []
        for _ in range(131):
            rows.append({"nash_status__nash_vs_nafl": "nash"})
        for _ in range(41):
            rows.append({"nash_status__nash_vs_nafl": "no_nash"})
        for _ in range(8):
            rows.append({"nash_status__nash_vs_nafl": "excluded"})
        return rows

    def test_the_excluded_group_is_dropped_not_folded_into_the_negative(self) -> None:
        labels, keep, audit = build_labels(
            endpoint_rows=self._rows(), mapping=MAPPING
        )
        self.assertEqual(len(labels), 172)
        self.assertEqual(len(keep), 172)
        self.assertEqual(audit["class_counts"], {"no_nash": 41, "nash": 131})
        self.assertEqual(audit["excluded_participants"], 8)
        self.assertAlmostEqual(audit["training_prevalence_nash"], 131 / 172)

    def test_an_unregistered_source_value_is_rejected(self) -> None:
        rows = self._rows()
        rows[0]["nash_status__nash_vs_nafl"] = "borderline"
        with self.assertRaises(SourceNashFitError):
            build_labels(endpoint_rows=rows, mapping=MAPPING)

    def test_a_single_class_source_is_rejected(self) -> None:
        rows = [{"nash_status__nash_vs_nafl": "nash"} for _ in range(10)]
        with self.assertRaises(SourceNashFitError):
            build_labels(endpoint_rows=rows, mapping=MAPPING)

    def test_the_positive_class_is_nash(self) -> None:
        self.assertEqual(CLASSES, ("no_nash", "nash"))
        self.assertEqual(POSITIVE_CLASS, "nash")


class NullTests(unittest.TestCase):
    def test_a_constant_scorer_gets_a_narrower_null_than_a_continuous_one(self) -> None:
        """This is the whole reason selection uses excess over each config's own null."""

        rng = np.random.default_rng(11)
        labels = np.asarray([1] * 40 + [0] * 60)
        constant = np.zeros(100)
        continuous = rng.normal(size=100)
        flat = average_precision_null(labels, constant, replicates=400, seed=5)
        wide = average_precision_null(labels, continuous, replicates=400, seed=5)
        self.assertLess(flat["null_p95"], wide["null_p95"])

    def test_the_null_mean_is_not_prevalence(self) -> None:
        rng = np.random.default_rng(7)
        labels = np.asarray([1] * 15 + [0] * 22)
        scores = rng.normal(size=37)
        null = average_precision_null(labels, scores, replicates=2000, seed=13)
        self.assertGreater(null["null_mean"], float(np.mean(labels)) + 0.02)


class PipelineTests(unittest.TestCase):
    CONFIG = {
        "classifier": "elastic_net",
        "reducer": "none",
        "representation": "gene_median",
    }

    def _fixture(self):
        counts, labels = toy_source()
        log2 = log2_cpm_complete_axis(counts)
        representation = build_representation(name="gene_median", log2_shared=log2)
        return counts, labels, representation

    def test_feature_selection_uses_only_the_fitting_partition(self) -> None:
        counts, _, representation = self._fixture()
        gene_ids = [f"ENSG{index:011d}" for index in range(counts.shape[1])]
        first = select_feature_indices(
            eligibility_matrix=counts[:30],
            representation=representation[:30],
            gene_ids=gene_ids,
            feature_count=10,
        )
        second = select_feature_indices(
            eligibility_matrix=counts[:30],
            representation=representation[:30],
            gene_ids=gene_ids,
            feature_count=10,
        )
        different = select_feature_indices(
            eligibility_matrix=counts[30:],
            representation=representation[30:],
            gene_ids=gene_ids,
            feature_count=10,
        )
        self.assertTrue(np.array_equal(first, second))
        self.assertFalse(np.array_equal(first, different))

    def test_apply_pipeline_returns_probabilities(self) -> None:
        counts, labels, representation = self._fixture()
        gene_ids = [f"ENSG{index:011d}" for index in range(counts.shape[1])]
        state = fit_pipeline(
            config=self.CONFIG,
            representation=representation,
            eligibility=counts,
            gene_ids=gene_ids,
            fitting=np.arange(len(labels)),
            labels=labels,
            feature_count=10,
            pca_components=3,
            hyperparameters={"c": 1.0, "l1_ratio": 0.0},
            seed=1701,
        )
        state.pop("_converged")
        state.pop("_design_columns")
        state.pop("_design")
        probabilities = apply_pipeline(
            config=self.CONFIG, state=state, representation=representation
        )
        self.assertEqual(probabilities.shape, (len(labels),))
        self.assertTrue(np.all((probabilities >= 0.0) & (probabilities <= 1.0)))

    def test_an_all_zero_coefficient_configuration_is_ineligible(self) -> None:
        """L1 at a tiny C erases the model; it must not be selectable."""

        counts, labels, representation = self._fixture()
        gene_ids = [f"ENSG{index:011d}" for index in range(counts.shape[1])]
        folds = np.arange(len(labels)) % 5
        chosen, rows = select_hyperparameters(
            model_id="gene_median_elastic_net",
            config=self.CONFIG,
            grid=[{"c": 1e-8, "l1_ratio": 1.0}],
            representation=representation,
            eligibility=counts,
            gene_ids=gene_ids,
            outer_folds=folds,
            labels=labels,
            feature_count=10,
            pca_components=3,
            seed=1701,
            null_replicates=50,
            full_indices=np.arange(len(labels)),
        )
        self.assertIsNone(chosen)
        self.assertEqual(rows[0]["full_fit_nonzero_coefficients"], 0)
        self.assertFalse(rows[0]["eligible"])

    def test_out_of_fold_scores_cover_every_participant(self) -> None:
        counts, labels, representation = self._fixture()
        gene_ids = [f"ENSG{index:011d}" for index in range(counts.shape[1])]
        folds = np.arange(len(labels)) % 5
        scores, converged, uses_data = out_of_fold_scores(
            config=self.CONFIG,
            representation=representation,
            eligibility=counts,
            gene_ids=gene_ids,
            outer_folds=folds,
            labels=labels,
            feature_count=10,
            pca_components=3,
            hyperparameters={"c": 1.0, "l1_ratio": 0.0},
            seed=1701,
        )
        self.assertEqual(len(scores), len(labels))
        self.assertTrue(np.all(np.isfinite(scores)))
        self.assertTrue(uses_data)
        self.assertIsInstance(converged, bool)


class SeparabilityTests(unittest.TestCase):
    def test_separable_and_inseparable_designs_are_distinguished(self) -> None:
        labels = np.asarray([1, 1, 0, 0])
        separable = np.asarray([[3.0], [4.0], [-3.0], [-4.0]])
        inseparable = np.asarray([[1.0], [-1.0], [1.0], [-1.0]])
        self.assertTrue(classes_are_linearly_separable(separable, labels))
        self.assertFalse(classes_are_linearly_separable(inseparable, labels))


class WeightTests(unittest.TestCase):
    def test_class_weights_are_inverse_frequency(self) -> None:
        labels = np.asarray([1] * 131 + [0] * 41)
        weights = class_weights(labels)
        self.assertAlmostEqual(weights[1], 172 / (2 * 131))
        self.assertAlmostEqual(weights[0], 172 / (2 * 41))

    def test_a_single_class_partition_is_rejected(self) -> None:
        with self.assertRaises(SourceNashFitError):
            class_weights(np.asarray([1, 1, 1]))


class RosterTests(unittest.TestCase):
    def test_the_learned_roster_matches_the_taskspec_molecular_baselines(self) -> None:
        self.assertEqual(
            sorted(LEARNED_MODEL_IDS),
            [
                "gene_median_elastic_net",
                "gene_median_linear_svm",
                "gene_median_pca_elastic_net",
                "per_array_rank_elastic_net",
            ],
        )


if __name__ == "__main__":
    unittest.main()

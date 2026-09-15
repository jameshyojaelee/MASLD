"""The separability analysis must agree with the package it claims to reuse.

The analysis vectorizes Spearman over 42,163 genes because a pure-Python loop
over that universe is not practical. Vectorizing is where a rank correlation
quietly stops handling ties, so these tests pin the vectorized path against
``masld_bench.evaluators.metrics`` rather than against a fresh reimplementation
of the same idea.

They also pin the arithmetic of the criteria themselves: the participation
ratio's closed form and, with it, the fact that C2 binds harder than C1 for
homogeneous correlations. That is a property of the frozen thresholds worth
knowing before the answer arrives, not after.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np

from masld_bench.evaluators.metrics import _average_ranks, spearman_correlation


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "scripts" / "evaluate_gse267145_aspect_separability.py"
FREEZER = ROOT / "scripts" / "freeze_aspect_separability_prespecification.py"
ENDPOINTS = (
    ROOT
    / "executions"
    / "model-data-061-21079623"
    / "activation"
    / "participant_endpoints.tsv"
)


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


analysis = _load(MODULE, "aspect_separability_analysis_tested")
freezer = _load(FREEZER, "aspect_separability_freezer_tested")


class ThresholdsMatchTheFreezerTests(unittest.TestCase):
    """The analysis must expect exactly the literals the freezer writes."""

    def test_expected_threshold_prose_is_the_frozen_prose(self) -> None:
        criteria = freezer.build()["criteria"]
        self.assertEqual(
            {name: criterion["threshold"] for name, criterion in criteria.items()},
            analysis.FROZEN_THRESHOLDS,
        )

    def test_numeric_thresholds_agree_with_the_prose(self) -> None:
        self.assertIn("0.80", analysis.FROZEN_THRESHOLDS["c1_labels_are_not_redundant"])
        self.assertEqual(analysis.C1_THRESHOLD, 0.80)
        self.assertIn(
            "2.0", analysis.FROZEN_THRESHOLDS["c2_effective_dimensionality"]
        )
        self.assertEqual(analysis.C2_THRESHOLD, 2.0)
        self.assertIn(
            "0.80", analysis.FROZEN_THRESHOLDS["c3_aspects_rank_genes_differently"]
        )
        self.assertEqual(analysis.C3_THRESHOLD, 0.80)

    def test_necrosis_is_not_one_of_the_three_aspects(self) -> None:
        self.assertEqual(
            analysis.ASPECTS,
            ("steatosis", "ballooning", "lobular_inflammation"),
        )
        self.assertNotIn("lobular_necrosis", analysis.ASPECTS)
        self.assertNotIn("fibrosis", analysis.ASPECTS)
        self.assertIn("lobular_necrosis", analysis.CONTEXT_AXES)
        self.assertIn("fibrosis", analysis.CONTEXT_AXES)

    def test_the_three_pairs_are_the_three_unordered_pairs(self) -> None:
        self.assertEqual(len(analysis.PAIRS), 3)
        self.assertEqual(
            {frozenset(pair) for pair in analysis.PAIRS},
            {
                frozenset(("steatosis", "ballooning")),
                frozenset(("steatosis", "lobular_inflammation")),
                frozenset(("ballooning", "lobular_inflammation")),
            },
        )


class ParticipationRatioTests(unittest.TestCase):
    def test_uncorrelated_axes_span_three_dimensions(self) -> None:
        ratio, _ = analysis.participation_ratio(np.eye(3))
        self.assertAlmostEqual(ratio, 3.0, places=12)

    def test_perfectly_correlated_axes_span_one(self) -> None:
        ratio, _ = analysis.participation_ratio(np.ones((3, 3)))
        self.assertAlmostEqual(ratio, 1.0, places=12)

    def test_the_closed_form_holds(self) -> None:
        for correlations in ((0.1, 0.2, 0.3), (0.5, 0.5, 0.5), (0.9, 0.8, 0.7)):
            r12, r13, r23 = correlations
            matrix = np.array(
                [[1.0, r12, r13], [r12, 1.0, r23], [r13, r23, 1.0]]
            )
            ratio, _ = analysis.participation_ratio(matrix)
            closed = 9.0 / (3.0 + 2.0 * (r12**2 + r13**2 + r23**2))
            with self.subTest(correlations=correlations):
                self.assertAlmostEqual(ratio, closed, places=12)

    def test_c2_binds_at_half_for_homogeneous_correlations(self) -> None:
        """C2 at 2.0 is stricter than C1 at 0.80, and the boundary is 0.5."""

        matrix = np.full((3, 3), 0.5)
        np.fill_diagonal(matrix, 1.0)
        ratio, _ = analysis.participation_ratio(matrix)
        self.assertAlmostEqual(ratio, analysis.C2_THRESHOLD, places=12)
        stricter = np.full((3, 3), 0.55)
        np.fill_diagonal(stricter, 1.0)
        tighter, _ = analysis.participation_ratio(stricter)
        self.assertLess(tighter, analysis.C2_THRESHOLD)
        self.assertLess(0.55, analysis.C1_THRESHOLD)


class VectorizedRankAgreementTests(unittest.TestCase):
    """Ties are where a vectorized rank correlation goes wrong silently."""

    def setUp(self) -> None:
        generator = np.random.default_rng(11)
        # Deliberately coarse, so most values are tied, like these endpoints.
        self.matrix = generator.integers(0, 4, size=(99, 40)).astype(float)
        self.matrix[:, 0] = generator.random(99)

    def test_average_ranks_match_the_package(self) -> None:
        ranked = analysis.average_ranks(self.matrix)
        for column in range(self.matrix.shape[1]):
            reference = np.asarray(
                _average_ranks(list(self.matrix[:, column])), dtype=float
            )
            with self.subTest(column=column):
                np.testing.assert_allclose(ranked[:, column], reference, atol=0.0)

    def test_validate_rank_agreement_reports_exact_agreement(self) -> None:
        report = analysis.validate_rank_agreement(self.matrix, range(10))
        self.assertTrue(report["agrees_with_masld_bench_average_ranks"])
        self.assertEqual(report["max_absolute_rank_difference"], 0.0)

    def test_association_vector_matches_the_package_spearman(self) -> None:
        outcome = np.asarray([0, 1, 2] * 33, dtype=float)
        standardized = analysis.standardize(analysis.average_ranks(self.matrix))
        vector = analysis.association_vector(standardized, outcome)
        for column in range(self.matrix.shape[1]):
            expected = spearman_correlation(
                self.matrix[:, column].tolist(), outcome.tolist()
            )
            with self.subTest(column=column):
                self.assertAlmostEqual(float(vector[column]), expected, places=12)

    def test_spearman_of_vectors_matches_the_package(self) -> None:
        generator = np.random.default_rng(5)
        left = generator.random(500)
        right = 0.6 * left + 0.4 * generator.random(500)
        self.assertAlmostEqual(
            analysis.spearman_of_vectors(left, right),
            spearman_correlation(left.tolist(), right.tolist()),
            places=12,
        )

    def test_a_constant_column_is_refused_rather_than_silently_zeroed(self) -> None:
        constant = np.ones((99, 2))
        with self.assertRaises(analysis.SeparabilityError):
            analysis.standardize(analysis.average_ranks(constant))


class TieStructureTests(unittest.TestCase):
    def test_a_known_vector(self) -> None:
        report = analysis.tie_structure([0] * 4 + [1] * 2, "toy")
        self.assertEqual(report["distinct_values"], 2)
        self.assertEqual(report["largest_tied_block"], 4)
        self.assertEqual(report["total_pairs"], 15)
        self.assertEqual(report["tied_pairs"], 6 + 1)

    def test_the_deposited_profile_ties(self) -> None:
        rows = analysis.read_table(ENDPOINTS)
        report = analysis.joint_tie_structure(rows)
        self.assertEqual(report["n_participants"], 99)
        self.assertEqual(report["all_zero_profile_count"], 24)
        self.assertGreaterEqual(report["largest_tied_block"], 24)


class NASIdentityTests(unittest.TestCase):
    def test_the_component_sum_excludes_necrosis(self) -> None:
        rows = analysis.read_table(ENDPOINTS)
        report = analysis.verify_nas_identity(rows)
        self.assertEqual(report["component_sum_equals_three_aspects"], 99)
        self.assertEqual(report["adding_necrosis_breaks_the_identity_for"], 19)
        self.assertFalse(report["lobular_necrosis_is_a_nas_component"])

    def test_a_row_count_is_not_a_unit_count_here_but_is_checked(self) -> None:
        rows = analysis.read_table(ENDPOINTS)
        self.assertEqual(len(rows), 99)
        self.assertEqual(len({row["participant_id"] for row in rows}), 99)


class RecordedMarginalGuardTests(unittest.TestCase):
    def test_a_perturbed_marginal_aborts_the_run(self) -> None:
        rows = analysis.read_table(ENDPOINTS)
        prespec = freezer.build()
        self.assertIsNotNone(analysis.verify_recorded_marginals(rows, prespec))
        damaged = [dict(row) for row in rows]
        damaged[0]["steatosis"] = "3" if damaged[0]["steatosis"] != "3" else "0"
        with self.assertRaises(analysis.SeparabilityError):
            analysis.verify_recorded_marginals(damaged, prespec)


if __name__ == "__main__":
    unittest.main()

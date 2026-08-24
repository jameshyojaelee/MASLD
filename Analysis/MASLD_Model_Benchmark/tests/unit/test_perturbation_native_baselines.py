from __future__ import annotations

import unittest

from masld_bench.perturbation_baselines import (
    CellCountRow,
    PerturbationBaselineError,
    additive_features,
    bilinear_features,
    control_mean,
    fit_ridge,
    perturbed_mean,
    predict_linear,
    pseudobulk_counts,
    strict_training_reference,
    systema_matching_mean,
)


class PerturbationNativeBaselineTest(unittest.TestCase):
    def setUp(self) -> None:
        self.rows = [
            CellCountRow("r1", "control", "nt1", "heparg", True, (1, 2)),
            CellCountRow("r1", "control", "nt2", "heparg", True, (2, 1)),
            CellCountRow("r1", "A", "a1", "heparg", False, (4, 2)),
            CellCountRow("r1", "A", "a2", "heparg", False, (1, 3)),
            CellCountRow("r1", "B", "b1", "heparg", False, (3, 6)),
            CellCountRow("r2", "control", "nt1", "heparg", True, (3, 3)),
            CellCountRow("r2", "A", "a1", "heparg", False, (5, 5)),
            CellCountRow("r2", "B", "b1", "heparg", False, (7, 1)),
        ]
        self.bulk = pseudobulk_counts(self.rows)

    def test_pseudobulk_precedes_means_and_cells_do_not_add_replication(self) -> None:
        self.assertEqual(len(self.bulk), 6)
        r1_control = next(row for row in self.bulk if row.biological_unit == "r1" and row.is_control)
        self.assertEqual(r1_control.counts, (3.0, 3.0))
        self.assertEqual(r1_control.cells, 2)
        self.assertEqual(r1_control.guide_families, ("nt1", "nt2"))
        self.assertEqual(control_mean(self.bulk, "heparg"), (3.0, 3.0))
        self.assertEqual(perturbed_mean(self.bulk, "heparg"), (5.0, 4.25))

    def test_missing_matched_control_fails_closed(self) -> None:
        with self.assertRaisesRegex(PerturbationBaselineError, "structurally missing"):
            control_mean(self.bulk, "absent_context")

    def test_ridge_additive_fit_and_prediction_are_deterministic(self) -> None:
        design = [
            additive_features((0,), (0,)),
            additive_features((0,), (1,)),
            additive_features((1,), (0,)),
            additive_features((1,), (1,)),
        ]
        response = [(1, 2), (4, 1), (3, 6), (6, 5)]
        coefficients = fit_ridge(design, response, penalty=1e-9)
        predicted = predict_linear(design, coefficients)
        for observed, expected in zip(predicted, response, strict=True):
            for left, right in zip(observed, expected, strict=True):
                self.assertAlmostEqual(left, right, places=6)

    def test_bilinear_feature_includes_ordered_interaction(self) -> None:
        self.assertEqual(bilinear_features((2, 3), (5, 7)), (1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 14.0, 15.0, 21.0))

    def test_systema_single_held_target_equals_perturbed_mean_once(self) -> None:
        population = (5.0, 4.25)
        self.assertEqual(
            systema_matching_mean(["held"], {"A": (2, 3)}, population, held_single_component=True),
            population,
        )
        self.assertEqual(
            systema_matching_mean(["A", "unseen"], {"A": (3, 5)}, population),
            (4.0, 4.625),
        )

    def test_strict_reference_is_training_translation_only(self) -> None:
        self.assertEqual(strict_training_reference((8, 3), (5, 1)), (3.0, 2.0))


if __name__ == "__main__":
    unittest.main()

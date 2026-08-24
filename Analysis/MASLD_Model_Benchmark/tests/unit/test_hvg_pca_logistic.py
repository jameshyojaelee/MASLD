from __future__ import annotations

from collections import defaultdict
import unittest

from masld_bench.adapters.hvg_pca_logistic import (
    ScientificAdapterError,
    donor_class_weights,
    fold_index,
    join_hash,
)
from masld_bench.evaluators.cell_state_development import score_rows
from masld_bench.evaluators.metrics import (
    donor_class_balanced_weights as reference_weights,
    multiclass_brier_score,
    weighted_macro_f1,
)


class HvgPcaLogisticContractTests(unittest.TestCase):
    def test_fold_assignment_is_deterministic_and_unit_grouped(self) -> None:
        donors = ["GSE202379_P73", "GSE202379_P73", "GSE244832_D12"]
        first = [fold_index(item, seed=20260821, outer_folds=5) for item in donors]
        second = [fold_index(item, seed=20260821, outer_folds=5) for item in donors]
        self.assertEqual(first, second)
        self.assertEqual(first[0], first[1])
        self.assertTrue(all(0 <= item < 5 for item in first))
        with self.assertRaises(ScientificAdapterError):
            fold_index("donor", seed=1, outer_folds=1)

    def test_join_hashes_separate_rows_units_and_namespaces(self) -> None:
        row = join_hash("atlas:v1", "row", "cell-1")
        unit = join_hash("atlas:v1", "unit", "cell-1")
        other = join_hash("atlas:v2", "row", "cell-1")
        self.assertEqual(len(row), 64)
        self.assertEqual(len({row, unit, other}), 3)

    def test_training_weights_balance_classes_then_donors(self) -> None:
        donors = ["d1", "d1", "d2", "d3", "d3", "d3"]
        labels = ["a", "a", "a", "b", "b", "b"]
        weights = donor_class_weights(donors, labels)
        self.assertAlmostEqual(sum(weights), len(weights))
        mass_by_class: dict[str, float] = defaultdict(float)
        mass_by_donor_class: dict[tuple[str, str], float] = defaultdict(float)
        for donor, label, weight in zip(donors, labels, weights, strict=True):
            mass_by_class[label] += weight
            mass_by_donor_class[(donor, label)] += weight
        self.assertAlmostEqual(mass_by_class["a"], mass_by_class["b"])
        self.assertAlmostEqual(
            mass_by_donor_class[("d1", "a")],
            mass_by_donor_class[("d2", "a")],
        )

    def test_independent_evaluator_rederives_registered_metrics(self) -> None:
        observed = ["a", "a", "b", "b"]
        predicted = ["a", "b", "b", "b"]
        donors = ["d1", "d2", "d1", "d3"]
        probabilities = [
            {"a": 0.8, "b": 0.2},
            {"a": 0.4, "b": 0.6},
            {"a": 0.1, "b": 0.9},
            {"a": 0.2, "b": 0.8},
        ]
        result = score_rows(
            observed=observed,
            predicted=predicted,
            donors=donors,
            probabilities=probabilities,
            class_roster=("a", "b"),
        )
        weights = reference_weights(donors, observed)
        macro_f1, per_class = weighted_macro_f1(observed, predicted, weights)
        self.assertAlmostEqual(result["donor_class_balanced_macro_f1"], macro_f1)
        self.assertEqual(result["per_class_f1"], per_class)
        self.assertAlmostEqual(
            result["multiclass_brier_score"],
            multiclass_brier_score(observed, probabilities, weights),
        )


if __name__ == "__main__":
    unittest.main()
